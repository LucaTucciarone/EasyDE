# =============================================================================
# EasyDE — Snakefile
# =============================================================================

import yaml
import csv
import os
import re


# -----------------------------------------------------------------------------
# Config and helpers
# -----------------------------------------------------------------------------

PIPELINE_CONFIG = config.get("pipeline_config", "config/config.yaml")

with open(PIPELINE_CONFIG) as f:
    cfg = yaml.safe_load(f)

RESULTS_DIR = cfg["outputs"]["results_dir"]
LOGS_DIR    = cfg["outputs"]["logs_dir"]
SCRIPTS_DIR = "workflow/scripts"


def sanitize(name):
    """Mirror R sanitize_names(): replace non-alphanumeric chars with underscores."""
    return re.sub(r"[^A-Za-z0-9]", "_", name)


def parse_contrasts(contrasts_file):
    """
    Parse contrasts CSV.
    Returns:
        contrasts: { contrast_id: [ original_stratum, ... ] }
        safe_map:  { (contrast_id, stratum_safe): original_stratum }
    """
    contrasts = {}
    safe_map  = {}
    with open(contrasts_file, newline="") as f:
        lines = [l for l in f if not l.strip().startswith("#") and l.strip()]
    reader = csv.DictReader(lines)
    for row in reader:
        cid    = row["contrast_id"].strip()
        strata = [s.strip() for s in row["strata"].split("|") if s.strip()]
        contrasts[cid] = strata
        for s in strata:
            safe_map[(cid, sanitize(s))] = s
    return contrasts, safe_map


CONTRASTS, SAFE_MAP = parse_contrasts(cfg["inputs"]["contrasts_file"])

# Use sanitized stratum names as Snakemake wildcards (safe for file paths)
ALL_CONTRAST_STRATUM = [
    (cid, sanitize(s))
    for cid, strata in CONTRASTS.items()
    for s in strata
]


def original_stratum(wildcards):
    """Recover the original stratum name from sanitized wildcard."""
    return SAFE_MAP[(wildcards.contrast, wildcards.stratum_safe)]


def inter(contrast, stratum_safe, filename):
    return os.path.join(RESULTS_DIR, contrast, stratum_safe, "intermediates", filename)


def final(contrast, stratum_safe, filename):
    return os.path.join(RESULTS_DIR, contrast, stratum_safe, "finals", filename)


# -----------------------------------------------------------------------------
# Target rule
# -----------------------------------------------------------------------------

rule all:
    input:
        expand(
            os.path.join(RESULTS_DIR, "{contrast}", "summary", "contrast_summary.csv"),
            contrast=list(CONTRASTS.keys())
        )


# -----------------------------------------------------------------------------
# Rule 01 — Validate inputs (once per contrast)
# -----------------------------------------------------------------------------

rule validate:
    input:
        config    = PIPELINE_CONFIG,
        meta      = cfg["inputs"]["sample_metadata"],
        contrasts = cfg["inputs"]["contrasts_file"],
    output:
        touch(os.path.join(LOGS_DIR, "01_validate.done"))
    log:
        os.path.join(LOGS_DIR, "01_validate.log")
    shell:
        """
        Rscript {SCRIPTS_DIR}/01_prepare_inputs.R \
            --config {input.config} \
            > {log} 2>&1
        """


# -----------------------------------------------------------------------------
# Rule 02 — Prepare coldata
# -----------------------------------------------------------------------------

rule prepare_coldata:
    input:
        done      = os.path.join(LOGS_DIR, "01_validate.done"),
        config    = PIPELINE_CONFIG,
        meta      = cfg["inputs"]["sample_metadata"],
        contrasts = cfg["inputs"]["contrasts_file"],
    output:
        coldata = inter("{contrast}", "{stratum_safe}", "coldata.csv"),
        counts  = inter("{contrast}", "{stratum_safe}", "counts_filtered.csv"),
        namemap = inter("{contrast}", "{stratum_safe}", "name_map.csv"),
    params:
        stratum = original_stratum
    log:
        os.path.join(LOGS_DIR, "{contrast}", "{stratum_safe}_02.log")
    shell:
        """
        Rscript {SCRIPTS_DIR}/02_prepare_coldata.R \
            --config {input.config} \
            --contrast {wildcards.contrast} \
            --stratum "{params.stratum}" \
            > {log} 2>&1
        """


# -----------------------------------------------------------------------------
# Rule 03 — Base DESeq2
# -----------------------------------------------------------------------------

rule deseq_base:
    input:
        config  = PIPELINE_CONFIG,
        coldata = inter("{contrast}", "{stratum_safe}", "coldata.csv"),
        counts  = inter("{contrast}", "{stratum_safe}", "counts_filtered.csv"),
    output:
        results = final("{contrast}", "{stratum_safe}", "results_base.tsv"),
    params:
        stratum = original_stratum
    log:
        os.path.join(LOGS_DIR, "{contrast}", "{stratum_safe}_03.log")
    shell:
        """
        Rscript {SCRIPTS_DIR}/03_run_deseq_base.R \
            --config {input.config} \
            --contrast {wildcards.contrast} \
            --stratum "{params.stratum}" \
            > {log} 2>&1
        """


# -----------------------------------------------------------------------------
# Rule 04 — RUVseq
# -----------------------------------------------------------------------------

rule ruvseq:
    input:
        config  = PIPELINE_CONFIG,
        counts  = inter("{contrast}", "{stratum_safe}", "counts_filtered.csv"),
        results = final("{contrast}", "{stratum_safe}", "results_base.tsv"),
    output:
        coldata = inter("{contrast}", "{stratum_safe}", "ruvseq_best_coldata.csv"),
        summary = inter("{contrast}", "{stratum_safe}", "ruvseq_summary.csv"),
    params:
        stratum = original_stratum
    log:
        os.path.join(LOGS_DIR, "{contrast}", "{stratum_safe}_04.log")
    shell:
        """
        Rscript {SCRIPTS_DIR}/04_run_ruvseq.R \
            --config {input.config} \
            --contrast {wildcards.contrast} \
            --stratum "{params.stratum}" \
            > {log} 2>&1
        """


# -----------------------------------------------------------------------------
# Rule 05 — Final DESeq2
# -----------------------------------------------------------------------------

rule deseq_final:
    input:
        config  = PIPELINE_CONFIG,
        counts  = inter("{contrast}", "{stratum_safe}", "counts_filtered.csv"),
        coldata = inter("{contrast}", "{stratum_safe}", "ruvseq_best_coldata.csv"),
        summary = inter("{contrast}", "{stratum_safe}", "ruvseq_summary.csv"),
    output:
        results = final("{contrast}", "{stratum_safe}", "results_ruv.tsv"),
    params:
        stratum = original_stratum
    log:
        os.path.join(LOGS_DIR, "{contrast}", "{stratum_safe}_05.log")
    shell:
        """
        Rscript {SCRIPTS_DIR}/05_run_deseq_final.R \
            --config {input.config} \
            --contrast {wildcards.contrast} \
            --stratum "{params.stratum}" \
            > {log} 2>&1
        """


# -----------------------------------------------------------------------------
# Rule 06 — fGSEA
# -----------------------------------------------------------------------------

rule fgsea:
    input:
        config       = PIPELINE_CONFIG,
        results_base = final("{contrast}", "{stratum_safe}", "results_base.tsv"),
        results_ruv  = final("{contrast}", "{stratum_safe}", "results_ruv.tsv"),
        gene_sets    = cfg["fgsea"]["gene_sets_file"],
    output:
        touch(inter("{contrast}", "{stratum_safe}", "fgsea.done"))
    params:
        stratum = original_stratum
    log:
        os.path.join(LOGS_DIR, "{contrast}", "{stratum_safe}_06.log")
    shell:
        """
        Rscript {SCRIPTS_DIR}/06_run_fgsea.R \
            --config {input.config} \
            --contrast {wildcards.contrast} \
            --stratum "{params.stratum}" \
            > {log} 2>&1
        """


# -----------------------------------------------------------------------------
# Rule 07 — Aggregate results (once per contrast)
# -----------------------------------------------------------------------------

rule aggregate:
    input:
        config = PIPELINE_CONFIG,
        done   = lambda wc: [
            inter(wc.contrast, sanitize(s), "fgsea.done")
            for s in CONTRASTS[wc.contrast]
        ]
    output:
        summary = os.path.join(
            RESULTS_DIR, "{contrast}", "summary", "contrast_summary.csv"
        )
    log:
        os.path.join(LOGS_DIR, "{contrast}", "07_aggregate.log")
    shell:
        """
        Rscript {SCRIPTS_DIR}/07_aggregate_results.R \
            --config {input.config} \
            --contrast {wildcards.contrast} \
            > {log} 2>&1
        """
