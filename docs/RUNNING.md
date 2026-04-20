# Running the Pipeline

## Two modes of use

EasyDE can be run in two ways depending on your setup:

| Mode | When to use |
|------|-------------|
| **In-place** | Development, testing, or a single study — run from the pipeline directory itself |
| **External working directory** | Production / multi-study deployments — pipeline code stays in one place, each study or tissue runs from its own directory |

---

## In-place runs (pipeline directory as working directory)

The simplest setup. `cd` into the EasyDE directory and run Snakemake directly.
Snakemake state (`.snakemake/`) and any relative output paths land inside the
pipeline directory.

### Run all contrasts

```bash
cd /path/to/EasyDE
micromamba activate EasyDE

snakemake --profile profiles/local \
    --config pipeline_config=config/config_traits.yaml
```

This runs **all** contrasts defined in your contrasts CSV across **all**
their strata. For 20 contrasts × 11 cell types, that is ~1,100 parallel jobs.

### Run one contrast

```bash
snakemake --profile profiles/local \
    --config pipeline_config=config/config_traits.yaml \
    results/Diabetic_vs_ND/summary/contrast_summary.csv
```

### Dry run (see what would execute)

```bash
snakemake --profile profiles/local \
    --config pipeline_config=config/config_traits.yaml -n
```

### Run treatment (paired) analysis

```bash
snakemake --profile profiles/local \
    --config pipeline_config=config/config_treatments.yaml
```

In paired mode the pipeline flow is **01 → 02 → 03 → 06 → 07 → 08 → 09** (steps 04
and 05 are automatically skipped — see [Methods](METHODS.md#paired-analysis)).
Benchmarking (step 07) is optional and only runs when `benchmarking.signature_file`
is configured.

### Clean and re-run

```bash
# Remove results for one contrast
rm -rf results/Diabetic_vs_ND/

# Remove all results
rm -rf results/ logs/

# Then re-run
snakemake --profile profiles/local \
    --config pipeline_config=config/config_traits.yaml
```

---

## External working directory (recommended for multi-study deployments)

When running multiple studies or tissues in parallel from a shared pipeline
installation, each run should have its own working directory. Snakemake writes
its state to `.snakemake/` relative to the current directory — so `cd`-ing
into a per-study directory before invoking Snakemake is all that is needed to
keep runs fully isolated. The pipeline code itself stays in one place and is
never written to during a run — `Snakefile`, `workflow/`, and `resources/` are
read-only from each run's perspective.

```
/path/to/pipeline/EasyDE/        ← shared, never touched during runs
/path/to/runs/study_A/           ← cd here to run study A
    .snakemake/                  ← isolated state, created automatically
    results/
    logs/
    snakemake.log
/path/to/runs/study_B/           ← cd here to run study B
    .snakemake/
    ...
```

### Run from an external directory

```bash
PIPELINE=/path/to/EasyDE
WORKDIR=/path/to/runs/study_A
CONFIG=/path/to/runs/study_A/input/config.yaml

cd "$WORKDIR"
micromamba activate EasyDE

snakemake --snakefile  "${PIPELINE}/Snakefile" \
          --profile    "${PIPELINE}/profiles/local" \
          --config     pipeline_config="${CONFIG}" \
          --cores      10
```

Key points:
- `--snakefile` points at the pipeline's `Snakefile` by absolute path
- `--profile` points at the pipeline's profiles directory by absolute path
- `pipeline_config` must be an absolute path (or relative to `$WORKDIR`)
- `.snakemake/`, `results/`, and `logs/` will be created inside `$WORKDIR`
- **Paths inside your config YAML** (counts dir, metadata file, contrasts CSV, GMT files) resolve relative to `$WORKDIR`, not the pipeline dir. Use absolute paths in the config, or stage all inputs under `$WORKDIR/input/`
- **`resources/`** (GMT files, signatures) lives inside the pipeline dir. Reference it with the absolute path `${PIPELINE}/resources/...` in your config, or symlink it into `$WORKDIR`

### Running multiple studies in parallel (one screen per study)

This pattern is used for production multi-tissue runs. Each study gets its own
`screen` session and its own working directory; Snakemake state never collides.

```bash
PIPELINE=/path/to/EasyDE
RUNS=/path/to/runs
MAMBA=/path/to/micromamba
MAMBA_ROOT=/path/to/micromamba_root
CORES=10

for study in study_A study_B study_C; do
    workdir="${RUNS}/${study}"
    config="${workdir}/input/config.yaml"
    logfile="${workdir}/snakemake.log"

    mkdir -p "$workdir"

    screen -dmS "easyde_${study}" bash -c "
export MAMBA_EXE='${MAMBA}'
export MAMBA_ROOT_PREFIX='${MAMBA_ROOT}'
eval \"\$(${MAMBA} shell hook --shell bash --root-prefix \${MAMBA_ROOT} 2>/dev/null)\"
micromamba activate EasyDE
cd ${workdir}
echo '=== Starting ${study} at '\$(date)' ===' > ${logfile}
snakemake --snakefile  ${PIPELINE}/Snakefile \
          --profile    ${PIPELINE}/profiles/local \
          --config     pipeline_config=${config} \
          --cores      ${CORES} \
          2>&1 | tee -a ${logfile}
echo '=== Done: ${study} at '\$(date)' ===' >> ${logfile}
"
    echo "Launched: easyde_${study}"
done
```

### Rerun incomplete jobs

```bash
cd "$WORKDIR"
snakemake --snakefile "${PIPELINE}/Snakefile" \
          --profile   "${PIPELINE}/profiles/local" \
          --config    pipeline_config="${CONFIG}" \
          --cores     10 \
          --rerun-incomplete
```

---

## Running on SLURM

For cluster execution, use the SLURM profile. Edit
`profiles/slurm/config.yaml` to set your account, partition, and QOS:

```yaml
# profiles/slurm/config.yaml (edit these)
default-resources:
  slurm_account: your_account
  slurm_partition: your_partition
  slurm_qos: your_qos
  mem_mb: 8000
  time: 240      # minutes
  cpus: 1
```

Then run (use absolute paths — SLURM jobs may not inherit `$PWD`):

```bash
PIPELINE=/path/to/EasyDE

snakemake --snakefile "${PIPELINE}/Snakefile" \
          --profile   "${PIPELINE}/profiles/slurm" \
          --config    pipeline_config="${PIPELINE}/config/config_traits.yaml"
```

> **Note:** SLURM is only available on HPC clusters where `sbatch` is
> installed. On standalone servers (no job scheduler), always use
> `profiles/local`.

---

## Running Step by Step

For debugging or educational purposes, you can run each script manually.
All commands assume you are in the pipeline root directory.

```bash
micromamba activate EasyDE

CONFIG="config/config_traits.yaml"
CONTRAST="Diabetic_vs_ND"
STRATUM="Beta"

# Step 01 — Validate (once per run)
Rscript workflow/scripts/01_prepare_inputs.R --config "$CONFIG"

# Step 02 — Prepare coldata
Rscript workflow/scripts/02_prepare_coldata.R \
    --config "$CONFIG" --contrast "$CONTRAST" --stratum "$STRATUM"

# Step 03 — Base DESeq2
Rscript workflow/scripts/03_run_deseq_base.R \
    --config "$CONFIG" --contrast "$CONTRAST" --stratum "$STRATUM"

# Step 04 — RUVseq (skipped automatically in paired mode)
Rscript workflow/scripts/04_run_ruvseq.R \
    --config "$CONFIG" --contrast "$CONTRAST" --stratum "$STRATUM"

# Step 05 — Final DESeq2 with W factors (skipped automatically in paired mode)
Rscript workflow/scripts/05_run_deseq_final.R \
    --config "$CONFIG" --contrast "$CONTRAST" --stratum "$STRATUM"

# Step 06 — fGSEA
Rscript workflow/scripts/06_run_fgsea.R \
    --config "$CONFIG" --contrast "$CONTRAST" --stratum "$STRATUM"

# Step 07 — Benchmark signatures (optional — requires benchmarking.signature_file in config)
Rscript workflow/scripts/07_benchmark_signatures.R \
    --config "$CONFIG" --contrast "$CONTRAST" --stratum "$STRATUM"

# Step 08 — Aggregate (once per contrast, no --stratum needed)
Rscript workflow/scripts/08_aggregate_results.R \
    --config "$CONFIG" --contrast "$CONTRAST"

# Step 09 — Pipeline summary (once after ALL contrasts complete)
Rscript workflow/scripts/09_pipeline_summary.R --config "$CONFIG"
```
