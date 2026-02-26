# DESeq_RUVseq


pipeline/
├── Snakefile                        ← orchestrator (the DAG)
├── config/
│   └── config.yaml                  ← all hardcoded values move here
├── contrasts/
│   └── pankbase_new_contrasts.csv   ← your contrast definitions
├── workflow/
│   ├── rules/                       ← one .smk file per logical step
│   │   ├── fetch_metadata.smk
│   │   ├── preprocess.smk
│   │   ├── ruvseq.smk
│   │   ├── deseq2.smk
│   │   └── fgsea.smk
│   └── scripts/                     ← modular R scripts (no more monolith)
│       ├── utils/
│       │   └── clean_utils.R
│       ├── preprocess.R
│       ├── run_ruvseq.R
│       ├── run_deseq2.R
│       └── run_fgsea.R
├── resources/                       ← metadata lives here, downloaded once
│   └── .gitkeep
├── logs/                            ← structured per-rule logs
└── results/                        ← all outputs, organized by contrast
    └── {contrast_id}/
        ├── intermediates/
        ├── summaries/
        └── plots/





workflow/scripts/
├── utils/
│   ├── io_utils.R          ← metadata loading, URL building, file I/O
│   ├── filter_utils.R      ← gene + sample filtering (from 0_clean_utils.R)
│   ├── plot_utils.R        ← all plotting functions
│   └── stats_utils.R       ← labs.function, correlation helpers
│
├── 01_fetch_metadata.R     ← downloads & caches donor/biosample metadata
├── 02_prepare_coldata.R    ← per cell-type sample filtering & metadata merge
├── 03_run_deseq_base.R     ← initial DESeq2
├── 04_run_ruvseq.R         ← RUVseq + k selection
├── 05_run_deseq_final.R    ← final DESeq2 with best RUV formula
├── 06_run_fgsea.R          ← fGSEA + pathway plots
└── 07_aggregate_stats.R    ← collects de_stats across cell types





One count file per celltype
"|" as a separator - or what do u suggest? I would like to use "a saperator" but definitely the _ is among the worst choices
As for the togglable: RuvSeq, fGSEA, Paired analysis (subsets for only donors that are in both conditions we are comparing - something that I had to manually do in the csv file I gave you). Certain actions are kind of dependencies of choices, for example: if you choose to run RuvSeq I would just output all the relevant information (correlation plots included) that help understand what happened.
So overall I would make the user decide: Base deseq2 + ruvseq + FGSEA (on either) and we can discuss and reason what are the best summarizations of the results AND OF THE MESSAGE/ERRORS - which is one of the WORST things about this pipelines, it is unbearbly hard figure out what job failed because of what
