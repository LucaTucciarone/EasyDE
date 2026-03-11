# Installation

## 1. Install micromamba

```bash
# ========= Install Micromamba (user-level) =========
wget -O ~/micromamba.tar.bz2 https://micro.mamba.pm/api/micromamba/linux-64/latest
tar -xvjf ~/micromamba.tar.bz2 -C ~/
rm ~/micromamba.tar.bz2

# ========= Initialize shell integration =========
~/bin/micromamba shell init -s bash -r ~/micromamba
source ~/.bashrc

  

# ========= Verify installation =========
which micromamba
micromamba info

# ========= Configure global settings =========
micromamba config append channels bioconda
micromamba config append channels conda-forge
micromamba config append channels defaults
micromamba config set channel_priority strict
micromamba config set always_yes true
micromamba config set auto_activate_base false
micromamba config set repodata_fallback true

# ========= Verify config =========

micromamba config list

```

## 2. Create the environment

```bash
cd /path/to/EasyDE
micromamba env create -f installation/EasyDE_install.yml
```

This installs R, DESeq2, RUVSeq, fgsea, Snakemake, and all dependencies.
Takes ~5-10 minutes on first run.

## 3. Activate

```bash
micromamba activate EasyDE
```

Activate this environment every time before running any pipeline command.

## 4. Verify

```bash
Rscript -e "cat('R ok\n')"               # should print "R ok" and exit 0
snakemake --version                        # should print >=7.x
```

> **Using conda/mamba?** Replace `micromamba` with `conda` or `mamba` — the
> syntax is identical.

<details>
<summary><strong>Known issue: Conda R exit code 255</strong></summary>

Some conda-installed R builds return exit code 255 on completion regardless
of success. Test with `Rscript -e "cat('hello\n')"` — the exit code should
be 0. If not, either reinstall R (`mamba install -c conda-forge r-base
--force-reinstall`) or point the Snakefile to your system R.

</details>
