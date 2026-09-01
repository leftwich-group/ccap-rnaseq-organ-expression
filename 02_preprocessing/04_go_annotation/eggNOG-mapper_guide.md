# eggNOG-mapper setup and run guide

This is the primary GO annotation source used in the downstream R analysis
(`analysis/Ccap_Bulk_Transcriptomics.Rmd`). It documents how eggNOG-mapper
was installed and run on the HPC cluster to annotate the short-read and
long-read proteomes with GO terms.

Throughout, `$PROJECT_ROOT` refers to your own working directory on the
cluster (e.g. `/path/to/your/project`).

## 1. Install eggNOG-mapper (conda/mamba route)

```bash
module load python/anaconda/2023.07/3.11.4

# install mamba
conda install -n base -c conda-forge mamba

# create environment for eggNOG-mapper
mamba create -n eggnog_env python=3.10

# activate the environment
source activate eggnog_env

# install dependencies and eggNOG-mapper
mamba install -c bioconda -c conda-forge eggnog-mapper

# verify installation
emapper.py -h
emapper.py --version
diamond --version
```

### Download the eggNOG reference database

```bash
mkdir -p ${PROJECT_ROOT}/GO_Terms/eggNOG-mapper/eggnog_data
download_eggnog_database.py --data_dir ${PROJECT_ROOT}/GO_Terms/eggNOG-mapper/eggnog_data
```

On this cluster the conda/mamba route ran into dependency-solver issues
partway through, so the install was redone with micromamba instead (see
below). If `conda`/`mamba` works cleanly for you, you can skip straight to
running eggNOG-mapper once the database has downloaded.

## 2. Install eggNOG-mapper (micromamba route)

```bash
cd ~
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba

# add micromamba to path
export PATH=~/bin:$PATH

# initialise shell
micromamba shell init -s bash
source ~/.bashrc

# set root prefix
export MAMBA_ROOT_PREFIX=~/micromamba
echo 'export MAMBA_ROOT_PREFIX=~/micromamba' >> ~/.bashrc
source ~/.bashrc

# create environment
micromamba create -n eggnog_env -c conda-forge -c bioconda eggnog-mapper python=3.10

# activate environment
micromamba activate eggnog_env

# verify installation
emapper.py -h
emapper.py --version
diamond --version
```

At this point `emapper.py` will report that its reference data hasn't been
downloaded yet.

### Download the eggNOG reference database

```bash
mkdir -p ${PROJECT_ROOT}/GO_Terms/eggNOG-mapper/eggnog_data
download_eggnog_data.py --data_dir ${PROJECT_ROOT}/GO_Terms/eggNOG-mapper/eggnog_data
```

The automatic downloader is known to fail intermittently for this database;
if it does, download the files manually instead:

```bash
DATA_DIR=${PROJECT_ROOT}/GO_Terms/eggNOG-mapper/eggnog_data
mkdir -p "$DATA_DIR"
cd "$DATA_DIR"

wget http://eggnog6.embl.de/download/emapperdb-5.0.2/eggnog.db.gz
wget http://eggnog6.embl.de/download/emapperdb-5.0.2/eggnog.taxa.tar.gz
wget http://eggnog6.embl.de/download/emapperdb-5.0.2/eggnog_proteins.dmnd.gz

gunzip eggnog.db.gz
tar -xzf eggnog.taxa.tar.gz
rm eggnog.taxa.tar.gz
gunzip eggnog_proteins.dmnd.gz
```

Verify:

```bash
emapper.py --data_dir ${PROJECT_ROOT}/GO_Terms/eggNOG-mapper/eggnog_data --version
```

## 3. Run eggNOG-mapper

```bash
sbatch ${PROJECT_ROOT}/GO_Terms/eggNOG-mapper/scripts/run_eggnog_long_short.sh
```

**Note:** if eggNOG-mapper was installed via micromamba (as above), the job
log will contain a spurious "data not found" error early on that can be
ignored — the run proceeds normally regardless.

The `run_eggnog_long_short.sh` sbatch script referenced above runs
`emapper.py` against both the short-read and long-read proteomes and is not
included in this repository — reconstruct it by pointing `emapper.py` at
each proteome FASTA with `--data_dir` set to the downloaded database above,
or consult the [eggNOG-mapper
documentation](https://github.com/eggnogdb/eggnog-mapper) for the standard
invocation. The script is also available on request from the repository
maintainer.
