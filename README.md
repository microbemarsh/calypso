# calypso
Reimplementation of [ATLAS](https://github.com/marbl/ATLAS) designed for use with [minimap2](https://github.com/lh3/minimap2)
instead of [BLAST](https://www.ncbi.nlm.nih.gov/books/NBK279690/)

## Installation
Please make sure to have an active python base environment to which you can invoke venv.

and

[minimap2](https://github.com/lh3/minimap2)

```
# Clone the respository
git clone https://github.com/microbemarsh/calypso.git

# Make environment
python -m venv calypso_env

# Activate environment
source ~/calypso_env/bin/activate

# Install dependencies
pip install mappy numpy scipy biopython networkx python-louvain pysam

# Test if you minimap2 in your path (if not repo is above)
minimap2 -h

# And test it out
python calypso.py -h
```

## Usage

```
# Activate the conda environment
conda activate calypso

# And run the python script
python calypso.py \
    --query_dir fasta_dir \
    --ref ref_16S.fasta \
    --tax_file tax_map.tsv \
    --outdir results \
    [--qc 0.9] [--pid 99] [--raiseto 2.7] [--type map-ont] [--threads 8]
```

## Why the name?
Calypso is the daughter of Atlas in greek mythology.
