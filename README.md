# calypso
Reimplementation of [ATLAS](https://github.com/marbl/ATLAS) designed for use with [minimap2](https://github.com/lh3/minimap2)
and [aln_stats](https://github.com/microbemarsh/aln_stats).

## Dependencies

* [minimap2](https://github.com/lh3/minimap2)
* [samtools](https://www.htslib.org/)
* [numpy](https://numpy.org/)
* [pysam](http://pysam.readthedocs.io/en/latest/)
and some more but these are all installed for you within the method below. 

## Installation

Please make sure to have a working conda installation first. If you don't, please select the appropriate version found [here](https://github.com/conda-forge/miniforge).

```
# Clone the respository
git clone https://github.com/microbemarsh/calypso.git

# Go to repo
cd calypso

# Create conda environment
conda env create -n calypso -f calypso_env.yml

# Activate conda environment
conda activate calypso
```

## Usage

```
# Activate the conda environment
conda activate calypso

# And run the python script
python calypso.py \
    --query queries.fasta \
    --ref ref_16S.fasta \
    --tax_file tax_map.tsv \
    --outdir results
```

## Why the name?
Calypso is the daughter of Atlas in greek mythology.
