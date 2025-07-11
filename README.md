# calypso
Reimplementation of [ATLAS](https://github.com/marbl/ATLAS) designed for use with [parasail](https://github.com/jeffdaily/parasail)
instead of [BLAST](https://www.ncbi.nlm.nih.gov/books/NBK279690/).

Calypso aligns each query sequence to the reference database (16S sequences from [GTDB 220](https://gtdb.ecogenomic.org/stats/r220)), filters out low-quality hits,
uses an entropy-based MSA test to weed out stray alignments, groups hits by co-occurrence, assigns reads to graph clusters,
calls an LCA on both the outliers and clusters, and finally picks the deepest taxonomic level supported by the best percent identity.
Then it merges everything into a single report.

The database was created in a similar fashion to how [PICRUSt2](https://github.com/picrust/picrust2) makes it's database, however,
we choose the longest 16S sequence when multiple 16S copies are [barrnap](https://github.com/tseemann/barrnap)'d from a singular taxid. 
This leads to one 16S sequence being from one genomic reference from GTDB r220, making it possible to use GTDB taxonomies, without making consensus sequences.

## Installation
Please make sure to have an active python base environment to which you can invoke [venv](https://docs.python.org/3/library/venv.html).

```
# Clone the respository
git clone https://github.com/microbemarsh/calypso.git

# Make environment
python -m venv calypso_env

# Activate environment
source ~/calypso_env/bin/activate

# Install dependencies
pip install numpy scipy biopython networkx python-louvain parasail

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
    [--qc 0.9] [--pid 99] [--raiseto 2.7]
```

## Why the name?
Calypso is the daughter of Atlas in greek mythology.
