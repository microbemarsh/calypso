#!/bin/sh

#SBATCH --job-name=calypso
#SBATCH --account=commons
#SBATCH --partition=debug
#SBATCH --time=00:30:00
#SBATCH --mail-user=email@email.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --threads-per-core=1
#SBATCH --mem=128GB
#SBATCH --export=ALL

source /path/to/calypso_pyenv/bin/activate

DB_DIR="/path/to/calypso/utils/gtdb220_16S_emudb"

python /path/to/calypso/calypso.py \
  --query_dir /path/to/calypso/example/data \
  --ref $DB_DIR/species_taxid.fasta \
  --tax_file $DB_DIR/taxonomy.tsv \
  --outdir /path/to/calypso/example/example_results \
  --qc 0.9 --pid 99 --raiseto 2.7 --threads 8
