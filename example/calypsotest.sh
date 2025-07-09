#!/bin/sh

#SBATCH --job-name=calypso
#SBATCH --account=commons
#SBATCH --partition=debug
#SBATCH --time=00:30:00
#SBATCH --mail-user=am503@rice.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --threads-per-core=1
#SBATCH --mem=128GB
#SBATCH --export=ALL

source /home/am503/Documents/calypso_pyenv/bin/activate

DB_DIR="/home/am503/Documents/calypso/utils/gtdb220_16S_emudb"

python /home/am503/Documents/calypso/calypso.py \
  --query_dir /home/am503/Documents/calypso/example/data \
  --ref $DB_DIR/species_taxid.fasta \
  --tax_file $DB_DIR/taxonomy.tsv \
  --outdir /home/am503/Documents/calypso/example/example_results \
  --qc 0.9 --pid 99 --raiseto 2.7 --threads 8
