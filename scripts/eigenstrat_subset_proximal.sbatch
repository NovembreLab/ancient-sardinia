#!/bin/sh

#SBATCH --job-name=subset
#SBATCH --output=log/subset.%J.out
#SBATCH --error=log/subset.%J.err
#SBATCH --time=02:00:00
#SBATCH --partition=jnovembre
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10g

module load gsl
module load openblas

CONVERTF='/project/jnovembre/jhmarcus/ancient-sardinia/bin/convertf'
PAR=$(mktemp merge_par.XXXXXX)

PROJECT='/project/jnovembre/jhmarcus/ancient-sardinia'
PREFIX="${PROJECT}/output/eigenstrat/full/qpall"
OUT="${PROJECT}/output/eigenstrat/subset/qpall_proximal_subset"

echo "genotypename: $PREFIX.geno" > $PAR
echo "snpname: $PREFIX.snp" >> $PAR
echo "indivname: ${PREFIX}_proximal_uniqID.ind" >> $PAR
echo "outputformat: PACKEDANCESTRYMAP" >> $PAR
echo "genotypeoutname: $OUT.geno" >> $PAR
echo "snpoutname: $OUT.snp" >> $PAR
echo "indivoutname: $OUT.ind" >> $PAR
echo "poplistname: ${PROJECT}/output/eigenstrat/subset/proximal_populations.txt" >> $PAR
echo "hashcheck: NO" >> $PAR

$CONVERTF -p $PAR
rm $PAR
echo "Subset Complete!"
