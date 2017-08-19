#!/bin/bash

#SBATCH --job-name=freebayes # Job name
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=600
#SBATCH --mem=40000 # Memory pool for all cores (see also --mem-per-cpu)
#    SBATCH --reservation=workshop


begin=`date +%s`

BAMLIST=$1
REF=../ref/chr18.fa

module load freebayes

freebayes -L $BAMLIST -f $REF --vcf freebayes.chr18.all.vcf

end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: $elapsed

