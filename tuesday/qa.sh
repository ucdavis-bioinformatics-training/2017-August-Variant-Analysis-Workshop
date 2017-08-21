#!/bin/bash

#SBATCH --job-name=trim # Job name
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=180
#SBATCH --mem=2000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --reservation=workshop


begin=`date +%s`

sample=$1
R1=${sample}.chr18.R1.fastq
R2=${sample}.chr18.R2.fastq

echo $sample

module load scythe
module load sickle

scythe -a adapters.fasta -q sanger -o ${sample}.chr18.R1.scythe.fastq $R1
scythe -a adapters.fasta -q sanger -o ${sample}.chr18.R2.scythe.fastq $R2
sickle pe -f ${sample}.chr18.R1.scythe.fastq -r ${sample}.chr18.R2.scythe.fastq -t sanger -o ${sample}.chr18.R1.sickle.fastq -p ${sample}.chr18.R2.sickle.fastq -s ${sample}.chr18.singles.fastq

end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: $elapsed

