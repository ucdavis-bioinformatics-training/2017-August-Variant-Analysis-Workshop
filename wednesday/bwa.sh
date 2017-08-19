#!/bin/bash

#SBATCH --job-name=align # Job name
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=600
#SBATCH --mem=20000 # Memory pool for all cores (see also --mem-per-cpu)
#    SBATCH --reservation=workshop


begin=`date +%s`

sample=$1
R1=${sample}.chr18.R1.sickle.fastq
R2=${sample}.chr18.R2.sickle.fastq
SINGLES=${sample}.chr18.singles.fastq
REF=../ref/chr18.fa

echo $sample

module load bwa
module load samtools

bwa mem -t 4 -R "@RG\tID:$sample\tSM:$sample" $REF $R1 $R2 > ${sample}.chr18.paired.sam
bwa mem -t 4 -R "@RG\tID:$sample\tSM:$sample" $REF $SINGLES > ${sample}.chr18.singles.sam

samtools view -@ 4 -b ${sample}.chr18.paired.sam > ${sample}.chr18.paired.bam
samtools view -@ 4 -b ${sample}.chr18.singles.sam > ${sample}.chr18.singles.bam

samtools sort -@ 4 -o ${sample}.chr18.paired.sorted.bam ${sample}.chr18.paired.bam
samtools sort -@ 4 -o ${sample}.chr18.singles.sorted.bam ${sample}.chr18.singles.bam

samtools merge -@ 4 ${sample}.chr18.all.bam ${sample}.chr18.paired.sorted.bam ${sample}.chr18.singles.sorted.bam
samtools index ${sample}.chr18.all.bam

end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: $elapsed

