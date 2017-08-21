#!/bin/bash

#SBATCH --job-name=gatk # Job name
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1440
#SBATCH --mem=40000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --reservation=workshop

begin=`date +%s`

SAMPID=$1

# modules
module load samtools
module load picard-tools
module load gatk

# KNOWN_SNPS must be in reference order
# REF must end in .fa or .fasta
REF=../ref/chr18.fa
KNOWN_SNPS=../ref/chr18.vcf

echo "marking duplicates"
# mark duplicates, assume sorted, has problems with unmapped reads with MQ!=0, so validation needs to be lenient
picard MarkDuplicates VALIDATION_STRINGENCY=LENIENT AS=true REMOVE_DUPLICATES=true I=${SAMPID}.chr18.all.bam O=${SAMPID}.chr18.markdup.bam M=${SAMPID}.chr18.metrics

# create a dictionary file for the ref
if [ ! -e ${REF%.*}.dict ]
    then
        echo "creating dictionary file for reference"
        picard CreateSequenceDictionary R=$REF O=${REF%.*}.dict
fi

echo "indexing the mark duplicates bam file"
# index the mark duplicates bam file
samtools index ${SAMPID}.chr18.markdup.bam

# need to index the reference
if [ ! -e $REF.fai ]
  then
    echo "Indexing reference"
    samtools faidx $REF
fi


echo "replacing read group and platform info"
# have to replace read groups to add platform info
picard AddOrReplaceReadGroups VALIDATION_STRINGENCY=LENIENT I=${SAMPID}.chr18.markdup.bam O=${SAMPID}.chr18.rg.bam RGID=$SAMPID RGLB=$SAMPID RGPL=illumina RGPU=$SAMPID RGSM=$SAMPID

echo "index new bam file"
# index new bam file
samtools index ${SAMPID}.chr18.rg.bam

echo "creating an intervals file for IndelRealigner"
# create an intervals file for IndelRealigner
gatk -T RealignerTargetCreator -R $REF -I ${SAMPID}.chr18.rg.bam -o ${SAMPID}.chr18.intervals

echo "running IndelRealigner"
# run IndelRealigner, can't use more than one processor
gatk -T IndelRealigner -R $REF -I ${SAMPID}.chr18.rg.bam -targetIntervals ${SAMPID}.chr18.intervals -o ${SAMPID}.chr18.realigned.bam

# had to rename chromosomes and rearrange vcf to match reference order
# cat Canis_familiaris.vcf | perl -ne 'if (/^\d+/ || /^X/ || /^M/ || /^Un/) {print "chr$_";} else {print}' > Canis_familiaris.with_chr.vcf
# perl rearrange_vcf.pl Canis_familiaris.with_chr.vcf Canis_familiaris.with_chr.rearranged.vcf

echo "running BaseRecalibrator"
# run base quality recalibrator to get table for next step
gatk -T BaseRecalibrator -I ${SAMPID}.chr18.realigned.bam -R $REF -knownSites $KNOWN_SNPS -o ${SAMPID}.chr18.recal_data.grp

echo "creating recalibrated BAM file"
# create recalibrated BAM using recalibration data from BaseRecalibrator step
gatk -T PrintReads -R $REF -I ${SAMPID}.chr18.realigned.bam -BQSR ${SAMPID}.chr18.recal_data.grp -o ${SAMPID}.chr18.recalibrated.bam

echo "running haplotypecaller"
gatk -T HaplotypeCaller -R $REF -ERC GVCF -I ${SAMPID}.chr18.recalibrated.bam -o ${SAMPID}.chr18.g.vcf

end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: $elapsed

