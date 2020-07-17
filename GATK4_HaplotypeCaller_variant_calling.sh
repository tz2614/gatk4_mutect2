#!/bin/bash

# Use GATK4 version 4.0.7.0 to call bam files in /users/tz1/sandbox/gatk4

reference_fa="/mnt/lustre/references/hg19/hg19_validated.fa"

# test reference fasta
#reference_fa="/users/tz1/git/GATK4/refs/ucsc.hg19.HODS.fasta"

germline_vars="/users/tz1/git/GATK4/refs/af-only-gnomad.raw.sites.hg19.vcf.gz"
gatk="gatk"
vcf_format=".snps.indels.haplotypeCaller.vcf"
log_format=".log"

data_dir="/users/tz1/git/GATK4/MYH7/"

# test data directory
#data_dir="/users/tz1/git/GATK4/test_data/"
#interval_list="/users/tz1/git/GATK4/refs/hotspot_region_Ion_AmpliSeq_CHPv2_Somatic_Hotspots_GATK4_formatted3_output.interval_list"
interval_list="/mnt/lustre/references/targetedReseq/Cardiac/ALL/current/ALL.inhouse.bed"

#cmd="$gatk Mutect2 -R $reference_fa -L $interval_list -XL chrM -I $bam -O $vcf"

for bam in $(find $data_dir -name "*.bam")
    do
        #vcf=${bam/".bam"/$vcf_format}
        unfiltered_vcf=${bam/".bam"/$vcf_format}
        unfiltered_log=${bam/".bam"/$log_format}
        echo "variant calling $bam using GATK4 haplotypecaller"
        $gatk HaplotypeCaller -R $reference_fa -L $interval_list -I $bam -O $unfiltered_vcf | tee $unfiltered_log
        #$gatk Mutect2 -R $reference_fa ---germline-resource $germline_vars -L $interval_list -I $bam -O $unfiltered_vcf | tee $unfiltered_log &&
    done

