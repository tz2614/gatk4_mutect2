#!/bin/bash

# Use GATK4 vrsion 4.1.4.1 to call bam files in /users/tz1/sandbox/gatk4

reference_fa="/mnt/lustre/references/hg19/hg19_validated.fa"
germline_vars="/users/tz1/git/GATK4_Mutect2/refs/af-only-gnomad.raw.sites.hg19.vcf.gz"
gatk="/users/tz1/apps/GATK4/gatk-4.1.4.1/gatk"
#vcf_format=".vcf.gz"
#log_format=".log"
data_dir="/users/tz1/sandbox/gatk4/"
#test_data_dir="/users/tz1/GATK4_Mutect2/test_data/"
#interval_list="/users/tz1/GATK4_Mutect2/refs/hotspot_region_Ion_AmpliSeq_CHPv2_Somatic_Hotspots_GATK4_formatted3_output.interval_list"
interval_list="/mnt/lustre/references/targetedReseq/Cardiac/ALL/current/ALL.inhouse.bed"

#cmd="$gatk Mutect2 -R $reference_fa -L $interval_list -XL chrM -I $bam -O $vcf"

for bam in $(find $data_dir -name "*sorted.nodup.bam")
    do
        #vcf=${bam/".bam"/$vcf_format}
        unfiltered_vcf=${bam/".sorted.nodup.bam"/".snps.indels.haplotypeCaller.vcf"}
        unfiltered_log=${bam/".sorted.nodup.bam"/.log}
        echo "variant calling $bam using GATK4 haplotypecaller"
        $gatk HaplotypeCaller -R $reference_fa -L $interval_list -I $bam -O $unfiltered_vcf | tee $unfiltered_log
        #$gatk Mutect2 -R $reference_fa ---germline-resource $germline_vars -L $interval_list -I $bam -O $unfiltered_vcf | tee $unfiltered_log &&
    done

