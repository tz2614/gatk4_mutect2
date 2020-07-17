#!/bin/bash

# Use GATK4 vrsion 4.1.4.1 to call AmpliSeq based bam files in /users/tz1/GATK4_Mutect2/

reference_fa="/users/tz1/git/GATK4/refs/ucsc.hg19.HODS.fasta"
germline_vars="/users/tz1/git/GATK4/refs/af-only-gnomad.raw.sites.hg19.vcf.gz"
gatk="/users/tz1/apps/GATK4/gatk-4.1.4.1/gatk"
#vcf_format=".vcf.gz"
#log_format=".log"
unfiltered_vcf_format=".unfiltered.vcf.gz"
unzipped_unfiltered_vcf_format=".unfiltered.vcf"
filtered_log_format1=".filtered1.log"
filtered_vcf_format1=".filtered1.vcf"
data_dir="/users/tz1/git/GATK4/"
#test_data_dir="/users/tz1/git/GATK4/test_data/"
bed_interval_list="/users/tz1/git/GATK4/refs/hotspot_region_Ion_AmpliSeq_CHPv2_Somatic_Hotspots_GATK4_formatted3_output.interval_list"
interval_list="chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr7 -L chr9 -L chr10 -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20"

#cmd="$gatk Mutect2 -R $reference_fa -L $interval_list -XL chrM -I $bam -O $vcf"

for bam in $(find $data_dir -name "*.bam")
    do
        #vcf=${bam/".bam"/$vcf_format}
        unfiltered_vcf=${bam/".bam"/$unfiltered_vcf_format}
        unzipped_unfiltered_vcf=${bam/".bam"/$unzipped_unfiltered_vcf_format}
        unfiltered_log=${bam/".bam"/$log_format}
        filtered_log1=${bam/".bam"/$filtered_log_format1}
        filtered_vcf1=${bam/".bam"/$filtered_vcf_format1}
        echo "variant calling $bam using Mutect2"

        # Initial variant calling without filtering
        $gatk Mutect2 -R $reference_fa -L $interval_list -L $bed_interval_list -I $bam -O $unfiltered_vcf | tee $unfiltered_log &&

        #$gatk Mutect2 -R $reference_fa ---germline-resource $germline_vars -L $interval_list -I $bam -O $unfiltered_vcf | tee $unfiltered_log &&
        zcat $unfiltered_vcf > $unzipped_unfiltered_vcf &&

        # Filtering based using FilterMutectCalls
        echo "variant filtering using FilterMutectCalls" &&
        $gatk FilterMutectCalls -R $reference_fa -V $unfiltered_vcf -O $filtered_vcf1 | tee $filtered_log1 &&
        echo "$filtered_vcf1 created"
    done
