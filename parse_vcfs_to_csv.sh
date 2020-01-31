#!/bin/bash

# Parse all VCF variants to CSV for analysis

# load modules
module "add" "apps/python/2.7.8/gcc-4.8.5"
module "add" "apps/virtualenv/16.0.0/python-2.7.8"

# assign variables
filtered_vcf_format1=".filtered1.vcf"
parse_log_format=".filtered1.parse.log"
#test_data_dir="/users/tz1/GATK_Mutect2/test_data"
data_dir="/users/tz1/GATK4_Mutect2/"
script=$data_dir"variants_from_Mutect2.py"
#echo $script

# activate venv containing pyvcf 0.6.8
cd $data_dir
source "venv/bin/activate"

echo "Parsing variants to CSV file"

for vcf in $(find $data_dir -name "*.filtered1.vcf")
	do
		parse_log=${vcf/$filtered_vcf_format1/$parse_log_format}
		echo "Parsing variants from $vcf to CSV"
		python2 $script $vcf | tee $parse_log
		echo "Parsing COMPLETE for $vcf"
	done
