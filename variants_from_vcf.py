#!usr/bin/python2

from __future__ import print_function
import sys
import os
import datetime
import vcf
import csv

"""generate a CSV file for a given VCF"""

# author: Tengyue Zheng
# date: 28/01/2020

filedate = str(datetime.date.today())
filedate = "".join(filedate.split("-"))
print (filedate)

def get_variants_in_vcf(annotated_vcf):
	"""
	Args:
		vcf(str): VCF filepath
	Returns:
		csv(str): CSV filepath
	"""
	dir_path = os.path.dirname(annotated_vcf)
	error_log = os.path.join(dir_path, "vcf_parse_error_log.txt")
	filename = annotated_vcf.split(".vcf")[0]
	csv_output = filename + "_" + filedate + ".csv" 
	print (csv_output)
	header = ["Chromosome", "ID", "Position", "Reference Allele", "Alternative Allele", "Quality Metric", "Filter Metrics", "Annotation Info"]

	"""open the annotated vcf and search records that contain the gene of interest, return the record as a pyvcf Record object"""

	with open (annotated_vcf, "r") as vcf_file:

		"""handle errors generated using a try, except statement and record some of the common errors in error log, 
		if no errors are generated then no log is recorded."""

		try: 
			vcf_reader = vcf.Reader(vcf_file)

		except (TypeError, RuntimeError, NameError, ValueError) as e:
			with open(error_log, "a") as err_log:
				err_log.writelines(e)

		print ("START parsing VCF file")

		with open (csv_output, "w") as csv_fh:
			csv_w = csv.writer(csv_fh, delimiter=",", quotechar='"', quoting=csv.QUOTE_ALL)
			print ("writing header to csv: {0}".format(csv_output))
			csv_w.writerow(header)

			for record in vcf_reader:
				print (record)
				chrom = str(record.CHROM)
				ids = str(record.ID)
				pos = int(record.POS)
				ref = str(record.REF)
				alt = ",".join([str(x) for x in record.ALT])
				qual = str(record.QUAL)
				sieve = str(record.FILTER)
				mpos = ",".join([str(x) for x in record.INFO['MPOS']])
				mmq = ",".join([str(x) for x in record.INFO['MMQ']])
				ecnt = str(record.INFO['ECNT'])
				tlod = ",".join([str(x) for x in record.INFO['TLOD']])
				popaf = ",".join([str(x) for x in record.INFO['POPAF']])
				mfrl = ",".join([str(x) for x in record.INFO['MFRL']])
				mbq = ",".join([str(x) for x in record.INFO['MBQ']])
				dp = str(record.INFO['DP'])
				info = "MPOS:{0}, MMQ:{1}, ECNT:{2}, TLOD:{3}, POPAF:{4}, MFRL:{5}, MBQ:{6}, DP:{7}".format(mpos, mmq, ecnt, tlod, popaf, mfrl, mbq, dp)
				csv_w.writerow([chrom, ids, pos, ref, alt, qual, sieve, info])

	print ("parsing COMPLETE")
	print ("new csv file created: {0}".format(csv_output))
	return csv_output

def main(annotated_vcf):

	# check that the annotated_vcf as a input
	assert os.path.exists(annotated_vcf), "{} DO NOT exists".format(annotated_vcf)

	# create a CSV file contaiuning VCF variants
	csv_output = get_variants_in_vcf(annotated_vcf)

if __name__ == "__main__":
	main(sys.argv[1])

