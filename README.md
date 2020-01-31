## GATK4 Mutect2 Variant Calling
  - By Tengyue Zheng
  - 24/01/2020
  email: tengyue.zheng@mft.nhs.uk

## Description
  
  Generate GATK4 Mutect2 VCF files for all BAM files in a given directory

## Getting Started

  1. Install GATK4 on the cluster
  ```Bash
   wget --no-check-certificate https://github.com/broadinstitute/gatk/releases/download/4.1.4.1/gatk-4.1.4.1.zip
   
   unzip gatk-4.1.4.1.zip -d /users/tz1/apps/GATK4/
  ```

## Prerequisites
  1. Make sure you have all the dependencies for GATK4 beforehand
  check the following website for details:
  - https://github.com/broadinstitute/gatk/blob/master/README.md#requirements

  2. Perform a thorough code review for the script before executing the variant calling

## Main scripts:
  - GATK4_mutect2_variant_calling.sh
  - variants_from_Mutect2.py
  - parse_vcfs_to_csv.sh

## Reference files
  - ucsc.hg19.HODS.fasta
  - ucsc.hg19.HODS.dict
  - af-only-gnomad.raw.sites.hg19.vcf.gz
  - hotspot_region_Ion_AmpliSeq_CHPv2_Somatic_Hotspots_GATK4_formatted3_output.interval_list

## User Requirements:

  1. Generate VCF files using GATK4 and use FilterMutectCalls to remove variants of low quality or low likelihood of being somatic variants

## Instructions:
  
  2. Load dependencies

  ```Bash
  cd /users/tz1/GATK4_Mutect2/

  module add apps/virtualenv/16.0.0/python-2.7.8
  
  source venv/bin/activate
  ``` 
  3. Test following script on a test bam file

  ```Bash
  nano GATK4_mutect2_variant_calling.sh
  ```
  uncomment /users/tz1/GATK_Mutect2/test_data/

  ```Bash
  sh GATK_mutect2_variant_calling.sh
  ```
  4. Check the VCF files generated
  
  If the VCF files have the correct format and content, then perform the same step as step3, EXCEPT uncomment data_dir="/users/tz1/GATK_Mutect2/"

  5. Check the output, it should display the information on terminal:

  ```Bash
  variant calling /users/tz1/GATK4_Mutect2/TP53/Detected/MA5245_NGS190206_198.IonXpress_009.bam using Mutect2
  Using GATK jar /users/tz1/apps/GATK4/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar
  Running:
      java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /users/tz1/apps/GATK4/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar Mutect2 -R /users/tz1/GATK4_Mutect2/refs/ucsc.hg19.HODS.fasta -L /users/tz1/GATK4_Mutect2/refs/hotspot_region_Ion_AmpliSeq_CHPv2_Somatic_Hotspots_GATK4_formatted3_output.interval_list -XL chrM -I /users/tz1/GATK4_Mutect2/TP53/Detected/MA5245_NGS190206_198.IonXpress_009.bam -O /users/tz1/GATK4_Mutect2/TP53/Detected/MA5245_NGS190206_198.IonXpress_009.vcf.gz
  19:14:24.354 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/users/tz1/apps/GATK4/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar!/com/intel/gkl/native/libgkl_compression.so
  Jan 27, 2020 7:14:24 PM shaded.cloud_nio.com.google.auth.oauth2.ComputeEngineCredentials runningOnComputeEngine
  INFO: Failed to detect whether we are running on Google Compute Engine.
  19:14:24.563 INFO  Mutect2 - ------------------------------------------------------------
  19:14:24.564 INFO  Mutect2 - The Genome Analysis Toolkit (GATK) v4.1.4.1
  19:14:24.564 INFO  Mutect2 - For support and documentation go to https://software.broadinstitute.org/gatk/
  19:14:24.564 INFO  Mutect2 - Executing as tz1@login01.pri.nextgen3.alces.network on Linux v3.10.0-862.3.3.el7.x86_64 amd64
  19:14:24.564 INFO  Mutect2 - Java runtime: OpenJDK 64-Bit Server VM v1.8.0_171-b10
  19:14:24.564 INFO  Mutect2 - Start Date/Time: 27 January 2020 19:14:24 GMT
  19:14:24.564 INFO  Mutect2 - ------------------------------------------------------------
  19:14:24.564 INFO  Mutect2 - ------------------------------------------------------------
  19:14:24.565 INFO  Mutect2 - HTSJDK Version: 2.21.0
  19:14:24.565 INFO  Mutect2 - Picard Version: 2.21.2
  19:14:24.565 INFO  Mutect2 - HTSJDK Defaults.COMPRESSION_LEVEL : 2
  19:14:24.565 INFO  Mutect2 - HTSJDK Defaults.USE_ASYNC_IO_READ_FOR_SAMTOOLS : false
  19:14:24.565 INFO  Mutect2 - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_SAMTOOLS : true
  19:14:24.565 INFO  Mutect2 - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_TRIBBLE : false
  19:14:24.565 INFO  Mutect2 - Deflater: IntelDeflater
  19:14:24.565 INFO  Mutect2 - Inflater: IntelInflater
  19:14:24.565 INFO  Mutect2 - GCS max retries/reopens: 20
  19:14:24.565 INFO  Mutect2 - Requester pays: disabled
  19:14:24.565 INFO  Mutect2 - Initializing engine
  19:14:25.028 INFO  FeatureManager - Using codec IntervalListCodec to read file file:///users/tz1/GATK4_Mutect2/refs/hotspot_region_Ion_AmpliSeq_CHPv2_Somatic_Hotspots_GATK4_formatted3_output.interval_list
  19:14:25.187 INFO  IntervalArgumentCollection - Initial include intervals span 2232 loci; exclude intervals span 16569 loci
  19:14:25.187 INFO  IntervalArgumentCollection - Excluding 0 loci from original intervals (0.00% reduction)
  19:14:25.187 INFO  IntervalArgumentCollection - Processing 2232 bp from intervals
  19:14:25.203 INFO  Mutect2 - Done initializing engine
  19:14:25.221 INFO  NativeLibraryLoader - Loading libgkl_utils.so from jar:file:/users/tz1/apps/GATK4/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar!/com/intel/gkl/native/libgkl_utils.so
  19:14:25.223 INFO  NativeLibraryLoader - Loading libgkl_pairhmm_omp.so from jar:file:/users/tz1/apps/GATK4/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar!/com/intel/gkl/native/libgkl_pairhmm_omp.so
  19:14:25.270 INFO  IntelPairHmm - Using CPU-supported AVX-512 instructions
  19:14:25.270 INFO  IntelPairHmm - Flush-to-zero (FTZ) is enabled when running PairHMM
  19:14:25.270 INFO  IntelPairHmm - Available threads: 2
  19:14:25.270 INFO  IntelPairHmm - Requested threads: 4
  19:14:25.270 WARN  IntelPairHmm - Using 2 available threads, but 4 were requested
  19:14:25.270 INFO  PairHMM - Using the OpenMP multi-threaded AVX-accelerated native PairHMM implementation
  19:14:25.313 INFO  ProgressMeter - Starting traversal
  19:14:25.313 INFO  ProgressMeter -        Current Locus  Elapsed Minutes     Regions Processed   Regions/Minute
  19:14:30.834 INFO  Mutect2 - 12759 read(s) filtered by: (((((((((((MappingQualityReadFilter AND MappingQualityAvailableReadFilter) AND MappingQualityNotZeroReadFilter) AND MappedReadFilter) AND NotSecondaryAlignmentReadFilter) AND NotDuplicateReadFilter) AND PassesVendorQualityCheckReadFilter) AND NonChimericOriginalAlignmentReadFilter) AND NonZeroReferenceLengthAlignmentReadFilter) AND ReadLengthReadFilter) AND GoodCigarReadFilter) AND WellformedReadFilter)
    12759 read(s) filtered by: ((((((((((MappingQualityReadFilter AND MappingQualityAvailableReadFilter) AND MappingQualityNotZeroReadFilter) AND MappedReadFilter) AND NotSecondaryAlignmentReadFilter) AND NotDuplicateReadFilter) AND PassesVendorQualityCheckReadFilter) AND NonChimericOriginalAlignmentReadFilter) AND NonZeroReferenceLengthAlignmentReadFilter) AND ReadLengthReadFilter) AND GoodCigarReadFilter)
        12759 read(s) filtered by: (((((((((MappingQualityReadFilter AND MappingQualityAvailableReadFilter) AND MappingQualityNotZeroReadFilter) AND MappedReadFilter) AND NotSecondaryAlignmentReadFilter) AND NotDuplicateReadFilter) AND PassesVendorQualityCheckReadFilter) AND NonChimericOriginalAlignmentReadFilter) AND NonZeroReferenceLengthAlignmentReadFilter) AND ReadLengthReadFilter)
            12610 read(s) filtered by: ((((((((MappingQualityReadFilter AND MappingQualityAvailableReadFilter) AND MappingQualityNotZeroReadFilter) AND MappedReadFilter) AND NotSecondaryAlignmentReadFilter) AND NotDuplicateReadFilter) AND PassesVendorQualityCheckReadFilter) AND NonChimericOriginalAlignmentReadFilter) AND NonZeroReferenceLengthAlignmentReadFilter)
                12610 read(s) filtered by: (((((((MappingQualityReadFilter AND MappingQualityAvailableReadFilter) AND MappingQualityNotZeroReadFilter) AND MappedReadFilter) AND NotSecondaryAlignmentReadFilter) AND NotDuplicateReadFilter) AND PassesVendorQualityCheckReadFilter) AND NonChimericOriginalAlignmentReadFilter)
                    12610 read(s) filtered by: ((((((MappingQualityReadFilter AND MappingQualityAvailableReadFilter) AND MappingQualityNotZeroReadFilter) AND MappedReadFilter) AND NotSecondaryAlignmentReadFilter) AND NotDuplicateReadFilter) AND PassesVendorQualityCheckReadFilter)
                        12610 read(s) filtered by: (((((MappingQualityReadFilter AND MappingQualityAvailableReadFilter) AND MappingQualityNotZeroReadFilter) AND MappedReadFilter) AND NotSecondaryAlignmentReadFilter) AND NotDuplicateReadFilter)
                            12610 read(s) filtered by: ((((MappingQualityReadFilter AND MappingQualityAvailableReadFilter) AND MappingQualityNotZeroReadFilter) AND MappedReadFilter) AND NotSecondaryAlignmentReadFilter)
                                12610 read(s) filtered by: (((MappingQualityReadFilter AND MappingQualityAvailableReadFilter) AND MappingQualityNotZeroReadFilter) AND MappedReadFilter)
                                    12610 read(s) filtered by: ((MappingQualityReadFilter AND MappingQualityAvailableReadFilter) AND MappingQualityNotZeroReadFilter)
                                        12610 read(s) filtered by: (MappingQualityReadFilter AND MappingQualityAvailableReadFilter)
                                            12610 read(s) filtered by: MappingQualityReadFilter 
            149 read(s) filtered by: ReadLengthReadFilter 
  
  19:14:30.835 INFO  ProgressMeter -       chr22:24176287              0.1                   793           8616.4
  19:14:30.835 INFO  ProgressMeter - Traversal complete. Processed 793 total regions in 0.1 minutes.
  19:14:30.860 INFO  VectorLoglessPairHMM - Time spent in setup for JNI call : 0.005268424000000001
  19:14:30.860 INFO  PairHMM - Total compute time in PairHMM computeLogLikelihoods() : 0.21107677500000002
  19:14:30.860 INFO  SmithWatermanAligner - Total compute time in java Smith-Waterman : 0.36 sec
  19:14:30.861 INFO  Mutect2 - Shutting down engine
  [27 January 2020 19:14:30 GMT] org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2 done. Elapsed time: 0.11 minutes.
  Runtime.totalMemory()=613941248
  Tool returned:
  SUCCESS
  /users/tz1/GATK4_Mutect2/TP53/Detected/MA5245_NGS190206_198.IonXpress_009.vcf.gz created
  ```
  6. Test following script on a test vcf file

  ```Bash
  nano parse_vcfs_to_csv.sh
  ```
  uncomment /users/tz1/GATK_Mutect2/test_data/

  ```Bash
  sh parse_vcfs_to_csv.sh
  ```
  
  7. Check the output, it should display the information on terminal:

  ```Bash
  Parsing variants to CSV file
  Parsing variants from /users/tz1/GATK4_Mutect2/<sample>.filtered1.vcf to CSV
  START parsing VCF file
  writing header to csv: /mnt/repository/Bioinformatics/tengyue_zheng_projects/GATK4_Mutect2/Mutect2_output_YYYYMMDD.csv
  var imported from vcf: /users/tz1/GATK4_Mutect2/<sample>.filtered1.vcf
  csv_output: /mnt/repository/Bioinformatics/tengyue_zheng_projects/GATK4_Mutect2/Mutect2_output_YYYYMMDD.csv
  Parsing COMPLETE for /users/tz1/GATK4_Mutect2/<sample>.filtered1.vcf
  ```

  8. When the job is completed, check that log files have been generated.
  
  If in doubt consult Tengyue Zheng or senior member of the bioinformatics team.

  ```Bash
  $ less /users/tz1/<run_folder>/<sample>.filtered1.vcf.log
  $ less /users/tz1/<run_folder>/<sample>.log
  ```
   You should see the outputs as above
