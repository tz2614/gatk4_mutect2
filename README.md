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
  Using GATK jar /users/tz1/apps/GATK4/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar
  Running:
      java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /users/tz1/apps/GATK4/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar Mutect2 -R /users/tz1/GATK4_Mutect2/refs/ucsc.hg19.HODS.fasta -L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr7 -L chr9 -L chr10 -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 -I /users/tz1/GATK4_Mutect2/test_data/test.bam -O /users/tz1/GATK4_Mutect2/test_data/test.unfiltered.vcf.gz
  12:00:33.756 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/users/tz1/apps/GATK4/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar!/com/intel/gkl/native/libgkl_compression.so
  Feb 04, 2020 12:00:34 PM shaded.cloud_nio.com.google.auth.oauth2.ComputeEngineCredentials runningOnComputeEngine
  INFO: Failed to detect whether we are running on Google Compute Engine.
  12:00:34.008 INFO  Mutect2 - ------------------------------------------------------------
  12:00:34.008 INFO  Mutect2 - The Genome Analysis Toolkit (GATK) v4.1.4.1
  12:00:34.008 INFO  Mutect2 - For support and documentation go to https://software.broadinstitute.org/gatk/
  12:00:34.009 INFO  Mutect2 - Executing as tz1@login01.pri.nextgen3.alces.network on Linux v3.10.0-862.3.3.el7.x86_64 amd64
  12:00:34.009 INFO  Mutect2 - Java runtime: OpenJDK 64-Bit Server VM v1.8.0_171-b10
  12:00:34.009 INFO  Mutect2 - Start Date/Time: 04 February 2020 12:00:33 GMT
  12:00:34.009 INFO  Mutect2 - ------------------------------------------------------------
  12:00:34.009 INFO  Mutect2 - ------------------------------------------------------------
  12:00:34.009 INFO  Mutect2 - HTSJDK Version: 2.21.0
  12:00:34.009 INFO  Mutect2 - Picard Version: 2.21.2
  12:00:34.009 INFO  Mutect2 - HTSJDK Defaults.COMPRESSION_LEVEL : 2
  12:00:34.009 INFO  Mutect2 - HTSJDK Defaults.USE_ASYNC_IO_READ_FOR_SAMTOOLS : false
  12:00:34.009 INFO  Mutect2 - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_SAMTOOLS : true
  12:00:34.009 INFO  Mutect2 - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_TRIBBLE : false
  12:00:34.009 INFO  Mutect2 - Deflater: IntelDeflater
  12:00:34.009 INFO  Mutect2 - Inflater: IntelInflater
  12:00:34.010 INFO  Mutect2 - GCS max retries/reopens: 20
  12:00:34.010 INFO  Mutect2 - Requester pays: disabled
  12:00:34.010 INFO  Mutect2 - Initializing engine
  12:00:34.382 INFO  IntervalArgumentCollection - Processing 2464119736 bp from intervals
  12:00:34.388 INFO  Mutect2 - Done initializing engine
  12:00:34.407 INFO  NativeLibraryLoader - Loading libgkl_utils.so from jar:file:/users/tz1/apps/GATK4/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar!/com/intel/gkl/native/libgkl_utils.so
  12:00:34.409 INFO  NativeLibraryLoader - Loading libgkl_pairhmm_omp.so from jar:file:/users/tz1/apps/GATK4/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar!/com/intel/gkl/native/libgkl_pairhmm_omp.so
  12:00:34.449 INFO  IntelPairHmm - Using CPU-supported AVX-512 instructions
  12:00:34.450 INFO  IntelPairHmm - Flush-to-zero (FTZ) is enabled when running PairHMM
  12:00:34.451 INFO  IntelPairHmm - Available threads: 2
  12:00:34.451 INFO  IntelPairHmm - Requested threads: 4
  12:00:34.451 WARN  IntelPairHmm - Using 2 available threads, but 4 were requested
  12:00:34.451 INFO  PairHMM - Using the OpenMP multi-threaded AVX-accelerated native PairHMM implementation
  12:00:34.487 INFO  ProgressMeter - Starting traversal
  12:00:34.487 INFO  ProgressMeter -        Current Locus  Elapsed Minutes     Regions Processed   Regions/Minute
  12:00:44.487 INFO  ProgressMeter -        chr1:36806494              0.2                122700         736200.0
  12:00:54.487 INFO  ProgressMeter -        chr1:85303288              0.3                284370         853110.0
  12:01:04.487 INFO  ProgressMeter -       chr1:141384149              0.5                471310         942620.0
  12:01:14.487 INFO  ProgressMeter -       chr1:197226679              0.7                657460         986190.0
  12:01:24.488 INFO  ProgressMeter -       chr1:248444238              0.8                828200         993840.0
  12:01:34.487 INFO  ProgressMeter -        chr2:45037656              1.0                981020         981020.0
  12:01:44.487 INFO  ProgressMeter -        chr2:99120385              1.2               1161300         995400.0
  12:01:54.487 INFO  ProgressMeter -       chr2:151564000              1.3               1336120        1002090.0
  12:02:04.487 INFO  ProgressMeter -       chr2:203075784              1.5               1507840        1005226.7
  12:02:14.487 INFO  ProgressMeter -         chr3:4462201              1.7               1656480         993888.0
  12:02:24.487 INFO  ProgressMeter -        chr3:54838643              1.8               1824420         995138.2
  12:02:34.487 INFO  ProgressMeter -       chr3:105963641              2.0               1994840         997420.0
  12:02:44.487 INFO  ProgressMeter -       chr3:159316155              2.2               2172690        1002780.0
  12:02:54.487 INFO  ProgressMeter -         chr4:3096915              2.3               2312060         990882.9
  12:03:04.487 INFO  ProgressMeter -        chr4:53689402              2.5               2480710         992284.0
  12:03:14.487 INFO  ProgressMeter -        chr4:99196455              2.7               2632430         987161.3
  12:03:24.701 INFO  ProgressMeter -       chr4:153250800              2.8               2812620         991441.4
  12:03:34.701 INFO  ProgressMeter -        chr5:13132383              3.0               2982750         993069.4
  12:03:44.701 INFO  ProgressMeter -        chr5:66890120              3.2               3161950         997387.2
  12:03:54.704 INFO  ProgressMeter -       chr5:118679786              3.3               3334600         999295.8
  12:04:04.704 INFO  ProgressMeter -       chr5:169928454              3.5               3505440        1000520.4
  12:04:14.704 INFO  ProgressMeter -        chr7:40723400              3.7               3677830        1002056.2
  12:04:24.706 INFO  ProgressMeter -        chr7:91093359              3.8               3845740        1002282.2
  12:04:34.709 INFO  ProgressMeter -       chr7:131130752              4.0               3979220         993885.7
  12:04:44.709 INFO  ProgressMeter -        chr9:20616190              4.2               4141320         993035.0
  12:04:54.709 INFO  ProgressMeter -        chr9:67137055              4.3               4296400         990631.1
  12:05:04.709 INFO  ProgressMeter -       chr9:119046219              4.5               4469440         992392.9
  12:05:14.709 INFO  ProgressMeter -       chr10:26636744              4.7               4632140         991815.1
  12:05:24.709 INFO  ProgressMeter -       chr10:76913929              4.8               4799750         992292.1
  12:05:34.709 INFO  ProgressMeter -      chr10:127710709              5.0               4969090         993083.1
  12:05:44.709 INFO  ProgressMeter -       chr11:44306348              5.2               5142880         994683.8
  12:05:54.709 INFO  ProgressMeter -       chr11:96490927              5.3               5316840         996216.4
  12:06:04.709 INFO  ProgressMeter -        chr12:4854154              5.5               5461440         992321.5
  12:06:14.709 INFO  ProgressMeter -       chr12:55781985              5.7               5631210         993094.5
  12:06:24.709 INFO  ProgressMeter -      chr12:106522789              5.8               5800350         993715.4
  12:06:34.709 INFO  ProgressMeter -       chr13:24430675              6.0               5972900         994869.8
  12:06:44.709 INFO  ProgressMeter -       chr13:65535826              6.2               6109940         990207.0
  12:06:54.709 INFO  ProgressMeter -         chr14:630601              6.3               6277500         990605.5
  12:07:04.709 INFO  ProgressMeter -       chr14:52650182              6.5               6450910         991883.1
  12:07:14.709 INFO  ProgressMeter -      chr14:105429925              6.7               6626860         993477.6
  12:07:24.709 INFO  ProgressMeter -       chr15:50444166              6.8               6801420         994791.1
  12:07:34.709 INFO  ProgressMeter -      chr15:100921339              7.0               6969690         995144.0
  12:07:44.709 INFO  ProgressMeter -       chr16:48825812              7.2               7137820         995460.9
  12:07:54.714 INFO  ProgressMeter -        chr17:7578266              7.3               7301530         995149.8
  12:08:04.714 INFO  ProgressMeter -       chr17:54491300              7.5               7457930         993889.3
  12:08:14.714 INFO  ProgressMeter -       chr18:25409688              7.7               7631660         994942.9
  12:08:24.714 INFO  ProgressMeter -       chr18:75603198              7.8               7798990         995135.1
  12:08:34.714 INFO  ProgressMeter -       chr19:39232423              8.0               7938030         991784.7
  12:08:44.714 INFO  ProgressMeter -       chr20:32370765              8.2               8112270         992879.2
  12:08:50.864 INFO  Mutect2 - 14863 read(s) filtered by: (((((((((((MappingQualityReadFilter AND MappingQualityAvailableReadFilter) AND MappingQualityNotZeroReadFilter) AND MappedReadFilter) AND NotSecondaryAlignmentReadFilter) AND NotDuplicateReadFilter) AND PassesVendorQualityCheckReadFilter) AND NonChimericOriginalAlignmentReadFilter) AND NonZeroReferenceLengthAlignmentReadFilter) AND ReadLengthReadFilter) AND GoodCigarReadFilter) AND WellformedReadFilter)
    14863 read(s) filtered by: ((((((((((MappingQualityReadFilter AND MappingQualityAvailableReadFilter) AND MappingQualityNotZeroReadFilter) AND MappedReadFilter) AND NotSecondaryAlignmentReadFilter) AND NotDuplicateReadFilter) AND PassesVendorQualityCheckReadFilter) AND NonChimericOriginalAlignmentReadFilter) AND NonZeroReferenceLengthAlignmentReadFilter) AND ReadLengthReadFilter) AND GoodCigarReadFilter)
         14863 read(s) filtered by: (((((((((MappingQualityReadFilter AND MappingQualityAvailableReadFilter) AND MappingQualityNotZeroReadFilter) AND MappedReadFilter) AND NotSecondaryAlignmentReadFilter) AND NotDuplicateReadFilter) AND PassesVendorQualityCheckReadFilter) AND NonChimericOriginalAlignmentReadFilter) AND NonZeroReferenceLengthAlignmentReadFilter) AND ReadLengthReadFilter)
            14635 read(s) filtered by: ((((((((MappingQualityReadFilter AND MappingQualityAvailableReadFilter) AND MappingQualityNotZeroReadFilter) AND MappedReadFilter) AND NotSecondaryAlignmentReadFilter) AND NotDuplicateReadFilter) AND PassesVendorQualityCheckReadFilter) AND NonChimericOriginalAlignmentReadFilter) AND NonZeroReferenceLengthAlignmentReadFilter)
                14635 read(s) filtered by: (((((((MappingQualityReadFilter AND MappingQualityAvailableReadFilter) AND MappingQualityNotZeroReadFilter) AND MappedReadFilter) AND NotSecondaryAlignmentReadFilter) AND NotDuplicateReadFilter) AND PassesVendorQualityCheckReadFilter) AND NonChimericOriginalAlignmentReadFilter)
                    14635 read(s) filtered by: ((((((MappingQualityReadFilter AND MappingQualityAvailableReadFilter) AND MappingQualityNotZeroReadFilter) AND MappedReadFilter) AND NotSecondaryAlignmentReadFilter) AND NotDuplicateReadFilter) AND PassesVendorQualityCheckReadFilter)
                        14635 read(s) filtered by: (((((MappingQualityReadFilter AND MappingQualityAvailableReadFilter) AND MappingQualityNotZeroReadFilter) AND MappedReadFilter) AND NotSecondaryAlignmentReadFilter) AND NotDuplicateReadFilter)
                            14635 read(s) filtered by: ((((MappingQualityReadFilter AND MappingQualityAvailableReadFilter) AND MappingQualityNotZeroReadFilter) AND MappedReadFilter) AND NotSecondaryAlignmentReadFilter)
                                14635 read(s) filtered by: (((MappingQualityReadFilter AND MappingQualityAvailableReadFilter) AND MappingQualityNotZeroReadFilter) AND MappedReadFilter)
                                    14635 read(s) filtered by: ((MappingQualityReadFilter AND MappingQualityAvailableReadFilter) AND MappingQualityNotZeroReadFilter)
                                        14635 read(s) filtered by: (MappingQualityReadFilter AND MappingQualityAvailableReadFilter)
                                            14635 read(s) filtered by: MappingQualityReadFilter 
            228 read(s) filtered by: ReadLengthReadFilter 

  12:08:50.864 INFO  ProgressMeter -       chr20:63024336              8.3               8214473         992931.5
  12:08:50.865 INFO  ProgressMeter - Traversal complete. Processed 8214473 total regions in 8.3 minutes.
  12:08:50.915 INFO  VectorLoglessPairHMM - Time spent in setup for JNI call : 0.029467037
  12:08:50.915 INFO  PairHMM - Total compute time in PairHMM computeLogLikelihoods() : 8.782496036000001
  12:08:50.915 INFO  SmithWatermanAligner - Total compute time in java Smith-Waterman : 3.95 sec
  12:08:50.915 INFO  Mutect2 - Shutting down engine
  [04 February 2020 12:08:50 GMT] org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2 done. Elapsed time: 8.29 minutes.
  Runtime.totalMemory()=770703360
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
