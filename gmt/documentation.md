# Running BreakDancer

***

To run BreakDancer, first use bam2cfg.pl to prepare the required per-invocation config file. See the below for parameter descriptions.

<p class='terminal' markdown='1'>
/usr/lib/breakdancer-max1.2/bam2cfg.pl bam_files breakdancer_options
</p>
Then run BreakDancer on the config

<p class='terminal' markdown='1'>
breakdancer-max config_file.cfg
</p>

# EXAMPLE PIPELINE

***

Please see below for detailed information about the *bam2cfg.pl* and *breakdancer-max* commands.

## STEP 1
Create a configuration file using bam2cfg.pl:

The precompiled Debian package will install this in /usr/lib/breakdancer-max{{page.version_suffix}} along with a few required Perl modules.

<p class='terminal' markdown='1'>
/usr/lib/breakdancer-max{{page.version_suffix}}/bam2cfg.pl -g -h tumor.bam normal.bam > BRC6.cfg
</p>

## STEP 2
Detect inter-chromosomal translocations:

<p class='terminal' markdown='1'>
breakdancer_max -t -q 10 -d BRC6.ctx BRC6.cfg > BRC6.ctx
</p>

The -d option dumps CTX supporting read pairs into fastq files (in this case BRC6.ctx) by library.

This step normally takes 12 hours or so for three bam files, 8 hours or so for two bam files for cpp version, around three days for perl version.

# bam2cfg

***

## NAME
bam2cfg - create a configuration file for BreakDancer

## SYNOPSIS
	bam2cfg.pl [options] input_bams > config_file

## NOTES
Manually view the insert size and flag distribution results in the output .cfg file to see if there are any data quality issue. Usually std/mean should be < 0.2 or 0.3 at most. The flag 32(x%), represents percent of chimeric insert, this number (x%) should usually be smaller than 3%.

View png files for the insert size distribution. You should usually see a normal distribution, a bimodal distribution is undesirable and it is not recommended to continue BreakDancerMax step with this situation existing.

## OPTIONS
<dl>
<dt>-q INT</dt>
<dd>Minimum mapping quality [default = 35]</dd>

<dt>-m</dt>
<dd>Using mapping quality instead of alternative mapping quality</dd>

<dt>-s</dt>
<dd>Minimal mean insert size [default = 50]</dd>

<dt>-C</dt>
<dd>Change default system from Illumina to SOLiD</dd>

<dt>-c FLOAT</dt>
<dd>Cut off in unit of standard deviation [default = 4]</dd>

<dt>-n INT</dt>
<dd>Number of observation required to estimate mean and s.d. insert size [10000]</dd>

<dt>-v FLOAT</dt>
<dd>Cutoff on coefficients of variation [default = 1]</dd>

<dt>-f STRING</dt>
<dd>A two column tab-delimited text file (RG, LIB) specify the RG=>LIB mapping, useful when BAM header is incomplete</dd>

<dt>-b INT</dt>
<dd>Number of bins in the histogram [default = 50]</dd>

<dt>-g</dt>
<dd>Output mapping flag distribution</dd>

<dt>-h</dt>
<dd>Plot insert size histogram for each BAM library</dd>
</dl>

# BreakDancerMax

***

## NAME
BreakDancerMax - SV detection

## SYNOPSIS

breakdancer-max [options] config_file

## OPTIONS
<dl>
<dt>-o STRING</dt>
<dd>operate on a single chromosome [default is all chromosomes]</dd>

<dt>-c INT</dt>
<dd>cutoff in unit of standard deviation [default = 3]</dd>

<dt>-t</dt>
<dd>only detect transchromosomal rearrangement</dd>

<dt>-p FLOAT</dt>
<dd>prior probability of SV [default = 0.001] (Not available in C++ version)</dd>

<dt>-e</dt>
<dd>learn parameters from data before applying to SV detection (Not available in C++ version)</dd>

<dt>-l</dt>
<dd>analyze Illumina long insert (mate-pair) library</dd>

<dt>-f</dt>
<dd>use Fisher's method to combine P values from multiple library (Not available in C++ version)</dd>

<dt>-d STRING</dt>
<dd>prefix of fastq files that SV supporting reads will be saved by library</dd>

<dt>-g _FILE_</dt>
<dd>dump SVs and supporting reads in BED format for GBrowse</dd>

<dt>-m INT</dt>
<dd>maximum SV size [default = 1000000]</dd>

<dt>-q INT</dt>
<dd>minimum alternative mapping quality [default = 35]</dd>

<dt>-s INT</dt>
<dd>minimum length of a region [default = 7]</dd>

<dt>-b INT</dt>
<dd>buffer size for building connection [default = 100]</dd>

<dt>-r INT</dt>
<dd>minimum number of read pairs required to establish a connection</dd>

<dt>-x INT</dt>
<dd>maximum threshold of haploid sequence coverage for regions to be ignored [default = 1000]</dd>

<dt>-a</dt>
<dd>print out copy number by bam file rather than library [defaults on]</dd>

<dt>-h</dt>
<dd>print out Allele Frequency column [defaults off]</dd>

<dt>-y INT</dt>
<dd>output score filter [default = 30]</dd>
</dl>

## DESCRIPTION
Most of these options are self-explanatory. It is convenient to use the -o option to parallelize SV detection for each chromosome. When -o is used, the detection of inter-chromosomal translocation is disabled. In that case, it may be convenient to use -t in a separate process to detect putative inter-chromosomal translocations without bothering to analyze read pairs that are mapped to the same chromosome.

BreakDancer only supports properly formatted bam files and has only been tested using bam files produced by BWA. To obtain the correct result, it is important to have readgroup (@RG) tag in both the header and each alignment in the bam files. 

The input to breakdancer-max is a set of map files produced by a front-end aligner such as MAQ, BWA, NovoAlign and Bfast, and a tab-delimited configuration file that specifies the locations of the map files, the detection parameters, and the sample information.

### CONFIGURATION
If your map files are in the sam/bam format, you can use the bam2cfg.pl in the released package to automatic generate a configuration file (bam2cfg.pl also has dependence on AlnParser.pm in the release package). If you have a single bam file that contains multiple libraries, make sure that the readgroup and library information are properly encoded in the sam/bam header, and in each alignment record, otherwise bam2cfg.pl may fail to produce a correct configuration file. Please follow instructions on [SourceForge](http://samtools.sourceforge.net) to properly format your bam files.

An example manual configuration file is like this 

<p class='terminal' markdown='1'>
map:1.map mean:219 std:18 readlen:36.00 sample:tA exe:maq-0.6.8 mapview -b 
map:2.map mean:220 std:19 readlen:36.00 sample:tB exe:maq-0.6.8 mapview -b 
map:3.map mean:219 std:18 readlen:36.00 sample:nA exe:maq-0.7.1 mapview -b 
map:4.map mean:219 std:18 readlen:36.00 sample:nB exe:maq-0.7.1 mapview -b 
</p>

An example configuration file produced by bam2cfg.pl look like this: 

<p class='terminal' markdown='1'>
readgroup:2825107881 platform:illumina map:tumor.bam readlen:75.00 lib:demolib1
num:10001 lower:86.83 upper:443.91 mean:315.09 std:43.92 exe:samtools view 
readgroup:2843249908 platform:illumina map:tumor.bam readlen:75.00 lib:demolib1
num:10001 lower:86.83 upper:443.91 mean:315.09 std:43.92 exe:samtools view 
readgroup:2843255910 platform:illumina map:normal.bam readlen:75.00 lib:demolib2
num:10001 lower:95.36 upper:443.31 mean:311.68 std:42.86 exe:samtools view 
readgroup:2843255906 platform:illumina map:normal.bam readlen:75.00 lib:demolib2
num:10001 lower:95.36 upper:443.31 mean:311.68 std:42.86 exe:samtools view 
</p>

Each row must contain at least 6 key:value pairs (separated by colon) that specify:

1.	the location of the map file
2.	the mean insert size
3.	the standard deviation insert size
4.	the average read length
5.	a unique identifier assigned to the map file (usually representing a PE library)
6.	a command line that can run by perl system calls to produce MAQ mapview alignment

Listing multiple map files in a single configuration file would automatically enable pooled analysis: reads from all the map files are jointly analyzed to find unified SV hypotheses across all the map files.

### SEPARATION THRESHOLDS
In addition to the above 6 keys: map, mean, std, readlen, sample, and exe, BreakDancerMax allows users to explicitly specify the separation thresholds using the keys: upper and lower. For example:

<p class='terminal' markdown='1'>
map:1.map upper:300 lower:100 readlen:36.00 sample:tA exe:maq-0.6.8 mapview -b
</p>

This will instruct BreakDancerMax to detect deletions using read pairs that are at least 300 bp apart (outer distance) and detect insertions using read pairs that are at most 100 bp apart.

The upper and the lower key:value pairs, when explicitly specified, take precedence over the upper and the lower thresholds computed from the mean, the std, and the user specified threshold in the unit of standard deviation.

upper: 

        mean + std * threshold specified by user option -c

lower: 

        mean - std * threshold specified by user option -c

### OPTION DESCRIPTIONS
The -c option by default equals to 3. Therefore, the upper and the lower separation threshold would be: mean + 3 std and mean - 3 std respectively. It is useful to explicitly specify the upper and the lower separation thresholds when the insert size distribution is not symmetric to the mean. 

The -o option enables per-chromosome/reference analysis and is much faster when the input files are in the bam format. Please index the bam file using "samtools index" to utilize this option. You need to specify the exact reference names as they are in the bam files. 

When -e is on, BreakDancerMax tries to estimate the mean and the standard deviation insert size from the data instead of relying on user's spec in the configuration file. Current implementation of this estimation process is slow. So it is recommended that users can specify the accurate thresholds in the configuration file. 

The -l option tell BreakDancerMax that the data is produced from Illumina long insert circularized library 

The -f option uses the Fisher's methods to summarize scores from multiple libraries. It is recommended when there are many libraries. It ensures that the scores are independent of the number of libraries (uniform distribution of the P values) 

The -q specifies the MAQ mapping quality threshold and can be used to skip reads that are not confidently mapped. 

The -s specifies the minimal required size of a SV anchoring region from which the anomalously mapped reads are found. This parameter has some small effects on the SV detection accuracy. Increasing -s improves the specificity but also reduces the sensitivity. The default 7 bp seemed to work well. 

The -b parameter specifies the number of anomalous regions resides in the RAM before SV hypotheses begin to form among these regions. The default works well in general. For dataset that is exceptionally large, it may be helpful to reduce it to cut the resident RAM usage. 

The -d specifies a fastq file where all SV supporting reads will be saved in the fastq format. These reads can be realigned by other aligners such as novoalign, and then reanalyzed by BreakDancer.

