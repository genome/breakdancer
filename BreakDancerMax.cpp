#include <iostream>
#include <fstream>
//#include <string>
#include <vector>
#include <stdio.h>
#include <cstdlib>
#include <list>
#include <iomanip>
#include <cmath>
#include <math.h>
#include <time.h>
#include <cstring>
#include <map>
#include "BreakDancerMax.h"
#include "AlnParser.h"
#include "sam.h"
#include "bam.h"

using namespace std;

// open the bam file
int main(int argc, char *argv)
{
	string version = "BreakDancerMax-0.0.1r81";
	int c;
	string bam_file;
	
	string chr; // chromosome
	int min_len = 7;
	int cut_sd = 3;
	int max_sd = 1000000;
	int min_map_qual = 35;
	int min_read_pair = 2;
	int buffer_size = 100;
	int learn_par = 0;//bool
	float prior_prob = 0.001;
	int transchr_rearrange = 0;//bool
	int fisher = 0;//bool
	string prefix_fastq;
	string dump_BED;
	int Illumina_long_insert = 0;// bool
	int Illumina_to_SOLiD = 0;// bool
	while((c = getopt(argc, argv, "o:s:c:m:q:r:b:ep:tfd:g:lC")) >= 0){
		switch(c) {
			case 'o': chr = strdup(optarg); break;
			case 's': min_len = atoi(optarg); break;
			case 'c': cut_sd = atoi(optarg); break;
			case 'm': max_sd = atoi(optarg); break;
			case 'q': min_map_qual = atoi(optarg); break;
			case 'r': min_read_pair = atoi(optarg); break;
			case 'b': buffer_size = atoi(optarg); break;
			case 'e': learn_par = atoi(optarg); break;
			case 'p': prior_prob = atof(optarg); break;
			case 't': transchr_rearrange = atoi(optarg); break;
			case 'f': fisher = atoi(optarg); break;
			case 'd': prefix_fastq = strdup(optarg); break;
			case 'g': dump_BED = strdup(optarg); break;
			case 'l': Illumina_long_insert = atoi(optarg); break;
			case 'C': Illumina_to_SOLiD = atoi(optarg); break;
			
			default: fprintf(stderr, "Unrecognized option '-%c'.\n", c);
			return 1;
		}
	}
	if(optind == argc){
		fprintf(stderr, "\n");
		fprintf(stderr, "BreakDancerMax.pl <analysis_config.lst>\n\n");
		fprintf(stderr, "Options: \n");
		fprintf(stderr, "	-o STRING	operate on a single chromosome [all chromosome]\n");
		fprintf(stderr, "	-s INT	minimum length of a region [%d]\n", min_len);		 
		fprintf(stderr, "	-c INT	cutoff in unit of standard deviation [%d]\n", cut_sd);		
		fprintf(stderr, "	-m INT	maximum SV size [%d]\n", max_sd);		 
		fprintf(stderr, "	-q INT	minimum alternative mapping quality [%d]\n", min_map_qual);	
		fprintf(stderr, "	-r INT	minimum number of read pairs required to establish a connection [%d]\n", min_read_pair);		 
		fprintf(stderr, "	-b INT	buffer size for building connection [%d]\n", buffer_size);		
		fprintf(stderr, "	-e INT	learn parameters from data before applying to SV detection [%d]\n", learn_par);		 
		fprintf(stderr, "	-p FLOAT	prior probability of SV [%f]\n", prior_prob);	
		fprintf(stderr, "	-t INT	only detect transchromosomal rearrangement [%d]\n", transchr_rearrange);		 
		fprintf(stderr, "	-f INT	use Fisher's method to combine P values from multiple library [%d]\n", fisher);		
		fprintf(stderr, "	-d STRING	prefix of fastq files that SV supporting reads will be saved by library\n");		 
		fprintf(stderr, "	-g INT	dump SVs and supporting reads in BED format for GBrowse[%d]\n", dump_BED);
		fprintf(stderr, "	-l INT	analyze Illumina long insert (mate-pair) library [%d]\n", Illumina_long_insert);		 
		fprintf(stderr, "	-C INT	change system default from Illumina to SOLiD [%d]\n", Illumina_to_SOLiD);
		fprintf(stderr, "Version: %s\n", version);
		fprintf(stderr, "\n");
		return 1;
	}
	
	// define the readgroup_platform map
	map<int, string> readgroup_platform;
    //readgroup_platform[*] = "*";
	
	/************************************************* read bam file *****************************************************/
	////////////////////////// suppose we already know bam file to open and parse (for all the chromosomes)
	samfile_t *in = 0;
   	char in_mode[5], *fn_list = 0, *fn_ref = 0, *fn_rg = 0;
   	strcpy(in_mode, "r");
	// generate the fn_list if necessary (not using because we don't have reference now)
	// if (fn_list == 0 && fn_ref) fn_list = samfaipath(fn_ref);
	// open file handlers
	if ((in = samopen(argv[optind], in_mode, fn_list)) == 0) {
		fprintf(stderr, "[main_samview] fail to open file for reading.\n");
		goto file_end;
	}
	if (in->header == 0) {
		fprintf(stderr, "[main_samview] fail to read the header.\n");
		goto file_end;
	}
	// convert/print the entire file
	bam2_t *b = bam_init2();
	int r;
	while ((r = samread(in, b)) >= 0) { // read one alignment from `in'
		//************************** note: in this loop, do everything (algorithm) inside ************************** //
		// Different from perl, in C/CPP, we are not going to read all the lines of bam file, store data, and then process
		// Instead, we'd better read one line, process this line, and another
		// That's why we have to put this samread code in the main function, also, some of the BreakDancer Algorithm would be better if written to another function
		// so that when repeat calling (at least four options: 1. SAM all chr; 2. SAM one chr; 3. MAQ all chr; 4. MAQ one chr), the code would not be repeated,
		// but only one function will be called. 
		// Of course, those code the same to all the four options can still remain in the main function.
		
		// pass b to AlnParser and let it handle the rest.
		AlnParser(b, format, readgroup_platform, alt);
		
		// process data in b		 
	}
	if (r < -1) fprintf(stderr, "[main_samview] truncated file.\n");
	bam_destroy1(b);
		
file_end:
	free(fn_list); free(fn_ref); free(fn_rg);
	samclose(in);
	
	return 0;
}			

// first of all, we need to use cpp since it has better alternative (map) for the hash in perl
// second, for the rest of the code, we should use the following correspondence, since the data structure changed. Left column: perl. Right column: cpp
/*
t.readname 		:		bam1_qname(b)  (length: b.core.l_qname)
t.flag			:		b.core.flag
t.chr			:		b.core.tid (this transit from char to int)
t.pos			:		b.core.pos
t.qual			:		b.core.qual
mchr			:		b.core.mtid (this transit from char to int)
mpos			:		b.core.mpos
t.dist			:		? don't know yet.
t.seq			:		bam1_seq(b)		(length: b.core.l_qseq)
t.basequal		:		bam1_qual(b)	(length: b.core.l_qseq)
t.ori			:		b.ori	(in samtools, can be derived by bam1_strand(b), bam1_mstrand(b))
*/
