#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdio.h>
#include <cstdlib>
#include <list>
#include <iomanip>
#include <cmath>
#include <math.h>
#include <time.h>
#include <map>
#include "BreakDancerMax.h"
#include "AlnParser.h"
#include "sam.h"
#include "bam.h"

using namespace std;

// open the bam file
int main(int argc, char *argv[])
{
	string version ("BreakDancerMax-0.0.1r81");
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
		fprintf(stderr, "	-g STRING	dump SVs and supporting reads in BED format for GBrowse\n");
		fprintf(stderr, "	-l INT	analyze Illumina long insert (mate-pair) library [%d]\n", Illumina_long_insert);		 
		fprintf(stderr, "	-C INT	change system default from Illumina to SOLiD [%d]\n", Illumina_to_SOLiD);
		//fprintf(stderr, "Version: %s\n", version);
		fprintf(stderr, "\n");
		return 1;
	}
	

	// define the map SVtype
	std::map<int, std::string> SVtype;
	if(Illumina_long_insert == 1){
		SVtype[1] = "INV";
		SVtype[3] = "INS";
		SVtype[4] = "DEL";
		SVtype[8] = "INV";
		SVtype[32] = "CTX";
	}
	else{
		SVtype[1] = "INV";
		SVtype[2] = "DEL";
		SVtype[3] = "INS";
		SVtype[4] = "ITX";
		SVtype[8] = "INV";
		SVtype[32] = "CTX";
	}
	
	//?AP
	int LZERO = -99;
	float ZERO = exp(LZERO);
	map<string,string> exes;
	map<string,string> fmaps;
	map<string,string> libmaps;
	map<string,float> mean_insertsize;
	map<string,float> std_insertsize;
	map<string,float> uppercutoff;
	map<string,float> lowercutoff;
	map<string,float> readlens;
	map<string,int> mapQual;// global
	int max_readlen = 0;
	map<Keys, int> x_readcounts;
	map<string,string> readgroup_library;
	// define the readgroup_platform map
	std::map<char *, std::string> readgroup_platform;	
	map<string,ofstream> ReadsOut;
	int d = 1e10;// global
	
	// configure file
	ifstream CONFIG(argv[optind]);
	string line; // each line of the config file
	if(CONFIG.is_open())
	{
		while(! CONFIG.eof())
		{
			getline(myfile,line);
			// analyze the line
			string fmap = get_from_line(line,"map",1);
			float mean = get_from_line(line,"mean",0);
			float std = get_from_line(line,"std",0);
			float readlen = get_from_line(line,"readlen",0);
			float upper = get_from_line(line,"upp",0);
			float lower = get_from_line(line,"low",0);
			int mqual = get_from_line_two(line,"map","qual",0,0);// -1 if none
			string lib = get_from_line(line,"lib",0);
			if(lib.compare("NA")==0)//may be wrong here to check if lib has been defined or not
				lib = get_from_line(line,"samp",0);
			
			string readgroup = get_from_line(line,"group",1);
			if(readgroup.compare("NA")==0)
				readgroup = lib;
			readgroup_library[readgroup] = lib;
			
			string platform = get_from_line(line,"platform",1);
			if(Illumina_to_SOLiD == 1)
				readgroup_platform[readgroup] = "solid";
			else
				readgroup_platform[readgroup] = "illumina";
			readgroup_platform[readgroup] = platform;
			
			string exe = get_from_line(line,"exe",0);
			if(prefix_fastq != ""){
				ofstream ReadsOut[lib.append("1")](prefix_fastq.append(lib).append(".1.fastq"));
				if(!ReadsOut[lib.append("1")].is_open())
					cout << "unable to open " << lib << ".1.fastq, check write permission\n";
				ofstream ReadsOut[lib.append("2")](prefix_fastq.append(lib).append(".2.fastq"));
				if(!ReadsOut[lib.append("2")].is_open())
					cout << "unable to open " << lib << ".2.fastq, check write permission\n";					
			}// need to figure out map ReadsOut standard			
			
			libmaps[lib] = fmap;
			if(mqual != -1)
				mapQual[lib] = mqual;
			fmaps[fmap] = lib;
			
			if(mean != -1 && std != -1 && upper == -1 && lower == -1){
				upper = mean + std*cut_sd;
				lower = mean - std*cut_sd;
				lower = lower > 0 ? lower:0;
			}
			
			max_readlen = max_readlen < readlen ? readlen:max_readlen;
			
			mean_insertsize[lib] = mean;
			std_insertsize[lib] = std;
			uppercutoff[lib] = upper;
			lowercutoff[lib] = lower;
			readlens[lib] = readlen;
			
			//skip exes[fmap]
			
			int tmp = mean - readlen*2;	// this determines the mean of the max of the SV flanking region
			d = d<tmp ? d:tmp;
		}
	}
	CONFIG.close();
		
	if(d < 50)
		d = 50;
	
	if(dump_BED != "")
		ofstream BED(dump_BED);
	vector<string> format; 
	vector<string>::iterator it_format;

	map<> cmds; 	//need to initialize
 	
	// go through the iteration of fmaps
	map<string,string>::iterator ii;
	for(ii=fmaps.begin(); ii!=fmaps.end(); ++ii)
	{
 	    string exe = exes[(*ii).first];
 	    cmds[exe] = 0;
 	    if(exe.find("maq")!=string::npos)
 	      	format.insert(format.end(),1,"maq");
 	    else
 	      	format.insert(format.end(),1,"sam");
 	}
 	        
 	if(learn_par == 1)    
 	    EstimatePriorParameters();//need to work on the i/o of this function
 	        
 	int reference_len = 1;
 	map<string, int> nreads;
 	int defined_all_readgroups = 1;
 	    
	samfile_t *in = 0;
	char in_mode[5], *fn_list = 0, *fn_ref = 0, *fn_rg = 0;
	strcpy(in_mode, "r");
 	for(ii=fmaps.begin(); ii!=fmaps.end(); ++ii)
 	{
 		int ref_len = 0;
 		cmds[exe] ++;
 		
 		int p_pos = 0;
 		char p_chr;
 		
 	   	/****************** read bam file *****************/
		// suppose we already know bam file to open and parse (for all the chromosomes)
		// generate the fn_list if necessary (not using because we don't have reference now)
		// if (fn_list == 0 && fn_ref) fn_list = samfaipath(fn_ref);
		// open file handlers
		if ((in = samopen((*ii).first, in_mode, fn_list)) == 0) {
			fprintf(stderr, "[main_samview] fail to open file for reading.\n");
			continue;
		}
		if (in->header == 0) {
			fprintf(stderr, "[main_samview] fail to read the header.\n");
			continue;
		}
		// convert/print the entire file
		bam1_t *b = bam_init1();
		int r;
		while ((r = samread(in, b)) >= 0) { // read one alignment from `in'
   		//************** note: in this loop, do everything (algorithm) inside ******** //
		// Different from perl, in C/CPP, we are not going to read all the lines of bam file, store data, and then process
		// Instead, we'd better read one line, process this line, and another
		// That's why we have to put this samread code in the main function, also, some of the BreakDancer Algorithm would be better if written to another function
		// so that when repeat calling (at least four options: 1. SAM all chr; 2. SAM one chr; 3. MAQ all chr; 4. MAQ one chr), the code would not be repeated,
		// but only one function will be called. 
		// Of course, those code the same to all the four options can still remain in the main function.
		
		// pass b to AlnParser and let it handle the rest.
		
		//// have to deal with a certain chromosome here;
			char *readgroup;
			string format;
			string alt;
			char ori = AlnParser(b, format, alt, readgroup, readgroup_platform);
			// skip one chrmosome at the moment
			if(b->core.tid == p_chr)
				ref_len += b->core.pos - p_pos;
			p_pos = b->core.pos;
			p_chr = b->core.tid;
			string lib;
			if(readgroup.length() != 0)
				lib = readgroup_library[readgroup];
			else{
				defined_all_readgroups = 0;
				lib = fmaps[*ii.first];
			}
			if(lib.length() == 0)
				continue;
				
			nreads[lib] ++;	// need to initialize
			if(mapQual[lib] != 0 && b->core.qual <= mapQual[lib])
				continue;
			else if(b->core.qual <= min_map_qual)
				continue;
			if(b->core.flag == 0)
				continue;
			if(transchr_rearrange && b->core.flag < 32 || b->core.flag >= 64)
				continue;
				
			if(Illumina_long_insert){
				if(abs(b->core.isize) > uppercutoff[lib] && b->core.flag == 20))
					b->core.flag = 4;
				if(abs(b->core.isize) < uppercutoff[lib] && b->core.flag == 4))
					b->core.flag = 20;
				if(abs(b->core.isize) < lowercutoff[lib] && b->core.flag == 20))
					b->core.flag = 3;
			}
			else{
				if(abs(b->core.isize) > uppercutoff[lib] && b->core.flag == 18))
					b->core.flag = 2;
				if(abs(b->core.isize) < uppercutoff[lib] && b->core.flag == 2))
					b->core.flag = 18;
				if(abs(b->core.isize) < lowercutoff[lib] && b->core.flag == 18))
					b->core.flag = 3;
			}
			
			if(b->core.flag == 18 || b->core.flag == 20 || b->core.flag == 130)
				continue;
			
			// here we need a two keys map
			Keys keys(b->core.flag, lib);
			x_readcounts[keys] ++;	// need to work on the initialization here
				
			// process data in b		 
		}
		if (r < -1) fprintf(stderr, "[main_samview] truncated file.\n");
		if(ref_len == 0)
			fprintf(stderr, "Unable to decode %s. Please check that you have the correct paths and the bam files are indexed.", *ii);
		if(reference_len < ref_len)
			reference_len = ref_len;
		bam_destroy1(b);
		samclose(in);		
	}
	
	int merge = (cmds.size()==1 && defined_all_readgroups)?1:0;
	
	float total_phy_cov = 0;
	float total_seq_cov = 0;
	
	// map has already been sorted by the key
	// recflags = x_readcounts[keys the first key]
	map<string,int>::iterator nreads_ii;
	for(nreads_ii=nreads.begin(); nreads_ii!=nreads.end(); ++nreads_ii)
	{
		string lib = *nreads_ii.first;
		float sequence_coverage = lib*readlen[lib]/reference_len;
		total_seq_cov += sequence_coverage;
		float physical_coverage = lib*mean_insertsize[lib]/(2*reference_len);
		total_phy_cov += physical_coverage;
		
		if(x_readcounts[2][lib])
			int nread_lengthDiscrepant = x_readcounts[2][lib];
		if(x_readcounts[3][lib])
			nread_lengthDiscrepant += x_readcounts[3][lib]
		float tmp = (__if_exists(nread_lengthDiscrepant) && nread_lengthDiscrepant > 0)?reference_len/nread_lengthDiscrepant:50;// still need to test if exists
		d = d<tmp?d:tmp;
		
		printf("#%s\tmean:%.3f\tstd:%.3f\tuppercutoff:%.3f\tlowercutoff:%.3f\treadlen:%.3f\tlibrary:%s\treflen:%d\tseqcov:%.3fx\tphycov:%.3fx", libmaps[$lib],mean_insertsize[lib],std_insertsize[lib],uppercutoff[lib],lowercutoff[lib],readlens[lib],lib,reference_len, sequence_coverage,physical_coverage);
		
		map<Keys,int>::iterator x_readcounts_ii;
		for(x_readcounts_ii = x_readcounts.begin(); x_readcounts_ii!=x_readcounts.end(); ++x_readcounts_ii){
			uint32_t t = *x_readcounts_ii.first.key1;// get the first key out, which is a member of recflags
			printf("\t%s:%d",t,*x_readcounts_ii.second());// don't know why || 0, need to ask
		}
		printf("\n");
	}
	printf("#Chr1\tPos1\tOrientation1\tChr2\tPos2\tOrientation2\tType\tSize\tScore\tnum_Reads\tnum_Reads_lib\tAllele_frequency\tVersion\tRun_Param\n");
	
	string begins;// global
	int beginc = -1;// global
	string lasts;// global (chr, should be int in samtools)
	int lastc = -1; // global
	map<int, vector<string>> regs;//global in analysis
	map<string, vector<int>> read;// global in analysis
	map<int,string> reg_name;// global in analysis
	vector<string> reg_seq; // global need to see if it's the key or value of one of the above global. should be a string
	vector<string>::iterator it_reg_seq; // global
	
	int idx_buff = 0;// global
	int final_buff = 0;
	int reg_idx = 0;// global
	int normal_switch = 0; // global
	int nnormal_reads = 0; // global
	
	if(merge && fmaps.size()>1 && !(format[0].compare("sam")) && chr.empty() ){
 /* open pipe, improvement made by Ben Oberkfell (boberkfe@genome.wustl.edu)
   samtools merge - in1.bam in2.bam in3.bam in_N.bam | samtools view - 
   maq mapmerge		*/
   		string merge_command_line; // may not need to write the string, but need samtools functions of merge bam files
   		// merge_fh; //file handle of the merge bam file
   		
   			if ((in = samopen(?(merged bam name), in_mode, fn_list)) == 0) {
				fprintf(stderr, "[main_samview] fail to open file for reading.\n");
				continue;
			}
			if (in->header == 0) {
				fprintf(stderr, "[main_samview] fail to read the header.\n");
				continue;
			}
			// convert/print the entire file
			int r;
			while ((r = samread(in, b)) >= 0) { // read one alignment from `in'
				char *readgroup;
				string format;
				string alt;
				char ori = AlnParser(b, format, alt, readgroup, readgroup_platform);
				string library = readgroup.empty()?readgroup_library[readgroup]:fmaps.
			  	if(chr.empty() || chr.compare(b->core.tid)!=0)
			  		continue;
			  	if(!library.empty())
			  		Analysis(library, b);
			}
			bam_destroy1(b);
			samclose(in);		
	}
	else{
		// define FHs
		map<> Idxs;
		map<int, string> buffer; // have to define it here because the two for loop are put together here
		int i = 0;
		for(ii=fmaps.begin(); ii!=fmaps.end(); ++ii){
			// work on each bam file one by one, rather than the whole
			if ((in = samopen(?(merged bam name), in_mode, fn_list)) == 0) {
				fprintf(stderr, "[main_samview] fail to open file for reading.\n");
				continue;
			}
			if (in->header == 0) {
				fprintf(stderr, "[main_samview] fail to read the header.\n");
				continue;
			}
			Idxs[i] = 1;
			i++;
		
			int r;
			while ((r = samread(in, b)) >= 0) { // read one alignment from `in'
				char *readgroup;
				string format;
				string alt;
				char ori = AlnParser(b, format, alt, readgroup, readgroup_platform);
				
				// need to get the line out, through samtools, don't know how now
				
				// really don't know the logic here. Why loop twice
				if(/*defined b*/ && !chr.empty() && !((b->core.chr).compare(chr)) || /*defined chr*/ || /*eof(fh*/){
					buffer[
				
				if(b->core.chr > minchr || b->core.chr == minchr && b->core.pos > minpos){// need to work on defined b
					continue;
				}
			
			// don't know why BreakDancer read each line of the bam file for nothing. # analyze only 1 chromosome
			int minchr = 255;// need to see this is the correspondence of char and int
			int minpos = 1e10;
			int minidx;
			int min_t;
			string library;
			
			// convert/print the entire file
			
			minchr = b->core.chr;
			minpos = b->core.pos;
			minidx = 
			
		}
	}
	free(fn_list); free(fn_ref); free(fn_rg);
	
 	return 0;
}			

class Keys{
	public:
		Keys(uint32_t k1, string k2) : key1(k1), key2(k2) {}
		uint32_t key1;
		string key2;
}

class Keys_int_int{
	public:
		Keys_int_int(int k1, int k2) : key1(k1), key2(k2) {}
		int key1;
		int key2;
}

void Analysis (string lib, bam1_t *b){

  //main analysis code
  //return if($t->{qual}<$opts{q} && $t->{flag}!=64 && $t->{flag}!=192);   #include unmapped reads, high false positive rate
  
	if(mapQual[lib] >= 0}{
		if(b->core.qual <= mapQual[lib])
			return;
	}
	else if(b->core.qual <= min_map_qual)
		return;
	if(b->core.chr == '*') // need to figure out how to compare a char and int //#ignore reads that failed to associate with a reference
		return;
	if(b->core.flag == 0)
		return; // return fragment reads
	if(transchr_rearrange && b->core.flag < 32 || b->core.flag >=64) // only care flag 32 for CTX
		return;
	// for long insert
	if(Illumina_long_insert){
		if(abs(b->core.isize) > uppercutoff[lib] && b->core.flag == 20))
			b->core.flag = 4;
		if(abs(b->core.isize) < uppercutoff[lib] && b->core.flag == 4))
			b->core.flag = 20;
		if(abs(b->core.isize) < lowercutoff[lib] && b->core.flag == 20))
			b->core.flag = 3;
		}
	else{
		if(abs(b->core.isize) > uppercutoff[lib] && b->core.flag == 18))
			b->core.flag = 2;
		if(abs(b->core.isize) < uppercutoff[lib] && b->core.flag == 2))
			b->core.flag = 18;
		if(abs(b->core.isize) < lowercutoff[lib] && b->core.flag == 18))
			b->core.flag = 3;
		if(b->core.flag == 20)
			b->core.flag = 4; // if it is RF orientation, then regardless of distance
	}
	if(b->core.flag == 8)
		b->core.flag = 1

	if(b->core.flag < 32 && abs(b->core.isize)>max_sd) // skip read pairs mapped too distantly on the same chromosome
		return;
	
	if(b->core.flag == 18 || b->core.flag == 20 || b->core.flag == 130){
		if(normal_switch && b->core.isize > 0)
			nnormal_reads++;
		return;
	}
	
	int do_break = (b->core.chr != lasts || b->core.pos - lastc > d)?1:0;
	
	if(do_break){ // breakpoint in the assembly
		if(lastc - beginc > min_len){ // skip short/unreliable flnaking supporting regions
			// register reliable region and supporting reads across gaps
			int k = reg_idx ++;
			sprintf(reg_name[k], "%s\t%d\t%s\t%d", begins, beginc, lastc, nnormal_reads);
			
			vector<string> p;
			for(it_reg_seq = reg_seq.begin(); it_reg_seq < reg_seq.end(); it_reg_seq ++){
				p.push_back(*it_reg_seq);
				string s = get_item_from_string(*it_reg_seq,0); // extract the ith item from the string
				read[s].push_back(k);
			}
			
			regs[k] = p;
			idx_buff++;
			if(idx_buff > buffer_size){
				buildConnection();
				idx_buff = 0;
			}
		}
		else{
			for(it_reg_seq = reg_seq.begin(); it_reg_seq < reg_seq.end(); it_reg_seq ++){
				string s = get_item_from_string(*it_reg_seq,0);
				read.erase(read.find(s));
			}
		}
		begins = b->core.chr;
		beginc = b->core.pos;
		reg_seq.clear();
		normal_switch = 0;
		nnormal_reads = 0;
	}
	bam1_qname(b) = substitue_string(); //need to figure out exactly what s/\/[12]$// is; // need to figure out the length is needed or not
	if(! prefix_fastq.empty() && ! bam1_seq(b).empty() && ! bam1_qual(b).empty()){
		string tmp_str1;
		sprintf(tmp_str1, "%s %d %d %c %d %d %s %d %s %s %s", bam1_qname(b), b->core.chr, b->core.pos, ori/*problem*/, b->core.isize, b->core.flag, b->core.qual, /* readlen */, lib, bam1_seq(b), bam1_qual(b));
		reg_seq.push_back(bam1_qname, tmp_str1);
	}
	else{
		string tmp_str1;
		sprintf(tmp_str1, "%s %d %d %c %d %d %s %d %s", bam1_qname(b), b->core.chr, b->core.pos, ori/*problem*/, b->core.isize, b->core.flag, b->core.qual, /* readlen */, lib);
	}
	if(reg_seq.size() == 0)
		normal_switch = 1;
	lasts = b->core.chr;
	lastc = b->core.pos;
}

void buildConnection(){
  // build connections
  // find paired regions that are supported by paired reads
  //warn("-- link regions\n");
	map<int, map<int, int>> link;
	map<string,vector>::iterator ii_read;
	for(ii_read = read.begin(); ii_read < read.end(); ii_read++){
		vector<int> p = *ii_read.second();
		if(p.size() != 2) // skip singleton read (non read pairs)
			continue;
		nt_int keys_int_int(p[0], p[1]);
		Keys_int_int keys_int_int_(p[1], p[0]);
		if(link.find(p[0]) < link.end() && link.find(p[0]).find(p[1]) < link[p[0]].end()){
			++link[p[0]][p[1]];
			++link[p[1]][p[0]];
		}
		else{
			link[p[0]][p[1]] = 1;
			link[p[1]][p[0]] = 1;
		}
	}
	map<int, map<int, int>> clink(link);// need to check initial of a map
	// segregate graph, find nodes that have connections
	map<> free_nodes;
	map<int, map<int, int>>::iterator ii_clink;
	for(ii_clink = clink.begin(); ii_clink < clink.end(); ii_clink++){
		s0 = *ii_clink.first();
		if(clink.find(s0) == clink.end())
			continue;
		// construct a subgraph
		vector<int> tails;
		tails.push_back(s0);
		while(tails.size() > 0){// need to figure out how to see tails is vacant since it's a class
			vector<int> newtails;
			vector<int>::iterator it_tails;
			for(it_tails = tails.begin(); it_tails < tails.end(); it_tails ++){
				if(clink.find(*it_tails) == clink.end())
					continue;
				if(reg_name.find(*it_tails) == reg_name.end())
					continue;
				vector<int> s1s;
				for(map<int, int>::iterator ii_clink_ = clink[tail].begin(); ii_clink_ < clink[tail].end(); ii_clink_++){
					s1s.push_back(*ii_clink_.first());
				}
				for(vector<int> ii_s1s = s1s.begin(); ii_s1s < s1s.end(); ii_s1s++){
					s1 = *s1s;
					vector<string> free_reads;// may not be true
					map<int,map<int,int>> nodepair;
					int nlinks = clink[tail][s1];
					if(nlinks<min_read_pair) // require sufficient number of pairs
						continue;
					if(nodepair.find(s1) < nodepair.end()) // a node only appear once in a pair
						continue;
					if(reg_name.find(s1) == reg_name.end()) // a node must be defined
						continue;
					nodepair[tail][s1] = clink[tail][s1];
					nodepair[s1][tail] = clink[s1][tail];
					clink[tail].erase(clink[tail].find(s1));	// use a link only once
					clink[s1].erase(clink[s1].find(tail));
					newtail.push_back(s1);
					
					// analysis a nodepair
					vector<int> snodes;
					for(map<int,map<int,int>>::iterator ii_nodepair = nodepair.begin(); ii_nodepair < nodepair.end(); ii_nodepair ++){
						snodes.push_back(*ii_nodepair);
					}
					int node1 = snodes[0];
					int node2;
					if(snodes.size() == 1)
						node2 = snodes[0];
					else
						node2 = snodes[1];
					if(nodepair.find(node1) == nodepair.size() || nodepair[node1].find(node2)== nodepair[node1].size())
						continue;
					
					int nread_pairs = 0;
					map<string,string> read_pair;
					map<int,int> type;
					map<int,map<string,int>> type_library_readcount;
					vector<map<int,int>> type_orient_counts;
					map<int,map<string,int>> type_library_meanspan;	//diff span distance;
					vector<string> support_reads;
					for(vector<int> ii_snodes = snodes.begin(); ii_snodes < snodes.end(); ii_snodes++){
						int node = *ii_snodes;
						map<int,int> orient_count;
						vector<string> nonsupportives;
						for(vector<string> ii_regs = regs[node].begin(); ii_regs < regs[node].end(); ii_regs++){
							string y = *ii_regs;
							if(get_item_from_string(y,0).size() == 0)
								continue;
							orient_count[get_item_from_string(y,3)]++;// where did I define and initialize this?
							
							if(read_pair[get_item_from_string(y,0)].size() == 0){
								read_pair[get_item_from_string(y,0)] = y;
								nonsupportives.push_back(y);
							}
							else{
								type[get_item_from_string(y,5)]++;// where did I initialize this
								type_library_readcount[get_item_from_string(y,5)][get_item_from_string(y,8)]++;//where I initialize this?
								type_library_meanspan[get_item_from_string(y,5)][get_item_from_string(y,8)]+=abs(get_item_from_string(y,4));
								nread_pairs++;
								free_reads.push_back(get_item_from_string(y,0));
								support_reads.push_back(y);
								support_reads.push_back(read_pair[get_item_from_string(y,0)]);
								read_pair.erase(read_pair.find(get_item_from_string(y,0)));
							}
						}
						regs[node] = nonsupportives;
						type_orient_counts.push_back(orient_count);
					}
					
					//clean out supportive reads
							
// first of all, we need to use cpp since it has better alternative (map) for the hash in perl
// second, for the rest of the code, we should use the following correspondence, since the data structure changed. Left column: perl. Right column: cpp
/*
t.readname 		:		bam1_qname(b)  (length: b.core.l_qname)
t.flag			:		b->core.flag
t.chr			:		b->core.tid (this transit from char to int)
t.pos			:		b->core.pos
t.qual			:		b->core.qual
mchr			:		b->core.mtid (this transit from char to int)
mpos			:		b->core.mpos
t.dist			:		b->core.isize
t.seq			:		bam1_seq(b)		(length: b.core.l_qseq)
t.basequal		:		bam1_qual(b)	(length: b.core.l_qseq)
t.ori			:		ori as a return in AlnParser	(in samtools, can be derived by bam1_strand(b), bam1_mstrand(b))
t.readgroup		:		readgroup as a pointer in AlnParser input
t.readlen		:		?
*/
