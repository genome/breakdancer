#include "BreakDancerMax.h"

using namespace std;

// open the bam file
int main(int argc, char *argv[])
{
	string version ("BreakDancerMax-0.0.1r81");
	int c;
	string bam_file;
	
	int chr = -1; // chromosome
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
	string prefix_fastq = "";
	string dump_BED = "";
	int Illumina_long_insert = 0;// bool
	int Illumina_to_SOLiD = 0;// bool
	while((c = getopt(argc, argv, "o:s:c:m:q:r:b:ep:tfd:g:lC")) >= 0){
		switch(c) {
			case 'o': chr = atoi(optarg); break;
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
		fprintf(stderr, "	-o INT	operate on a single chromosome [all chromosome]\n");
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
	
	string options;
	sprintf(options,"-o %d -s %d -c %d -m %d -q %d -r %d -b %d -e %d -p %d -p %f -t %d -f %d -d %s -g %s -l %d -C %d",chr, min_len, cut_sd, max_sd, min_map_qual, min_read_pair, buffer_size, learn_par, prior_prob, transchr_rearrange, fisher, prefix_fastq, dum_BED, Illumina_long_insert, Illumina_to_SOLiD);
	// define the map SVtype
	std::map<string, std::string> SVtype;
	if(Illumina_long_insert == 1){
		SVtype["1"] = "INV";
		SVtype["3"] = "INS";
		SVtype["4"] = "DEL";
		SVtype["8"] = "INV";
		SVtype["32"] = "CTX";
	}
	else{
		SVtype["1"] = "INV";
		SVtype["2"] = "DEL";
		SVtype["3"] = "INS";
		SVtype["4"] = "ITX";
		SVtype["8"] = "INV";
		SVtype["32"] = "CTX";
	}
	
	// AP; no AP
	int LZERO = -99;
	float ZERO = exp(LZERO);
	map<string,string> exes;
	map<string,string> fmaps;
	map<string,string> libmaps;
	map<string,float> mean_insertsize;//global
	map<string,float> std_insertsize;//global
	map<string,float> uppercutoff;//global
	map<string,float> lowercutoff;//global
	map<string,float> readlens;//global
	map<string,int> mapQual;// global
	int max_readlen = 0;
	map<uint32_t, map<string,int> > x_readcounts;
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
			string mean_ = get_from_line(line,"mean",0);
			string std_ = get_from_line(line,"std",0);
			string readlen_ = get_from_line(line,"readlen",0);
			string upper_ = get_from_line(line,"upp",0);
			string lower_ = get_from_line(line,"low",0);
			string mqual_ = get_from_line_two(line,"map","qual",0);
			string lib = get_from_line(line,"lib",0);
			float mean,std,readlen,upper,lower;
			int mqual;
			if(lib.compare("NA")==0)
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
			}		
			
			libmaps[lib] = fmap;
			if(mqual_.compare("NA")){
				mqual = atoi(mqual_);
				mapQual[lib] = mqual;
			}
			fmaps[fmap] = lib;
			
			if(mean.compare("NA") && std.compare("NA") && upper.compare("NA") && lower.compare("NA")){
				mean = atof(mean_);
				std = atof(stdt_);
				upper = atof(mean_) + atof(std_)*float(cut_sd);
				lower = atof(mean_) - atof(std_)*float(cut_sd);
				lower = atof(lower) > 0 ? atof(lower):0;
			}
			
			if(readlen_.compare("NA"))
				readlen = atof(readlen_);
			max_readlen = max_readlen < readlen ? readlen:max_readlen;
			
			mean_insertsize[lib] = mean;
			std_insertsize[lib] = std;
			uppercutoff[lib] = upper;
			lowercutoff[lib] = lower;
			readlens[lib] = readlen;
			
			if(exes.find(fmap) == exes.end())
				exes[fmap] = exe.compare("NA")?"cat":exe;
			else if(exes[fmap].compare(exe) != 0){
				cout << "Please use identical exe commands to open the same map file.\n";
				return;
			}
			
			int tmp = mean - readlen*2;	// this determines the mean of the max of the SV flanking region
			d = d<tmp ? d:tmp;
		}
	}

	CONFIG.close();
		
	if(d < 50)
		d = 50;
	
	if(dump_BED.empty())
		ofstream BED(dump_BED);
	vector<string> format; 
	vector<string>::iterator it_format;

	map<string,int> cmds; 	
 	
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
	 	EstimatePriorParameters(fmaps, readgroup_library, mean_insertsize, std_insertsize, uppercutoff, lowercutoff, readlens, chr);
	 	
 	int reference_len = 1;
 	map<string, int> nreads;
 	int defined_all_readgroups = 1;
 	    
	samfile_t *in = 0;
	char in_mode[5] = "", *fn_list = 0, *fn_ref = 0, *fn_rg = 0;
	strcpy(in_mode, "r");
	int i = 0;
	bam1_t *b = bam_init1();
	string format_ = "sam";
	string alt = "";	
 	for(ii=fmaps.begin(); ii!=fmaps.end(); ++ii)
 	{
 		int ref_len = 0;
 		string exe = exes[maps[i]];
 		i ++;
 		cmds[exe] ++;
 		
 		int p_pos = 0;
 		int p_chr;
 		
		/*// convert the entire bam file
		if ((in = samopen((*ii).first, in_mode, fn_list)) == 0) {
			fprintf(stderr, "[main_samview] fail to open file for reading.\n");
			continue;
		}
		if (in->header == 0) {
			fprintf(stderr, "[main_samview] fail to read the header.\n");
			continue;
		}
		bam1_t *b = bam_init1();
		int r;
		while ((r = samread(in, b)) >= 0) { */
		
		// read the bam file by a chromosome
		string bam_name = *ii.first;
		bam_index_t *idx;
		samfile_t *in;
		int tid, beg, end, n_off;
		string chr_str;
		sprintf(chr_str, "%d", chr);
		pair64_t n_off;
		bamFile fp = ReadBamChr_prep(chr_str, bam_name, &tid, &beg, &end, in, off, &n_off);
		uint64_t curr_off;
		int i, ret, n_seeks;
		n_seeks = 0; i = -1; curr_off = 0;
		 
		while(ReadBamChr(b, fp, tid, beg, end, &curr_off, &i, &n_seeks, off, n_off)){
			char *readgroup;

			char ori = AlnParser(b, format_, alt, readgroup, readgroup_platform);
		
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
				
			if(nreads.find(lib) == nreads.end())
				nreads[lib] = 1;
			else
				nreads[lib] ++;	
			if(mapQual[lib] != 0 && b->core.qual <= mapQual[lib])
				continue;
			else if(b->core.qual <= min_map_qual)
				continue;
			if(b->core.flag == 0)
				continue;
			if(transchr_rearrange && b->core.flag < 32 || b->core.flag >= 64)
				continue;
				
			if(Illumina_long_insert){
				if(abs(b->core.isize) > uppercutoff[lib] && b->core.flag == 20)
					b->core.flag = 4;
				if(abs(b->core.isize) < uppercutoff[lib] && b->core.flag == 4)
					b->core.flag = 20;
				if(abs(b->core.isize) < lowercutoff[lib] && b->core.flag == 20)
					b->core.flag = 3;
			}
			else{
				if(abs(b->core.isize) > uppercutoff[lib] && b->core.flag == 18)
					b->core.flag = 2;
				if(abs(b->core.isize) < uppercutoff[lib] && b->core.flag == 2)
					b->core.flag = 18;
				if(abs(b->core.isize) < lowercutoff[lib] && b->core.flag == 18)
					b->core.flag = 3;
			}
			
			if(b->core.flag == 18 || b->core.flag == 20 || b->core.flag == 130)
				continue;
			
			if(x_readcounts.find(b->core.flag) < x_readcounts.end() && x_readcounts[b->core.flag].find(lib) < x_readcounts[b->core.flag].end())
				x_readcounts[b->core.flag][lib] ++;	
			else
				x_readcounts[b->core.flag][lib] = 1;					 
		}
		
		if (r < -1) fprintf(stderr, "[main_samview] truncated file.\n");
		if(ref_len == 0)
			fprintf(stderr, "Unable to decode %s. Please check that you have the correct paths and the bam files are indexed.", *ii);
		if(reference_len < ref_len)
			reference_len = ref_len;
		samclose(in);		
	}
	bam_destroy1(b);
	int merge = (cmds.size()==1 && defined_all_readgroups)?1:0;
	
	float total_phy_cov = 0;
	float total_seq_cov = 0;
	
	// map has already been sorted by the key
	// recflags = x_readcounts[keys the first key]
	map<string,int>::iterator nreads_ii;
	for(nreads_ii=nreads.begin(); nreads_ii!=nreads.end(); ++nreads_ii)
	{
		string lib = *nreads_ii.first;
		float sequence_coverage = float(nreads[lib]*readlen[lib])/float(reference_len);
		total_seq_cov += sequence_coverage;
		float physical_coverage = float(nreads[lib]*mean_insertsize[lib])/float(2*reference_len);
		total_phy_cov += physical_coverage;
		
		int nread_lengthDiscrepant = -1;
		
		if(x_readcounts[2][lib])
			nread_lengthDiscrepant = x_readcounts[2][lib];
		if(x_readcounts[3][lib]){
			if(nread_lengthDiscrepant == -1)
				nread_lengthDiscrepant = 0;
			nread_lengthDiscrepant += x_readcounts[3][lib];
		}
		float tmp = (nread_lengthDiscrepant > 0)?(float)reference_len/(float)nread_lengthDiscrepant:50;
		d = d<tmp?d:tmp;
		
		printf("#%s\tmean:%.3f\tstd:%.3f\tuppercutoff:%.3f\tlowercutoff:%.3f\treadlen:%.3f\tlibrary:%s\treflen:%d\tseqcov:%.3fx\tphycov:%.3fx", libmaps[$lib],mean_insertsize[lib],std_insertsize[lib],uppercutoff[lib],lowercutoff[lib],readlens[lib],lib,reference_len, sequence_coverage,physical_coverage);
		
		map<uint32_t,map<string,int> >::iterator x_readcounts_ii;
		for(x_readcounts_ii = x_readcounts.begin(); x_readcounts_ii!=x_readcounts.end(); ++x_readcounts_ii){
			uint32_t t = *x_readcounts_ii.first;// get the first key out, which is a member of recflags
			if(*x_readcounts_ii.second.find(lib) == *x_readcounts_ii.second.end())
				printf("\t%d:0",t);
			else
				printf("\t%d:%d",t,*x_readcounts_ii.second[lib]);
		}
		printf("\n");
	}
	printf("#Chr1\tPos1\tOrientation1\tChr2\tPos2\tOrientation2\tType\tSize\tScore\tnum_Reads\tnum_Reads_lib\tAllele_frequency\tVersion\tRun_Param\n");
	
	int begins;// global (chr)
	int beginc = -1;// global
	int lasts;// global (chr, should be int in samtools)
	int lastc = -1; // global
	map<int, vector<vector<string> > > regs;//global in analysis
	map<string, vector<int> > read;// global in analysis
	map<int,vector<int> > reg_name;// global in analysis
	vector<vector<string> > reg_seq; // global need to see if it's the key or value of one of the above global. should be a string
	vector<vector<string> >::iterator it_reg_seq; // global
	
	int idx_buff = 0;// global
	int final_buff = 0;
	int reg_idx = 0;// global  ################# node index here ###################
	int normal_switch = 0; // global
	int nnormal_reads = 0; // global
	
	// first, need to merge the bam files into one big string seperated by blank, and return the number
	int n = fmaps.size();
	string *big_bam = new string[bam_number];
	int i_tmp = 0;
	for(map<string, string> ii_fmaps = fmaps.begin(); ii_fmaps < fmaps.end(); ii_fmaps++)
		big_bam[i_tmp++] = *ii_fmaps.first;
	
 	bamFile *fp;
	heap1_t *heap;
	fp = (bamFile*)calloc(n,sizeof(bamFile));
	heap = (heap1_t*)calloc(n,sizeof(heap1_t));
	//################# node index here ###################
	if(merge && fmaps.size()>1 && !(format[0].compare("sam")) && chr == -1 ){
 	/* open pipe, improvement made by Ben Oberkfell (boberkfe@genome.wustl.edu)
   samtools merge - in1.bam in2.bam in3.bam in_N.bam | samtools view - 
   maq mapmerge		*/
   		
   		// dig into merge samtools code and utilize what we needed
   		if(MergeBams_prep(big_bam, n, fp, heap)){
   		
	   		while(heap->pos != HEAP_EMPTY){
   				bam1_t *b = heap->b;
   			
   				char *readgroup;
				char ori = AlnParser(b, format_, alt, readgroup, readgroup_platform);
				string library = readgroup.empty()?readgroup_library[readgroup]:(*fmaps.begin().second);
			  	//if(chr.empty() || chr.compare(b->core.tid)!=0) //this statement actually does nothing
			  		//continue;
			  	if(!library.empty())
			  		Analysis (lib, b, reg_seq, reg_name, read, regs, begins, beginc, lasts, lastc, idx_buff, buffer_size, nnormal_reads, min_len, normal_switch, reg_idx, transchr_rearrange, min_map_qual, Illumina_long_insert, prefix_fastq, x_readcounts, reference_len, fisher);
		  		
   				if((j = bam_read1(fp[heap->i], b)) >= 0){
   					heap -> pos = ((uint64_t)b->core.tid<<32) | (uint32_t)b->core.pos << 1 | bam1_strand(b);
   					heap -> idx = idx++;
   				}
   				else if(j == -1){
   					heap->pos = HEAP_EMPTY;
   					free(heap->b->data);
   					free(heap->b);
   					heap->b = 0;
   				}
   				else
   					cout << "[bam_merge_core] " << big_bam[heap->i] << " is truncated. Continue anyway.\n";
   				ks_heapadjust(heap, 0, n, heap);
   			}
		}
	}
	//############### find the designated chromosome and put all of them together for all bam files #############
	else{//#customized merge & sort
	// totally different from the perl version; use customized samtools instead
		samfile_t **in;
		*in = new samfile_t[n];
		
		int *tid, *beg, *end, *n_off;
		tid = new int[n];
		beg = new int[n];
		end = new int[n];
		n_off = new int[n];
		
		string chr_str;
		sprintf(chr_str, "%d", chr);
		
		pair64_t **off;
		*off = new pair64_t[n];
		
		uint64_t *curr_off;
		curr_off = new uint64_t[n];
		int *i, *n_seeks;
		i = new int[n];
		n_seeks = new int[n];
		for(int k = 0; k < n; k++){
			n_seeks[k] = 0;
			i[k] = -1;
			curr_off[k] = 0;
		}
		
		if(MergeBamsChr_prep(string *fn, int n, bamFile *fp, heap1_t *heap, string chr_str, int *tid, int *beg, int *end, samfile_t **in, pair64_t **off, int *n_off)){
			while(heap->pos != HEAP_EMPTY){
   				bam1_t *b = heap->b;
   			
   				char *readgroup;
				char ori = AlnParser(b, format_, alt, readgroup, readgroup_platform);
				string library = readgroup.empty()?readgroup_library[readgroup]:(*fmaps.begin().second);
		  		//if(chr.empty() || chr.compare(b->core.tid)!=0) //this statement actually does nothing
		  			//continue;
		  		if(!library.empty())
		  			Analysis (lib, b, reg_seq, reg_name, read, regs, begins, beginc, lasts, lastc, idx_buff, buffer_size, nnormal_reads, min_len, normal_switch, reg_idx, transchr_rearrange, min_map_qual, Illumina_long_insert, prefix_fastq, x_readcounts, reference_len, fisher);
		  		
   				if(ReadBamChr(b, fp[heap->i], tid[heap->i], beg[heap->i], end[heap->i], curr_off[heap->i], i[heap->i], n_seeks[heap->i], off[heap->i], n_off[heap->i])>0){
   					heap -> pos = ((uint64_t)b->core.tid<<32) | (uint32_t)b->core.pos << 1 | bam1_strand(b);
   					heap -> idx = idx++;
   				}
   				else if(j == -1){
   					heap->pos = HEAP_EMPTY;
   					free(heap->b->data);
   					free(heap->b);
   					heap->b = 0;
   				}
   				else
   					cout << "[bam_merge_core] " << big_bam[heap->i] << " is truncated. Continue anyway.\n";
   				ks_heapadjust(heap, 0, n, heap);
   			}
		}
		
		delete []tid, []beg, []end, []n_off, []curr_off, []i, []n_seeks;
		for(k = 0; k < n; k++)
			samclose(in[k]);
	}
	
	for(k = 0; k!=n; k++)
   	bam_close(fp[k]);
   	free(fp);
   	free(heap);
	
	delete []big_bam;
	big_bam = NULL;
	free(fn_list); free(fn_ref); free(fn_rg);
	
 	return 0;
}

// this function is good
void Analysis (string lib, bam1_t *b, vector<vector<string> > &reg_seq, vector<int,vector<int> > &reg_name, map<string,vector<int> > &read, map<int, vector<vector<string> > > &regs, int &begins, int &beginc, int &lasts, int &lastc, int &idx_buff, int buffer_size, int &nnormal_reads, int min_len, int &normal_switch, int &reg_idx, int transchr_rearrange, int min_map_qual, int Illumina_long_insert, int prefix_fastq, map<uint32_t, map<string,int> > &x_readcounts, int reference_len, int fisher){

  //main analysis code
  //return if($t->{qual}<$opts{q} && $t->{flag}!=64 && $t->{flag}!=192);   #include unmapped reads, high false positive rate
  
	if(mapQual[lib] >= 0){
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
		if(abs(b->core.isize) > uppercutoff[lib] && b->core.flag == 20)
			b->core.flag = 4;
		if(abs(b->core.isize) < uppercutoff[lib] && b->core.flag == 4)
			b->core.flag = 20;
		if(abs(b->core.isize) < lowercutoff[lib] && b->core.flag == 20)
			b->core.flag = 3;
		}
	else{
		if(abs(b->core.isize) > uppercutoff[lib] && b->core.flag == 18)
			b->core.flag = 2;
		if(abs(b->core.isize) < uppercutoff[lib] && b->core.flag == 2)
			b->core.flag = 18;
		if(abs(b->core.isize) < lowercutoff[lib] && b->core.flag == 18)
			b->core.flag = 3;
		if(b->core.flag == 20)
			b->core.flag = 4; // if it is RF orientation, then regardless of distance
	}
	if(b->core.flag == 8)
		b->core.flag = 1;

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
			//sprintf(reg_name[k], "%s\t%d\t%d\t%d", begins, beginc, lastc, nnormal_reads);
			reg_name[k].push_back(begins);
			reg_name[k].push_back(beginc);
			reg_name[k].push_back(lastc);
			reg_name[k].push_back(nnormal_reads);
			
			vector<vector<string>> p;
			for(vector<vector<string>> it_reg_seq = reg_seq.begin(); it_reg_seq < reg_seq.end(); it_reg_seq ++){
				p.push_back(*it_reg_seq);
				///string s = get_item_from_string(*it_reg_seq,0); // extract the ith item from the string
				string s = (*it_reg_seq)[0];
				read[s].push_back(k);
			}
			
			regs[k] = p;
			idx_buff++;
			if(idx_buff > buffer_size){
				buildConnection(read, reg_name, regs, x_readcounts, reference_len, fisher);
				idx_buff = 0;
			}
		}
		else{
			for(vecotr<vector<string>> it_reg_seq = reg_seq.begin(); it_reg_seq < reg_seq.end(); it_reg_seq ++){
				///string s = get_item_from_string(*it_reg_seq,0);
				string s= (*it_reg_seq)[0];
				read.erase(read.find(s));
			}
		}
		begins = b->core.chr;
		beginc = b->core.pos;
		reg_seq.clear();
		normal_switch = 0;
		nnormal_reads = 0;
	}
	// deal with the name string
	string qname_tmp = bam1_qname(b);
	size_t found1 = qname_tmp.rfind("/1");
	size_t found2 = qname_tmp.rfind("/2");
	if(found1 != string::npos || found2 != string::npos){
		size_t found = (found1 == string::npos) ? found2 : found1;
		qname_tmp.replace(found,2,"");
	}
	bam1_qname(b) = qname_tmp; 
	
	if(! prefix_fastq.empty() && ! bam1_seq(b).empty() && ! bam1_qual(b).empty()){
		//string tmp_str1;
		//sprintf(tmp_str1, "%s %d %d %c %d %d %s %d %s %s %s", bam1_qname(b), b->core.chr, b->core.pos, ori/*problem*/, b->core.isize, b->core.flag, b->core.qual, /* readlen */, lib, bam1_seq(b), bam1_qual(b));
		//reg_seq.push_back(bam1_qname, tmp_str1);
		vector<string> tmp_reg_seq;
		tmp_reg_seq.push_back(bam1_qname(b));
		tmp_reg_seq.push_back(itoa(b->core.chr));
		tmp_reg_seq.push_back(itoa(b->core.pos));
		tmp_reg_seq.push_back(ori);
		tmp_reg_seq.push_back(itoa(b->core.isize));
		tmp_reg_seq.push_back(itoa(b->core.flag));
		tmp_reg_seq.push_back(b->core.qual);
		tmp_reg_seq.push_back(itoa(b->core.l_qseq));
		tmp_reg_seq.push_back(lib);
		tmp_reg_seq.push_back(bam1_seq(b));
		tmp_reg_seq.push_back(bam1_qual(b));
		reg_seq.push_back(tmp_reg_seq);
	}
	else{
		//string tmp_str1;
		//sprintf(tmp_str1, "%s %d %d %c %d %d %s %d %s", bam1_qname(b), b->core.chr, b->core.pos, ori/*problem*/, b->core.isize, b->core.flag, b->core.qual, /* readlen */, lib);
		vector<string> tmp_reg_seq;
		tmp_reg_seq.push_back(bam1_qname(b));
		tmp_reg_seq.push_back(itoa(b->core.chr));
		tmp_reg_seq.push_back(itoa(b->core.pos));
		tmp_reg_seq.push_back(ori);
		tmp_reg_seq.push_back(itoa(b->core.isize));
		tmp_reg_seq.push_back(itoa(b->core.flag));
		tmp_reg_seq.push_back(b->core.qual);
		tmp_reg_seq.push_back(itoa(b->core.l_qseq));
		tmp_reg_seq.push_back(lib);
		reg_seq.push_back(tmp_reg_seq);		
	}
	if(reg_seq.size() == 0)
		normal_switch = 1;
	lasts = b->core.chr;
	lastc = b->core.pos;
}

// this function is good. May miss the i/o
void buildConnection(map<string,vector<int> > &read, map<int,vector<int> > &reg_name, map<int,vector<vector<string> > > &regs, map<uint32_t, map<string,int> > &x_readcounts, int reference_len, int fisher){
  // build connections
  // find paired regions that are supported by paired reads
  //warn("-- link regions\n");
	map<int, map<int, int> > link;
	map<string,vector<int> >::iterator ii_read;
	for(ii_read = read.begin(); ii_read < read.end(); ii_read++){
		vector<int> p = *ii_read.second;
		if(p.size() != 2) // skip singleton read (non read pairs)
			continue;
		if(link.find(p[0]) < link.end() && link[p[0]].find(p[1]) < link[p[0]].end()){
			++link[p[0]][p[1]];
			++link[p[1]][p[0]];
		}
		else{
			link[p[0]][p[1]] = 1;
			link[p[1]][p[0]] = 1;
		}
	}
	map<int, map<int, int> > clink(link);
	// segregate graph, find nodes that have connections
	map<int,int> free_nodes;
	map<int, map<int, int> >::iterator ii_clink;
	for(ii_clink = clink.begin(); ii_clink < clink.end(); ii_clink++){
		int s0 = *ii_clink.first;
		if(clink.find(s0) == clink.end())
			continue;
		// construct a subgraph
		vector<int> tails;
		tails.push_back(s0);
		while(tails.size() > 0){
			vector<int> newtails;
			vector<int>::iterator it_tails;
			for(it_tails = tails.begin(); it_tails < tails.end(); it_tails ++){
				if(clink.find(*it_tails) == clink.end())
					continue;
				if(reg_name.find(*it_tails) == reg_name.end())
					continue;
				vector<int> s1s;
				for(map<int, int>::iterator ii_clink_ = clink[tail].begin(); ii_clink_ < clink[tail].end(); ii_clink_++){
					s1s.push_back(*ii_clink_.first);
				}
				for(vector<int> ii_s1s = s1s.begin(); ii_s1s < s1s.end(); ii_s1s++){
					s1 = *s1s;
					vector<string> free_reads;
					map<int,map<int,int> > nodepair;
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
					for(map<int,map<int,int> >::iterator ii_nodepair = nodepair.begin(); ii_nodepair < nodepair.end(); ii_nodepair ++){
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
					map<string,vector<string>> read_pair;
					map<uint32_t,int> type;//first one is flags
					map<uint32_t,map<string,int> > type_library_readcount;
					vector<map<char,int> > type_orient_counts;
					map<uint32_t,map<string,int> > type_library_meanspan;	//diff span distance;
					vector<vector<string>> support_reads;
					for(vector<int> ii_snodes = snodes.begin(); ii_snodes < snodes.end(); ii_snodes++){
						int node = *ii_snodes;
						map<char,int> orient_count;
						vector<string> nonsupportives;
						for(vector<vecotr<string>> ii_regs = regs[read_pairnode].begin(); ii_regs < regs[node].end(); ii_regs++){
							vector<string> y = *ii_regs;
							if(y[0].size() == 0)
								continue;
							orient_count[y[3]]++;// where did I define and initialize this?
							
							if(read_pair[y[0]].size() == 0){
								read_pair[y[0]] = y;
								nonsupportives.push_back(yrea);
							}
							else{
								// see if initialized 'type' or not
								if(type.find(y[5]) != type.end())
									type[y[5]]++;
								else
									type[y[5]] = 1;
								// see if initialized 'type_library_readcount' or not
								if(type_library_readcount.find(y[5]) != type.end() && type_library_readcount[y[5]].find(y[8]) != type[y[5]].end())
									type_library_readcount[y[5]][y[8]]++;
								else{
									map<string, int> tmp_type_library_readcount;
									tmp_type_library_readcount[y[8]] = 1;
									type_library_readcount[y[5]] = tmp_type_library_readcount;
								}
								// see if initialized 'type_library_meanspan' or not
								if(type_library_meanspan.find(y[5]) != type_library_meanspan.end() && type_library_meanspan[y[5]].find(y[8]) != type_library_meanspan[y[5]].end())
									type_library_meanspan[y[5]][y[8]]+=abs(atoi(y[4]));
								else{
									map<string,int> tmp_type_library_meanspan;
									tmp_type_library_meanspan[y[8]] = 1;
									type_library_meanspan[y[5]] = tmp_type_library_meanspan;
								}
								nread_pairs++;
								free_reads.push_back(y[0]);
								support_reads.push_back(y);
								support_reads.push_back(read_pair[y[0]]);
								read_pair.erase(read_pair.find(y[0]));
							}
						}
						regs[node] = nonsupportives;
						type_orient_counts.push_back(orient_count);
					}
					
					//clean out supportive reads
					for(vector<int> ii_snodes = snodes.begin(); ii_snodes < snodes.end(); ii_snodes++){
						int node = *ii_snodes;
						for(vector<vector<string>> ii_regs = regs[node].begin(); ii_regs < regs[node].end(); ii_regs++){
							vector<string> y = *ii_regs;
							if(read_pair.find(y[0]) = read_pair.end())
								continue;
							nonsupportive.push_back(y);	
						}
					}
					
					//float score;//don't know if float; no usage actually
					//int bestIndelSize;//don't know if int; no usage actually
					if(nread_pairs >= min_read_pair){
						float maxscore = 0;
						uint32_t flag = 0;
						map<string, int> diffspans;
						map<string, string> sptypes;
						for(map<string,int> ii_type = type.begin(); ii_type < type.end(); ii_type ++){
							uint32_t fl = *ii_type.first;
							float ptype = float(*ii_type.second())/float(nread_pairs);
							if(maxscore<ptype){
								maxscore = ptype;
								flag = fl;
							}
							string sptype;
							float diffspan = 0;
							for(map<string,int> ii_type_lib_rc = type_library_readcount[fl].begin(); ii_type_lib_rc < type_library_readcount[fl].end(); ii_type_lib_rc ++){
								string sp = *ii_type_lib_rc.first;
								string str_num_tmp;
								sprintf(str_num_tmp, "%s", *ii_type_lib_rc.second); 
								if(!sptype.empty())
									sptype.append(":").append(sp).append("|").append(str_num_tmp);
								else
									sptype = sp.append("|").append(str_num_tmp);
								diffspan += float(type_library_meanspan[fl][sp]) - float(type_library_readcount[fl][sp])*mean_insertsize[sp];
							}
							diffspans[fl] = int(diffspan/float(type[fl]) + 0.5);
							sptypes[fl] = sptype;
						}
						
						if(type[flag] >= min_read_pair){
							// print out result
							int sv_chr1 = 0, sv_pos1 = 0, sv_chr2 = 0, sv_pos2 = 0;
							string sv_ori1, sv_ori2;
							int normal_rp; 
							// find inner most positions
							for(vector<int> ii_snodes = snodes.begin(); ii_snodes < snodes.end(); ii_snodes ++){
								int node = *ii_snodes;
								int chr = reg_name[node][0];
								int start = reg_name[node][1];
								int end = reg_name[node][2];
								int nrp = reg_name[node][3];
								map<char,int> ori_readcount = type_orient_counts.front();
								type_orient_counts.erase(type_orient_counts.begin());
								if(sv_chr1 != 0 && sv_chr2 != 0){
									sv_chr1 = sv_chr2;
									sv_pos1 = sv_pos2;
									sv_chr2 = chr;
									sv_pos2 = start;
									string sv_ori2_tmp1 = "0";
									string sv_ori2_tmp2 = "0";
									if(ori_readcount.find('+') < ori_readcount.end())
										sprintf(sv_ori2_tmp1, "%s", ori_readcount['+']);
									if(ori_readcount_find('-') < ori_readcount.end())
										sprintf(sv_ori2_tmp2, "%s", ori_readcount['-']);
									sv_ori2 = sv_ori2_tmp.append("+").append(sv_ori2_tmp2).append("-");
								}
								else{
									sv_chr1 = chr;
									sv_chr2 = chr;
									sv_pos1 = start;
									sv_pos2 = end;
									string sv_ori2_tmp1 = "0";
									string sv_ori2_tmp2 = "0";
									if(ori_readcount.find('+') < ori_readcount.end())
										sprintf(sv_ori2_tmp1, "%s", ori_readcount['+']);
									if(ori_readcount_find('-') < ori_readcount.end())
										sprintf(sv_ori2_tmp2, "%s", ori_readcount['-']);
									sv_ori1 = sv_ori2_tmp.append("+").append(sv_ori2_tmp2).append("-");
									sv_ori2 = sv_ori1;
									normal_rp = nrp;
								}
							}
							
							float LogPvalue = ComputeProbScore(snode, type_library_readcount[flag], flag, x_readcounts, reference_len, fisher);//need to work on the function;
							float PhredQ_tmp = -10*LogPvalue/log(10);
							int PhredQ = PhredQ_tmp>99 ? 99:int(PhredQ_tmp+0.5);
							
							float AF = float(type[flag])/float(type[flag]+normal_rp);
							
							sv_pos1 += max_readlen - 5; // apply extra padding to the start coordinates
							
							string SVT = SVtype.find(flag)==SVtype.end()?"UN":SVtype[flag]; // UN stands for unknown
							printf("%s\t%d\t%s\t%s\t%d\t%s\t%s\t%d\t%d\t%d\t%s\t%.2f\t%s\t%s\n",sv_chr1,sv_pos1,sv_ori1,sv_chr2,sv_pos2,sv_ori2,SVT,diffspans[flag],PhredQ,type[flag],sptypes[flag],AF,version,options);// version and options should be figured out
							
							if(!prefix_fastq.empty()){ // print out supporting read pairs
								map<string,int> pairing;
								for(vector<vector<string>> ii_support_reads = support_read.begin(); ii_support_read < support_read.end(); ii_support_read ++){
									vector<string> y = *ii_support_reads;
									if(y.size() != 10 || y[5].compare(flag) != 0)
										continue;
									string fh = (pairing.find(y[0])< pairing.end())?ReadsOut[y[8]].append("1"):ReadsOut[y[8]].append("2");
									pairing[y[0]] = 1;
									sprintf(fh,"@%s\n",y[0]);
									sprintf(fh,"%s\n",y[9]);
									sprintf(fh,"+\n%s\n",y[10]);
								}
							}
							
							if(!dump_BED.empty()){	// print out SV and supporting reads in BED format
								ofstream fh_BED;
								fh_BED.open(BED);
								string trackname = itoa(sv_chr1).append(itoa(sv_pos1)).append(SVT).append(itoa(diffspans[flag]));
								fh_BED << "track name=" << trackname << "\tdescription=\"BreakDancer " << sv_chr1 << " " << sv_pos1 << " " << SVT << " " << diffspans[flag] << "\"\tuseScore=0\n"; //"track name=%s	description=\"BreakDancer %d %d %s %d\"	useScore=0\n", trackname, sv_chr1, sv_pos1, SVT, diffspans[flag]);// fh_BED is a file handle of BED
								for(vector<vector<string>> ii_support_reads = support_read.begin(); ii_support_read < support_read.end(); ii_support_read ++){
									vector<string> y = *ii_support_reads;
									if(y.size() < 8 || y[5].compare(flag) != 0)
										continue;
									if(y[1].find(chr)!=string::npos)
										y[1] = "";
									int aln_end = atoi(y[2]) + atoi(y[7]) - 1;
									string color = y[3].compare("+")?"255,0,0":"0,0,255";
									fh_BED << "chr" << y[1] << "\t" << y[2] << "\t" << aln_end << "\t1\t" << y[0] << "\t" << y[3] << "\t" << y[4] << "\t" << y[2] << "\t" << aln_end << "\t" << color;//sprintf(fh_BED, "chr%s\t%s\t%s\t%s\t1\t%s\t%s\t%s\t%d\t%s\n",y[1],y[2],aln_end,y[0],y[3],y[4],y[2],aln_end,color);
								}
								fh_BED.close();
							}
						}
						// free reads
						free_reads.clear();
						//record list of nodes that can be potentially freed
						free_nodes[node1] = 1;
						free_nodes[node2] = 1;
					}
				}
				clink.erase(clink.find(tail));
			}
			tails = newtails;
		}
	}
	
	// free nodes
	for(map<int,int> ii_free_nodes = free_nodes.begin(); ii_free_nodes < free_nodes.end(); ii_free_nodes++){
		// remove reads in the regions
		int node = *ii_free_nodes.first;
		vector<vector<string>> reads = regs[node];
		if(reads.size()+1 < min_read_pair){
			for(vector<vector<string>> ii_reads = reads.begin(); ii_reads < reads.end(); ii_reads++){
				vector<string> y = *ii_reads;
				string readname = y[0];
				read.erase(read.find(readname));
			}
			// remove regions
			regs.erase(regs.find(node));
			reg_name.erase(reg_name.find(node));
		}
	}
}

// this function is good
void EstimatePriorParameters(map<string,string> &fmaps, map<string,string> &readgroup_library, map<string, float> &mean_insertsize, map<string, float> &std_insertsize, map<string,float> &uppercutoff, map<string,float> &lowercutoff, map<string,float> &readlens, int chr){
	map<string,float> es_means;
	map<string,float> es_stds;
	map<string,float> es_readlens;
	map<string,float> es_uppercutoff;
	map<string,float> es_lowercutoff;
	map<string,vector<int> > insert_stat;
	map<string,vector<int> > readlen_stat;
	
	bam1_t *b = bam_init1();
	for(map<string,string> ii=fmaps.begin(); ii!=fmaps.end(); ++ii){
		// read the bam file by a chromosome
		string bam_name = *ii.first;
		bam_index_t *idx;
		samfile_t *in;
		int tid, beg, end, n_off;
		string chr_str;
		sprintf(chr_str, "%d", chr);
		pair64_t *off;
		bamFile fp = ReadBamChr_prep(chr_str, bam_name, &tid, &beg, &end, in, off, &n_off);
		uint64_t curr_off;
		int i, ret, n_seeks;
		n_seeks = 0; i = -1; curr_off = 0;
		while(ReadBamChr(b, fp, tid, beg, end, &curr_off, &i, &n_seeks, off, n_off)){
			char *readgroup;
			string format = "sam";
			string alt = "";
			char ori = AlnParser(b, format, alt, readgroup, readgroupreadG_platform)
		
			// analyze the bam file line by line
			string lib = readgroup.empty()?*ii.second():readgroup_library[readgroup];// when multiple libraries are in a BAM file
			if(lib.empty())
				continue;
			if(!chr.empty() && b->core.tid != chr)// analyze 1 chromosome
				continue;
			//if(readlen_stat.find(lib) == readlen_stat.end())	// don't need to issue a new stat
				//readlen_stat[lib] = ; // Statistics::Descriptive::Sparse->new() // don't need to issue a new stat
			readlen_stat[lib].push_back(b->core.isize);
			if(b->core.qual <= min_map_qual)	// skip low quality mapped reads
				continue;
			if(b->core.flag != 18 && b->core.flag != 20 || b->core.isize <= 0)
				continue;
			//if(insert_stat.find(lib) == insert_stat.end())	// don't need to issue a new stat
			insert_stat[lib].push_back(b->core.isize);
		}
		samclose(in);		
	}
	bam_destroy1(b);
	for(map<string,vector<int> > ii_readlen_stat = readlen_stat.begin(); ii_readlen_stat < readlen_stat.end(); ii_readlen_stat ++){
		//double res = accumulate(insert_stat[lib].begin(), insert_stat[lib].end(), 0);
		//double mean = res/insert_stat[lib].size();
		float mean_insert = mean(insert_stat[lib]);
		float std_insert = standard_deviation(insert_stat[lib],mean_insert);
		float uppercutoff = mean_insert + std_insert*cut_sd;
		float lowercutoff = mean_insert - std_insert*cut_sd;
		es_readlens[lib] = mean(readlen_stat[lib]);
		es_means[lib] = mean_insert;
		es_stds[lib] = std_insert;
		es_uppercutoff[lib] = uppercutoff;
		es_lowercutoff[lib] = lowercutoff;
	}
	mean_insertsize = es_means;
	std_insertsize = es_stds;
	uppercutoff = es_uppercutoff;
	lowercutoff = es_lowercutoff;
	readlens = es_readlens;
	return;
}

// compute the mean of a vector of int
float mean(vector<int> &stat){
	int all = 0;
	vector<int> ii_stat;
	for(ii_stat = stat.begin(); ii_stat < stat.end(); ii_stat++){
		all += *ii_stat;
	}
	return (float)all/(float)(stat.size());
}

// compute the standard deviation of a vector of int knowing the mean
float standard_deviation(vector<int> &stat, float mean){
	int all = 0;
	vector<int> ii_stat;
	for(ii_stat = stat.begin(); ii_stat < stat.end(); ii_stat ++){
		all += (*ii_stat)*(*ii_stat);
	}
	return sqrt((float)all/(float)(stat.size()) - mean*mean);
}	

float ComputeProbScore(vector<int> rnode, map<string,int> rlibrary_readcount, uint32_t type, map<uint32_t, map<string,int> > &x_readcounts, int reference_len, int fisher) {
	// rnode, rlibrary_readcount, type
	int total_region_size = PutativeRegion(rnode);
	
	float lambda;
	float logpvalue = 0;
	for(map<int,map<string,int> > ii_rlibrary_readcount = rlibrary_readcount.begin(); ii_rlibrary_readcount < rlibrary_readcount.end(); ii_rlibrary_readcount ++){
		string lib = *ii_rlibrary_readcount.first;
		lambda = float(total_region_size* x_readcounts[type][lib])/float(reference_len);
		logpvalue += LogPoissonTailProb(float((rlibrary_readcount[lib]).size()),lambda);
	}
	
	if(fisher && logpvalue < 0){
		// Fisher's Method
		float fisherP = 1 - // Math::CDF::pchisq(-2*logpvalue, 2*(rlibrary_readcount.size()+1));
		logpvalue = (fisherP > ZERO)?log(fisherP):LZERO;
	}
	return logpvalue;
}

// this function is good
int PutativeRegion(vector<int> rnode, map<int,vector<int> > &reg_name){
	int total_region_size = 0;
	for(vector<int> ii_node = rnode.begin(); ii_node < rnode.end(); ii_node++){
		int node = *ii_node;
		int clust_start = reg_name[node][1];
		int clust_end = reg_name[node][2];
		total_region_size += clust_end - clust_start + 1;
	}
	return total_region_size;
}

// prepare: read bam file by a chromosome
bamFile ReadBamChr_prep(string chr_str, string bam_name, int &tid, int &beg, int &end, samfile_t *in, pair64_t *off, int &n_off){
	char *in_mode, *fn_list = 0;
	strcpy(in_mode, "r");
	if ((in = samopen(bam_name, in_mode, fn_list)) == 0) {
		fprintf(stderr, "[main_samview] fail to open file for reading.\n");
		continue;
	}
	if (in->header == 0) {
		fprintf(stderr, "[main_samview] fail to read the header.\n");
		continue;
	}
	bam_index_t *idx;
	idx = bam_index_load(bam_name);// index
	bam_parse_region(in->header, chr_str, &tid, &beg, &end);// parse
	
	// return the file handle for handle
	bamFile fp = in->x.bam;
	off = get_chunk_coordinates(idx, tid, beg, end, &n_off);
	return fp;
}

// read bam file by a chromosome by one line; fp will track where we are
int ReadBamChr(bam1_t *b, bamFile fp, int tid, int beg, int end, uint64_t &curr_off, int &i, int &n_seeks, pair64_t *off, int n_off){

	if (off == 0) return 0;
		
	if (curr_off == 0 || curr_off >= off[i].v) { // then jump to the next chunk
		if (i == n_off - 1) return 0; // no more chunks
		if (i >= 0) assert(curr_off == off[i].v); // otherwise bug
		if (i < 0 || off[i].v != off[i+1].u) { // not adjacent chunks; then seek
			bam_seek(fp, off[i+1].u, SEEK_SET);
			curr_off = bam_tell(fp);
			++n_seeks;
		}
		++i;
	}
	if ((ret = bam_read1(fp, b)) > 0) {
		curr_off = bam_tell(fp);
		if (b->core.tid != tid || b->core.pos >= end) return 0; // no need to proceed
		else if (is_overlap(beg, end, b)) return 1;
	} 
	else 
		return 0;
	return 1;
}

// read bam files all together, and merge them
int MergeBams_prep(string *fn, int n, bamFile *fp, heap1_t *heap){
	for(i = 0; i!=n; ++i){
		heap1_t *h;
		fp[i] = bam_open(fn[i], "r");
		if(fp[i] == 0){
			int j;
			cout << "[bam_merge_core] fail to open file " << fn[i] << "\n";
			for(int j = 0; j < i; ++j)
				bam_close(fp[j]);
			free(fp);
			free(heap);
			return 0;
		}
		h = heap + i;
		h->i = i;
		h->b = (bam1_*)calloc(1,sizeof(bam1_t));
		if(bam_read1(fp[i],h->b) >= 0){
			h->pos = ((uint64_t)h->b->core.tid <<32) | (uint32_t)h->b->core.pos << 1 | bam1_strand(h->b);
			h->idx = idx++;
		}
		else h->pos = HEAP_EMPTY;
	}
	
	ks_heapmake(heap, n, heap);
	return 1;
}
	
// read bam files all together by one particular chromosome, and merge them
int MergeBamsChr_prep(string *fn, int n, bamFile *fp, heap1_t *heap, string chr_str, int *tid, int *beg, int *end, samfile_t **in, pair64_t **off, int *n_off){
	for(i = 0; i!=n; ++i){
		heap1_t *h;
		fp[i] = ReadBamChr_prep(chr_str, fn[i], tid[i], beg[i], end[i], in[i], off[i], n_off[i]);
		if(fp[i] == 0){
			int j;
			cout << "[bam_merge_core] fail to open file " << fn[i] << "\n";
			for(int j = 0; j < i; ++j)
				bam_close(fp[j]);
			free(fp);
			free(heap);
			return 0;
		}
		h = heap + i;
		h->i = i;
		h->b = (bam1_*)calloc(1,sizeof(bam1_t));
		if(bam_read1(fp[i],h->b) >= 0){
			h->pos = ((uint64_t)h->b->core.tid <<32) | (uint32_t)h->b->core.pos << 1 | bam1_strand(h->b);
			h->idx = idx++;
		}
		else h->pos = HEAP_EMPTY;
	}
	
	ks_heapmake(heap, n, heap);
	return 1;
}
	
	
	
// this function is good
string get_from_line(string line,string search,int flag){
	size_t pos = line.find(search);
	if(pos != string::npos){
		if(flag == 1){
			pos = pos + search.length();
			if(line.substr(pos,1).compare(":") == 0){
				size_t pos_end = line.find("\t",pos);
				if(pos_end == string::npos)
					pos_end = line.find("\0",pos);
				size_t n = pos_end - pos - 1;
				return line.substr(pos+1,n);
			}
			else{
				return "NA";
			}
		}
		else{
			pos = pos + search.length();
			size_t pos_begin = line.find(":", pos);
			if(pos_begin != string::npos){
				size_t pos_end = line.find("\t",pos_begin);
				if(pos_end == string::npos)
					pos_end = line.find("\0", pos_begin);
				size_t n = pos_end - pos_begin - 1;
				return line.substr(pos+1,n);
			}
			else
				return "NA";
		}
	}
	else
		return "NA";
	return "NA";
}					
	
// this function is good
string get_from_line_two(string line,string search1,string search2,int flag2){
	size_t pos = line.find(search1);
	size_t pos2;
	if(pos != string::npos){
		pos = pos + search1.length();
		pos2 = line.find(search2,pos);
		if(pos2 != string::npos){
			if(flag2 == 1){
				pos2 = pos2 + search2.length();
				if(line.substr(pos2,1).compare(":") == 0){
					size_t pos_end = line.find("\t",pos2);
					if(pos_end == string::npos)
						pos_end = line.find("\0",pos2);
					size_t n = pos_end - pos2 - 1;
					return line.substr(pos2+1, n);
				}				
				else
					return "NA";
				}
			else{
				pos2 = pos2 + search2.length();
				size_t pos_begin = line.find(":", pos2);
				if(pos2 != string::npos){
					size_t pos_end = line.find("\t",pos2);
					if(pos_end == string::npos)
						pos_end = line.find("\0",pos2);
					size_t n = pos_end - pos_begin - 1;
					return line.substr(pos_begin, n);
				}
				else
					return "NA";
			}
		}
		else
			return "NA";
	}
	else
		return "NA";
	return "NA";
}
					

							
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
t.readlen		:		b->core.l_qseq
*/
