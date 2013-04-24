#include "BreakDancerMax.h"

#include "breakdancer/BamReader.hpp"
#include "breakdancer/Region.hpp"
#include "breakdancer/Read.hpp"
#include "breakdancer/LegacyBamReader.hpp"
#include "breakdancer/BDConfig.hpp"

#include "version.h"
#include <stdio.h>
#include <stdlib.h>
#include <boost/shared_ptr.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <memory>
#include <stdexcept>
#include <assert.h>

#ifndef SCORE_FLOAT_TYPE
# define SCORE_FLOAT_TYPE double
#endif

KSORT_INIT(heap, heap1_t, heap_lt)

/*
# Data structure menagerie
## Preliminary counting
+ nread -> number of reads for each library
+ x_readcounts -> number of reads of each orientation, per library
## Analysis
+ reg_idx -> region id generating variable
+ reg_name -> vector; contains coordinates of the region and the number of normal reads within it.
+ reg_seq -> array of read information underlying a region
+ read -> hash of readnames with array of region ids. Stores which two regions are associated with a readpair
+ regs -> hash storing each region's associated reads
## buildConnection
+ link -> hash of node ids linking two nodes and their weights
+ clink -> a copy of link
+ free_nodes -> list of nodes to remove?
+ nodepair -> additional copy of links between two nodes
+ read_pair -> hash on read name containing read information
+ nonsupportives -> array of reads. The first time a read is found, it is pushed onto here. They are deleted as mates are found in the pair. Thus, after processing a graph, all that's left is reads not supporting that node pair.
+ snodes -> array of node ids from nodepair can contain the same node twice


*/

using boost::math::cdf;
using boost::math::complement;
using boost::math::poisson_distribution;
using boost::math::chi_squared;
using boost::shared_ptr;

using namespace std;

typedef SCORE_FLOAT_TYPE real_type;
real_type _max_kahan_err = 0.0;

real_type ComputeProbScore(
    vector<int> &rnode,
    map<string,int> &rlibrary_readcount,
    uint32_t type,
    map<uint32_t,
    map<string,int> > &x_readcounts,
    uint32_t reference_len,
    int fisher,
    map<int,
    vector<int> > &reg_name
);

namespace {
    auto_ptr<IBamReader> openBam(std::string const& path, Options const& opts) {
        if (opts.chr == "0")
            return auto_ptr<IBamReader>(new BamReader(path));
        else
            return auto_ptr<IBamReader>(new RegionLimitedBamReader(path, opts.chr.c_str()));
    }

    Options parseArguments(int argc, char** argv) {
        int c;
        Options opts;
        while((c = getopt(argc, argv, "o:s:c:m:q:r:x:b:ep:tfd:g:lCahy:")) >= 0) {
            switch(c) {
                case 'o': opts.chr = optarg; break;
                case 's': opts.min_len = atoi(optarg); break;
                case 'c': opts.cut_sd = atoi(optarg); break;
                case 'm': opts.max_sd = atoi(optarg); break;
                case 'q': opts.min_map_qual = atoi(optarg); break;
                case 'r': opts.min_read_pair = atoi(optarg); break;
                case 'x': opts.seq_coverage_lim = atoi(optarg); break;
                case 'b': opts.buffer_size = atoi(optarg); break;
                case 'e': opts.learn_par = true; break;
                case 'p': opts.prior_prob = atof(optarg); break;
                case 't': opts.transchr_rearrange = true; break;
                case 'f': opts.fisher = true; break;
                case 'd': opts.prefix_fastq = optarg; break;
                case 'g': opts.dump_BED = optarg; break;
                case 'l': opts.Illumina_long_insert = true; break;
                //case 'C': opts.Illumina_to_SOLiD = true; break;
                case 'a': opts.CN_lib = true; break;
                case 'h': opts.print_AF = true; break;
                case 'y': opts.score_threshold = atoi(optarg); break;
                default: fprintf(stderr, "Unrecognized option '-%c'.\n", c);
                    exit(1);
            }
        }

        // FIXME: instead of printing out defaults, this will print any partial options
        // specified. Let's try to get clearance to use boost::program_options or something
        // more reasonable.
        if(optind == argc) {
            fprintf(stderr, "\nbreakdancer-max version %s (commit %s)\n\n", __g_prog_version, __g_commit_hash);
            fprintf(stderr, "Usage: breakdancer-max <analysis.config>\n\n");
            fprintf(stderr, "Options: \n");
            fprintf(stderr, "       -o STRING       operate on a single chromosome [all chromosome]\n");
            fprintf(stderr, "       -s INT          minimum length of a region [%d]\n", opts.min_len);
            fprintf(stderr, "       -c INT          cutoff in unit of standard deviation [%d]\n", opts.cut_sd);
            fprintf(stderr, "       -m INT          maximum SV size [%d]\n", opts.max_sd);
            fprintf(stderr, "       -q INT          minimum alternative mapping quality [%d]\n", opts.min_map_qual);
            fprintf(stderr, "       -r INT          minimum number of read pairs required to establish a connection [%d]\n", opts.min_read_pair);
            fprintf(stderr, "       -x INT          maximum threshold of haploid sequence coverage for regions to be ignored [%d]\n", opts.seq_coverage_lim);
            fprintf(stderr, "       -b INT          buffer size for building connection [%d]\n", opts.buffer_size);
            //fprintf(stderr, "    -e INT    learn parameters from data before applying to SV detection [%d]\n", opts.learn_par);
            //fprintf(stderr, "    -p FLOAT    prior probability of SV [%f]\n", prior_prob);
            fprintf(stderr, "       -t              only detect transchromosomal rearrangement, by default off\n");
            //fprintf(stderr, "    -f INT    use Fisher's method to combine P values from multiple library [%d]\n", opts.fisher);
            fprintf(stderr, "       -d STRING       prefix of fastq files that SV supporting reads will be saved by library\n");
            fprintf(stderr, "       -g STRING       dump SVs and supporting reads in BED format for GBrowse\n");
            fprintf(stderr, "       -l              analyze Illumina long insert (mate-pair) library\n");
            fprintf(stderr, "       -a              print out copy number and support reads per library rather than per bam, by default off\n");
            fprintf(stderr, "       -h              print out Allele Frequency column, by default off\n");
            fprintf(stderr, "       -y INT          output score filter [%d]\n", opts.score_threshold);
            //fprintf(stderr, "    -C INT    change system default from Illumina to SOLiD [%d]\n", opts.Illumina_to_SOLiD);
            //fprintf(stderr, "Version: %s\n", version);
            fprintf(stderr, "\n");
            exit(1);
        }

        opts.platform = opts.Illumina_to_SOLiD ? "solid" : "illumina";

        return opts;
    }
}

// main function
int main(int argc, char *argv[]) {
    Options opts = parseArguments(argc, argv);

    // define the map SVtype
    map<string, string> SVtype;
    if(opts.Illumina_long_insert) {
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
    map<string, string> readgroup_platform;
    map<string,string> ReadsOut;
    int d = 1e8;// global

    // configure file
    ifstream CONFIG;
    CONFIG.open(argv[optind]);
    BDConfig __cfg(CONFIG);
    CONFIG.close();

    CONFIG.open(argv[optind]);
    //CONFIG.open("/gscuser/xfan/kdevelop/BreakDancerMax/debug/src/configure2.cfg");
    char line_[512]; // each line of the config file
    if(CONFIG.is_open())
    {
        while(CONFIG.good())
        {
            CONFIG.getline(line_, 512);
            //string line = char2str(line_);
            string line(line_);
            if(line.length()==0)
                break;
            // analyze the line
            string fmap = get_from_line(line,"map",1);
            string mean_ = get_from_line(line,"mean",0);
            string std_ = get_from_line(line,"std",0);
            string readlen_ = get_from_line(line,"readlen",0);
            string upper_ = get_from_line(line,"upp",0);
            string lower_ = get_from_line(line,"low",0);
            string mqual_ = get_from_line_two(line,"map","qual",0);
            string lib = get_from_line(line,"lib",0);
            float mean=0.0f,std=0.0f,readlen=0.0f,upper=0.0f,lower=0.0f;
            int mqual;
            if(lib.compare("NA")==0)
                lib = get_from_line(line,"samp",0);


            string readgroup = get_from_line(line,"group",1);
            if(readgroup.compare("NA")==0)
                readgroup = lib;
            readgroup_library[readgroup] = lib;

            string platform = get_from_line(line,"platform",1);
            if(opts.Illumina_to_SOLiD)
                readgroup_platform[readgroup] = "solid";
            else
                readgroup_platform[readgroup] = "illumina";
            readgroup_platform[readgroup] = platform;

            string exe = get_from_line(line,"exe",0);
            if(!opts.prefix_fastq.empty()) {
                //ofstream ReadsOut[lib.append("1")](prefix_fastq.append(lib).append(".1.fastq"), ios::out | ios::app | ios::binary);
                ReadsOut[lib + "1"] = opts.prefix_fastq + "." + lib + ".1.fastq";
                ofstream ReadsOutTmp;
                ReadsOutTmp.open(ReadsOut[lib+"1"].c_str());
                if(!ReadsOutTmp.is_open())
                    cout << "unable to open " << opts.prefix_fastq << "." << lib << ".1.fastq, check write permission\n";
                ReadsOutTmp.close();
                ReadsOut[lib + "2"] = opts.prefix_fastq + "." + lib + ".2.fastq";
                ReadsOutTmp.open(ReadsOut[lib+"2"].c_str());
                if(!ReadsOutTmp.is_open())
                    cout << "unable to open " << opts.prefix_fastq << "." << lib << ".2.fastq, check write permission\n";
                ReadsOutTmp.close();
            }

            libmaps[lib] = fmap;
            if(mqual_.compare("NA")){
                mqual = atoi(mqual_.c_str());
                mapQual[lib] = mqual;
            }
            fmaps[fmap] = lib;

            if(mean_.compare("NA") && std_.compare("NA")){
                mean = atof(mean_.c_str());
                std = atof(std_.c_str());
                if(!upper_.compare("NA") || !lower_.compare("NA")){
                    upper = mean + std*float(opts.cut_sd);
                    lower = mean - std*float(opts.cut_sd);
                    lower = lower > 0 ? lower : 0;//atof(lower_.c_str()) > 0 ? atof(lower_.c_str()):0; this is not related with lower_
                }
                else{
                    upper = atof(upper_.c_str());
                    lower = atof(lower_.c_str());
                }
            }

            if(readlen_.compare("NA")){
                readlen = atof(readlen_.c_str());
            }
            max_readlen = max_readlen < readlen ? readlen:max_readlen;

            mean_insertsize[lib] = mean;
            std_insertsize[lib] = std;
            uppercutoff[lib] = upper;
            lowercutoff[lib] = lower;
            readlens[lib] = readlen;
            if(exes.find(fmap) == exes.end())
                exes[fmap] = exe.compare("NA")?exe:"cat";
            else if(exes[fmap].compare(exe) != 0){
                cout << "Please use identical exe commands to open the same map file.\n";
                return 1;
            }

            int tmp = mean - readlen*2;    // this determines the mean of the max of the SV flanking region
            d = d<tmp ? d:tmp;
        }
    }

    CONFIG.close();

    if(d < 50)
        d = 50;

    vector<string> format;
    vector<string>::const_iterator it_format;

    map<string,int> cmds;

    // go through the iteration of fmaps
    map<string,string>::const_iterator ii;
    for(ii=fmaps.begin(); ii!=fmaps.end(); ++ii)
    {
        string exe = exes[(*ii).first];
        cmds[exe] = 0;
        if(exe.find("maq")!=string::npos)
            format.insert(format.end(),1,"maq");
        else
            format.insert(format.end(),1,"sam");
    }

    if(opts.learn_par) {
        EstimatePriorParameters(opts, fmaps, readgroup_library, mean_insertsize, std_insertsize, uppercutoff, lowercutoff, readlens, readgroup_platform);
    }

    uint32_t reference_len = 1;
    map<string, int> nreads;
    map<string, int> nreads_;    // only to compute density per bam

    int i = 0;
    bam1_t *b = bam_init1();
    string format_ = "sam";
    string alt = "";
    vector<string> maps;

     for(ii=fmaps.begin(); ii!=fmaps.end(); ++ii)
     {

        maps.push_back((*ii).first);
         uint32_t ref_len = 0;
         string exe = exes[maps[i]];
         i ++;
         cmds[exe] ++;

         int p_pos = 0;
         const char *p_chr = 0;
        string bam_name = (*ii).first;



        // no chromosome defined
        // convert the entire bam file
        // XXX: Note: the original code would 'continue' if a particular
        // bam could not be opened. It would also spew the bam name into
        // the output, which probably messes it up anyway. We are going
        // to crash. If it is determined that this is wrong, we should
        // declare a variable for the BamReader, heap allocate it (with
        // boost::scoped_ptr/std::unique_ptr/std::auto_ptr) in a try
        // block and continue on exceptions.
        auto_ptr<IBamReader> reader(openBam(bam_name, opts));

        while (reader->next(b) > 0) {
            int same_tid = (b->core.tid == b->core.mtid)? 1:0;
            breakdancer::Read aln2(b, format_, readgroup_platform, readgroup_library);
            vector<string> aln_return = AlnParser(b, format_, alt, readgroup_platform, same_tid, opts.platform);
            assert(b->core.flag == aln2.bdflag);

            string ori = aln_return[1];
            string readgroup = aln_return[0];

            //FIXME this could be filtered out on reading the BAM
            if(b->core.tid < 0){
                //cout << b->core.tid << "\tError: no correspondence of the chromosome to the header file!" << endl;
                continue;
            }

            //This seems weird below, but it could allow for merging of differently sorted bams
            //If we don't care if we support this then we would need to check that they are the same elsewhere
            //and then this could switch to a tid comparison...
            char const* new_seq_name = reader->header()->target_name[b->core.tid];
            if (p_chr) {
                if (strcmp(p_chr, new_seq_name) == 0) {
                    ref_len += b->core.pos - p_pos;
                }
            }
            p_pos = b->core.pos;
            p_chr = new_seq_name;

            //FIXME It would be better if this was done as part of the Read class
            //doing it in a way in which the read class doesn't have to know about how we are 
            //reading the files would be key though
            string lib;
            if(!(readgroup.empty()))
                lib = readgroup_library[readgroup];
            else{
                lib = fmaps[(*ii).first];
            }
            if(lib.length() == 0)
                continue;

            //Below seems like a bug to me. 18 means FR orientation plus paired end req met in MAQ
            //this would also allow 20 (RF and paired end req met) and 24 (RR and paired end req met)
            //the latter two are abnormal mappings and this would allow them but would skip 17 (FF and paired end req met)
            //I don't know if these were ever real world values but this oddness seems inconsistent with the intent
            // (^dlarson)
            //
            // This is the thing I was talking about that creeped me out
            // (relational cmps on bitvectors)
            // -tabbott
            //
            // Indeed. Having looked over the flags, my comment above does not apply...
            // -dlarson
            if(aln2.bdqual > opts.min_map_qual && (aln2.bdflag == breakdancer::NORMAL_FR || aln2.bdflag == breakdancer::NORMAL_RF)) {
                ++nreads[lib];
                if(opts.CN_lib == 0){
                    ++nreads_[libmaps[lib]];
                }
            }

            //XXX This seems weird to me as well. Why are we checking the 
            //mapping quality in two different places and performing some 
            //calculations before this is applied?
            //-dlarson
            if(mapQual.find(lib) != mapQual.end()){
                if(aln2.bdqual <= mapQual[lib])
                    continue;
            }
            else{
                if(aln2.bdqual <= opts.min_map_qual)
                    continue;
            }
            if(aln2.bdflag == breakdancer::NA)
                continue;
            if((opts.transchr_rearrange && aln2.bdflag != breakdancer::ARP_CTX) || aln2.bdflag == breakdancer::MATE_UNMAPPED || aln2.bdflag == breakdancer::UNMAPPED)
                continue;

            //It would be nice if this was pulled into the Read class as well
            //for now, let's just set the bdflag directly here since it is public
            if(opts.Illumina_long_insert){
                if(abs(b->core.isize) > uppercutoff[lib] && aln2.bdflag == breakdancer::NORMAL_RF) {
                    b->core.flag = 4;   //why isn't this ARP_FR_big_insert ie. 2?
                    aln2.bdflag = breakdancer::ARP_RF;
                }
                if(abs(b->core.isize) < uppercutoff[lib] && aln2.bdflag == breakdancer::ARP_RF) {
                    b->core.flag = 20;
                    aln2.bdflag = breakdancer::NORMAL_RF;
                }
                if(abs(b->core.isize) < lowercutoff[lib] && aln2.bdflag == breakdancer::NORMAL_RF) {
                    b->core.flag = 3;
                    aln2.bdflag = breakdancer::ARP_FR_small_insert; //FIXME this name doesn't make a whole lot of sense here
                }
            }
            else{
                if(abs(b->core.isize) > uppercutoff[lib] && aln2.bdflag == breakdancer::NORMAL_FR) {
                    b->core.flag = 2;
                    aln2.bdflag = breakdancer::ARP_FR_big_insert;
                }
                if(abs(b->core.isize) < uppercutoff[lib] && aln2.bdflag == breakdancer::ARP_FR_big_insert) {
                    b->core.flag = 18;
                    aln2.bdflag = breakdancer::NORMAL_FR;
                }
                if(abs(b->core.isize) < lowercutoff[lib] && aln2.bdflag == breakdancer::NORMAL_FR) {
                    b->core.flag = 3;
                    aln2.bdflag = breakdancer::ARP_FR_small_insert;
                }
            }

            if(b->core.flag == 18 || b->core.flag == 20 || b->core.flag == 130){
                continue;
            }

            if(x_readcounts.find(b->core.flag) != x_readcounts.end() && x_readcounts[b->core.flag].find(lib) != x_readcounts[b->core.flag].end())
                x_readcounts[b->core.flag][lib] ++;
            else
                x_readcounts[b->core.flag][lib] = 1;
        }
        reader.reset(); // free bam reader

        if(ref_len == 0)
            cout << (*ii).second << " does not contain legitimate paired end alignment. Please check that you have the correct paths and the map/bam files are properly formated and indexed.";

        if(reference_len < ref_len)
            reference_len = ref_len;
    }
    bam_destroy1(b);

    // need to read the total base

    float total_phy_cov = 0;
    float total_seq_cov = 0;

    cout << "#Software: " << __g_prog_version << endl;
    cout << "#Command: ";
    for(int i=0;i<argc;i++) {
        cout << argv[i] << " ";
    }
    cout << endl;
    cout << "#Library Statistics:" << endl;
    map<string,int>::const_iterator nreads_ii;
    map<string,float> read_density;
    for(nreads_ii=nreads.begin(); nreads_ii!=nreads.end(); ++nreads_ii)
    {
        string lib = (*nreads_ii).first;
        float sequence_coverage = float(nreads[lib]*readlens[lib])/float(reference_len);
        total_seq_cov += sequence_coverage;

        // compute read_density
        if(opts.CN_lib == 1){
            if(nreads.find(lib) != nreads.end())
                read_density[lib] = float(nreads[lib])/float(reference_len);
            else{
                read_density[lib] = 0.000001;
                cout << lib << " does not contain any normals" << endl;
            }
        }
        else{
            if(nreads_.find(libmaps[lib]) != nreads_.end())
                read_density[libmaps[lib]] = float(nreads_[libmaps[lib]])/float(reference_len);
            else{
                read_density[libmaps[lib]] = 0.000001;
                cout << lib << " does not contain any normals" << endl;
            }
        }

        float physical_coverage = float(nreads[lib]*mean_insertsize[lib])/float(reference_len)/2;
        total_phy_cov += physical_coverage;

        int nread_lengthDiscrepant = -1;

        if(x_readcounts.find(2) != x_readcounts.end() && x_readcounts[2].find(lib) != x_readcounts[2].end())
            nread_lengthDiscrepant = x_readcounts[2][lib];
        if(x_readcounts.find(3) != x_readcounts.end() && x_readcounts[3].find(lib) != x_readcounts[3].end()){
            if(nread_lengthDiscrepant == -1)
                nread_lengthDiscrepant = 0;
            nread_lengthDiscrepant += x_readcounts[3][lib];
        }
        float tmp = (nread_lengthDiscrepant > 0)?(float)reference_len/(float)nread_lengthDiscrepant:50;
        d = d<tmp?d:tmp;

        //printf("#%s\tmean:%.3f\tstd:%.3f\tuppercutoff:%.3f\tlowercutoff:%.3f\treadlen:%.3f\tlibrary:%s\treflen:%d\tseqcov:%.3fx\tphycov:%.3fx", libmaps[lib],mean_insertsize[lib],std_insertsize[lib],uppercutoff[lib],lowercutoff[lib],readlens[lib],lib,reference_len, sequence_coverage,physical_coverage);
        cout << "#" << libmaps[lib] << "\tmean:" << mean_insertsize[lib] << "\tstd:" << std_insertsize[lib] << "\tuppercutoff:" << uppercutoff[lib] << "\tlowercutoff:" << lowercutoff[lib] << "\treadlen:" << readlens[lib] << "\tlibrary:" << lib << "\treflen:" << reference_len << "\tseqcov:" << sequence_coverage << "\tphycov:" << physical_coverage;

        map<uint32_t,map<string,int> >::const_iterator x_readcounts_ii;
        for(x_readcounts_ii = x_readcounts.begin(); x_readcounts_ii!=x_readcounts.end(); ++x_readcounts_ii){
            uint32_t t = (*x_readcounts_ii).first;// get the first key out, which is a member of recflags

            map<string,int>::const_iterator found = x_readcounts_ii->second.find(lib);
            if (found != x_readcounts_ii->second.end())
                printf("\t%d:%d",t,found->second);
        }
        printf("\n");
    }

    printf("#Chr1\tPos1\tOrientation1\tChr2\tPos2\tOrientation2\tType\tSize\tScore\tnum_Reads\tnum_Reads_lib");
    if(opts.print_AF == 1)
        printf("\tAllele_frequency");
    if(opts.CN_lib == 0){
        for(vector<string>::const_iterator it_map = maps.begin(); it_map != maps.end(); it_map++){
            size_t tmp = (*it_map).rfind("/");
            if(tmp!=string::npos)
                cout << "\t" << (*it_map).substr(tmp + 1);
            else
                cout << "\t" << *it_map;
        }
    }

    cout << "\n";


    BreakDancerData bdancer;

    map<string, uint32_t > possible_fake_data;
    map<int, vector<vector<string> > > regs;//global in analysis
    map<string, vector<int> > read;// global in analysis
    map<int,vector<int> > reg_name;// global in analysis
    vector<vector<string> > reg_seq; // global need to see if it's the key or value of one of the above global. should be a string

    int idx_buff = 0;// global
    int reg_idx = 0;// global  ################# node index here ###################
    int normal_switch = 0; // global
    int nnormal_reads = 0; // global
    uint32_t ntotal_nucleotides = 0; // global

    max_readlen = 0;

    // first, need to merge the bam files into one big string seperated by blank, and return the number
    int n = fmaps.size();

    if(n == 0){
        cout << "wrong: no bam file!\n";
        return 0;
    }
    //################# node index here ###################
    else if(n == 1){
        auto_ptr<IBamReader> reader(openBam(maps[0], opts));
        int r;
        bam1_t *b = bam_init1();
        while ((r = reader->next(b)) >= 0) {
            if(b->core.tid < 0)
                continue;
            int same_tid = b->core.tid == b->core.mtid ? 1:0;
            breakdancer::Read aln2(b, format_, readgroup_platform, readgroup_library);
            vector<string> aln_return = AlnParser(b, format_, alt, readgroup_platform, same_tid, opts.platform);
            assert(b->core.flag == aln2.bdflag);
            string readgroup = aln_return[0];
            string ori = aln_return[1];
            string library = (!readgroup.empty())?readgroup_library[readgroup]:((*(fmaps.begin())).second);

            if(!library.empty()){
                Analysis(
                    opts,
                    bdancer,
                    library,
                    b,
                    reg_seq,
                    reg_name,
                    read,
                    regs,
                    &idx_buff,
                    &nnormal_reads,
                    &normal_switch,
                    &reg_idx,
                    x_readcounts,
                    reference_len,
                    ReadsOut,
                    mean_insertsize,
                    SVtype,
                    mapQual,
                    uppercutoff,
                    lowercutoff,
                    d,
                    &max_readlen,
                    ori,
                    reader->header(),
                    &ntotal_nucleotides,
                    read_density,
                    possible_fake_data,
                    libmaps,
                    maps
                );
            }
        }
        if(reg_seq.size() != 0){
            do_break_func(
                opts,
                bdancer,
                reg_seq,
                reg_name,
                read,
                regs,
                &idx_buff,
                opts.buffer_size,
                &nnormal_reads,
                opts.min_len,
                &reg_idx,
                opts.transchr_rearrange,
                opts.min_map_qual,
                opts.Illumina_long_insert,
                opts.prefix_fastq,
                x_readcounts,
                reference_len,
                opts.fisher,
                ReadsOut,
                mean_insertsize,
                SVtype,
                mapQual,
                uppercutoff,
                lowercutoff,
                opts.max_sd,
                d,
                opts.min_read_pair,
                opts.dump_BED,
                &max_readlen,
                reader->header(),
                opts.seq_coverage_lim,
                &ntotal_nucleotides,
                read_density,
                opts.CN_lib,
                libmaps,
                maps,
                opts.print_AF,
                opts.score_threshold
            );
        }

        buildConnection(
            bdancer,
            read,
            reg_name,
            regs,
            x_readcounts,
            reference_len,
            opts.fisher,
            opts.min_read_pair,
            opts.dump_BED,
            max_readlen,
            opts.prefix_fastq,
            ReadsOut,
            SVtype,
            mean_insertsize,
            reader->header(),
            read_density,
            opts.CN_lib,
            maps,
            opts.print_AF,
            opts.score_threshold,
            libmaps
        );
        bam_destroy1(b);
    }
    else{
        string *big_bam = new string[n];
        int i_tmp = 0;
        for(map<string, string>::const_iterator ii_fmaps = fmaps.begin(); ii_fmaps != fmaps.end(); ii_fmaps++)
            big_bam[i_tmp++] = (*ii_fmaps).first;

        heap1_t *heap;
        heap = (heap1_t*)calloc(n,sizeof(heap1_t));

        uint64_t idx = 0;

        if(opts.chr == "0") {
            samfile_t **in = new samfile_t *[n];
             // open pipe, improvement made by Ben Oberkfell (boberkfe@genome.wustl.edu) samtools merge - in1.bam in2.bam in3.bam in_N.bam | samtools view - maq mapmerge
               // dig into merge samtools code and utilize what we needed

            if(MergeBams_prep(big_bam, n, in, heap, &idx)){
                bam1_t *b = heap->b;
                int skip_previous = 0;
                while(heap->pos != HEAP_EMPTY){

                    if(skip_previous == 0){
                        int same_tid = b->core.tid == b->core.mtid ? 1:0;
                        breakdancer::Read aln2(b, format_, readgroup_platform, readgroup_library);
                        vector<string> aln_return = AlnParser(b, format_, alt, readgroup_platform, same_tid, opts.platform);
                        assert(b->core.flag == aln2.bdflag);
                        string ori = aln_return[1];
                        string readgroup = aln_return[0];
                        string library = (!readgroup.empty())?readgroup_library[readgroup]:((*(fmaps.begin())).second);
                        //if(chr.empty() || chr.compare(b->core.tid)!=0) //this statement actually does nothing
                        //continue;
                        if(!library.empty()){
                            Analysis(
                                opts,
                                bdancer,
                                library,
                                b,
                                reg_seq,
                                reg_name,
                                read,
                                regs,
                                &idx_buff,
                                &nnormal_reads,
                                &normal_switch,
                                &reg_idx,
                                x_readcounts,
                                reference_len,
                                ReadsOut,
                                mean_insertsize,
                                SVtype,
                                mapQual,
                                uppercutoff,
                                lowercutoff,
                                d,
                                &max_readlen,
                                ori,
                                in[heap->i]->header,
                                &ntotal_nucleotides,
                                read_density,
                                possible_fake_data,
                                libmaps,
                                maps
                            );
                        }
                        /*else{
                         if(b->core.pos >= 43803315 && b->core.pos <= 103860201){
                         count_no_lib ++;
                         //cout << b->core.pos + 1 << endl;
                         }
                         }*/

                    }
                    int j = samread(in[heap->i],heap->b);
                    //    j = -1;
                    if(j >= 0){
                        heap -> pos = ((uint64_t)b->core.tid<<32) | (uint32_t)b->core.pos << 1 | bam1_strand(b);
                        heap -> idx = idx++;
                        b = heap->b;
                        if(b->core.tid < 0){
                            skip_previous = 1;
                            continue;
                        }
                        if(skip_previous == 1)
                            skip_previous = 0;
                        ks_heapadjust(heap, 0, n, heap);
                        //if(strcmp(in[heap->i]->header->target_name[b->core.tid],"NT_113888") == 0)
                        //int k = 0;
                    }
                    else if(j == -1){
                        heap->pos = HEAP_EMPTY;
                        //cout << "heap empty" << endl;
                        free(heap->b->data);
                        free(heap->b);
                        heap->b = 0;
                        //                                                cout << "here" << endl;
                    }
                    else
                        cout << "[bam_merge_core] " << big_bam[heap->i] << " is truncated. Continue anyway.\n";


                }
            }
            //cout << "build connection:" << endl;
            if(reg_seq.size() != 0){
                do_break_func(
                    opts,
                    bdancer,
                    reg_seq,
                    reg_name,
                    read,
                    regs,
                    &idx_buff,
                    opts.buffer_size,
                    &nnormal_reads,
                    opts.min_len,
                    &reg_idx,
                    opts.transchr_rearrange,
                    opts.min_map_qual,
                    opts.Illumina_long_insert,
                    opts.prefix_fastq,
                    x_readcounts,
                    reference_len,
                    opts.fisher,
                    ReadsOut,
                    mean_insertsize,
                    SVtype,
                    mapQual,
                    uppercutoff,
                    lowercutoff,
                    opts.max_sd,
                    d,
                    opts.min_read_pair,
                    opts.dump_BED,
                    &max_readlen,
                    in[0]->header,
                    opts.seq_coverage_lim,
                    &ntotal_nucleotides,
                    read_density,
                    opts.CN_lib,
                    libmaps,
                    maps,
                    opts.print_AF,
                    opts.score_threshold
                    );
            }
            buildConnection(
                bdancer,
                read,
                reg_name,
                regs,
                x_readcounts,
                reference_len,
                opts.fisher,
                opts.min_read_pair,
                opts.dump_BED,
                max_readlen,
                opts.prefix_fastq,
                ReadsOut,
                SVtype,
                mean_insertsize,
                in[0]->header,
                read_density,
                opts.CN_lib,
                maps,
                opts.print_AF,
                opts.score_threshold,
                libmaps
            );

            for (int k = 0; k < n; ++k) {
                samclose(in[k]);
            }
            free(in);
        }

        //############### find the designated chromosome and put all of them together for all bam files #############
        else{//#customized merge & sort

            // FIXME: make this not use shared_ptr
            // c++11 unique_ptr would be better, but we're c++98/03ing it for now
            typedef vector<shared_ptr<RegionLimitedBamReader> > ReaderVecType;
            ReaderVecType sp_readers;

            vector<IBamReader*> readers;

            int tid = -1;
            int reg_end = -1;
            for(int i = 0; i!=n; ++i){
                shared_ptr<RegionLimitedBamReader> reader(
                    new RegionLimitedBamReader(big_bam[i], opts.chr.c_str()));

                // XXX: these are low detail err msgs, but this code is going away soon
                if (tid >= 0 && reader->tid() != tid) {
                    throw runtime_error("tid mismatch");
                }

                if (reg_end >= 0 && reader->end() != reg_end) {
                    throw runtime_error("region end mismatch");
                }

                tid = reader->tid();
                reg_end = reader->end();

                sp_readers.push_back(reader);
                readers.push_back(reader.get());

                heap1_t *h;
                h = heap + i;
                h->i = i;
                h->b = bam_init1();
                if (sp_readers.back()->next(h->b)) {
                    h->pos = ((uint64_t)h->b->core.tid <<32) | (uint32_t)h->b->core.pos << 1 | bam1_strand(h->b);
                    h->idx = idx++;
                }
                else h->pos = HEAP_EMPTY;
            }
            ks_heapmake(heap, n, heap);
            int merge_tmp = 1;

            if(merge_tmp){
                bam1_t *b = heap->b;
                int skip_previous = 0;
                int n_ava_heap = n;
                while(heap->pos != HEAP_EMPTY){

                    if(skip_previous == 0){
                        if( b->core.tid == tid && b->core.pos < reg_end) {
                            int same_tid = b->core.tid == b->core.mtid ? 1:0;
                            breakdancer::Read aln2(b, format_, readgroup_platform, readgroup_library);
                            vector<string> aln_return = AlnParser(b, format_, alt, readgroup_platform, same_tid, opts.platform);
                            assert(b->core.flag == aln2.bdflag);
                            string ori = aln_return[1];
                            string readgroup = aln_return[0];
                            string library = (!readgroup.empty())?readgroup_library[readgroup]:((*(fmaps.begin())).second);
                            //if(chr.empty() || chr.compare(b->core.tid)!=0) //this statement actually does nothing
                            //continue;
                            if(!library.empty()){
                                Analysis(
                                    opts,
                                    bdancer,
                                    library,
                                    b,
                                    reg_seq,
                                    reg_name,
                                    read,
                                    regs,
                                    &idx_buff,
                                    &nnormal_reads,
                                    &normal_switch,
                                    &reg_idx,
                                    x_readcounts,
                                    reference_len,
                                    ReadsOut,
                                    mean_insertsize,
                                    SVtype,
                                    mapQual,
                                    uppercutoff,
                                    lowercutoff,
                                    d,
                                    &max_readlen,
                                    ori,
                                    readers[0]->header(),
                                    &ntotal_nucleotides,
                                    read_density,
                                    possible_fake_data,
                                    libmaps,
                                    maps
                                );
                            }
                        }
                    }

                    int j = readers[heap->i]->next(b);
                    if(j > 0){
                        heap -> pos = ((uint64_t)b->core.tid<<32) | (uint32_t)b->core.pos << 1 | bam1_strand(b);
                        heap -> idx = idx++;
                        b = heap->b;
                        if(b->core.tid < 0){
                            skip_previous = 1;
                            continue;
                        }
                        if(skip_previous == 1)
                            skip_previous = 0;
                        ks_heapadjust(heap, 0, n_ava_heap, heap);
                    }
                    else {
                        n_ava_heap --;
                        //int n_ava_heap = n-1;
                        int jj = 0;
                        while(n_ava_heap > 0 && jj == 0){
                            heap1_t tmp = heap[0];
                            heap[0] = heap[n_ava_heap];
                            heap[n_ava_heap] = tmp;
                            ks_heapadjust(heap, 0, n_ava_heap, heap);
                            jj = readers[heap->i]->next(b);
                            if(jj == 0){
                                n_ava_heap --;
                            }
                            else if(jj > 0){

                                heap -> pos = ((uint64_t)b->core.tid<<32) | (uint32_t)b->core.pos << 1 | bam1_strand(b);
                                heap -> idx = idx++;
                                b = heap->b;
                                if(b->core.tid < 0){
                                    skip_previous = 1;
                                    // continue;
                                }
                                if(skip_previous == 1)
                                    skip_previous = 0;
                                ks_heapadjust(heap, 0, n_ava_heap, heap);
                            }

                        }
                        if(n_ava_heap == 0){
                            heap->pos = HEAP_EMPTY;
                            free(heap->b->data);
                            free(heap->b);
                            heap->b = 0;
                        }
                    }
                }
            }
            if(reg_seq.size() != 0){
                do_break_func(
                    opts,
                    bdancer,
                    reg_seq,
                    reg_name,
                    read,
                    regs,
                    &idx_buff,
                    opts.buffer_size,
                    &nnormal_reads,
                    opts.min_len,
                    &reg_idx,
                    opts.transchr_rearrange,
                    opts.min_map_qual,
                    opts.Illumina_long_insert,
                    opts.prefix_fastq,
                    x_readcounts,
                    reference_len,
                    opts.fisher,
                    ReadsOut,
                    mean_insertsize,
                    SVtype,
                    mapQual,
                    uppercutoff,
                    lowercutoff,
                    opts.max_sd,
                    d,
                    opts.min_read_pair,
                    opts.dump_BED,
                    &max_readlen,
                    readers[0]->header(),
                    opts.seq_coverage_lim,
                    &ntotal_nucleotides,
                    read_density,
                    opts.CN_lib,
                    libmaps,
                    maps,
                    opts.print_AF,
                    opts.score_threshold
                    );
            }
            buildConnection(
                bdancer,
                read,
                reg_name,
                regs,
                x_readcounts,
                reference_len,
                opts.fisher,
                opts.min_read_pair,
                opts.dump_BED,
                max_readlen,
                opts.prefix_fastq,
                ReadsOut,
                SVtype,
                mean_insertsize,
                readers[0]->header(),
                read_density,
                opts.CN_lib,
                maps,
                opts.print_AF,
                opts.score_threshold,
                libmaps
            );

        }
        //cout << "release memory\n";

        free(heap);

        delete []big_bam;
        big_bam = NULL;
    }

    cerr << "Max Kahan error: " << _max_kahan_err << "\n";

    return 0;
}

// to take the rest of the reads and trying to pair up
void do_break_func(
    Options const& opts,
    BreakDancerData& bdancer,
    vector<vector<string> > const& reg_seq,
    map<int, vector<int> >& reg_name,
    map<string, vector<int> >& read,
    map<int, vector<vector<string> > > &regs,
    int *idx_buff,
    int buffer_size,
    int *nnormal_reads,
    int min_len,
    int *reg_idx,
    int transchr_rearrange,
    int min_map_qual,
    int Illumina_long_insert,
    string prefix_fastq,
    map<uint32_t, map<string,int> > &x_readcounts,
    uint32_t reference_len,
    int fisher,
    map<string, string> &ReadsOut,
    map<string, float> &mean_insertsize,
    map<string, string> &SVtype,
    map<string, int> &mapQual,
    map<string, float> &uppercutoff,
    map<string, float> &lowercutoff,
    int max_sd,
    int d,
    int min_read_pair,
    string dump_BED,
    int *max_readlen,
    bam_header_t* bam_header,
    int seq_coverage_lim,
    uint32_t *ntotal_nucleotides,
    map<string, float>& read_density,
    int CN_lib,
    map<string, string> libmaps,
    vector<string> maps,
    int print_AF,
    int score_threshold
    )
{

    int& begins = bdancer.begins;
    int& beginc = bdancer.beginc;
    int& lastc = bdancer.lastc;

    float seq_coverage = *ntotal_nucleotides/float(lastc - beginc + 1 + *max_readlen);
    if(lastc - beginc > min_len && seq_coverage < seq_coverage_lim){ // skip short/unreliable flnaking supporting regions
        // register reliable region and supporting reads across gaps
        int k = (*reg_idx) ++;
        reg_name[k].push_back(begins);
        reg_name[k].push_back(beginc);
        reg_name[k].push_back(lastc);
        reg_name[k].push_back(*nnormal_reads);

        vector<vector<string> > p;
        for(vector<vector<string> >::const_iterator it_reg_seq = reg_seq.begin(); it_reg_seq != reg_seq.end(); it_reg_seq ++){
            p.push_back(*it_reg_seq);
            string s = (*it_reg_seq)[0];
            read[s].push_back(k);
        }

        //this should replace both regs and reg_name
        region::Region new_region(begins, beginc, lastc, *nnormal_reads, reg_seq);

        regs[k] = p;
        (*idx_buff)++;
        if(*idx_buff > buffer_size){
            //cout << "build connection:" << endl;
            buildConnection(
                bdancer,
                read,
                reg_name,
                regs,
                x_readcounts,
                reference_len,
                fisher,
                min_read_pair,
                dump_BED,
                *max_readlen,
                prefix_fastq,
                ReadsOut,
                SVtype,
                mean_insertsize,
                bam_header,
                read_density,
                CN_lib,
                maps,
                print_AF,
                score_threshold,
                libmaps
            );
            *idx_buff = 0;
        }
    }
    else{
        if(reg_seq.size()>0){
            for(vector<vector<string> >::const_iterator it_reg_seq = reg_seq.begin(); it_reg_seq != reg_seq.end(); it_reg_seq ++){
                ///string s = get_item_from_string(*it_reg_seq,0);
                string s= (*it_reg_seq)[0];
                if(read.find(s) != read.end())
                    read.erase(read.find(s));
                //cout << "read erase: " << s << endl;
            }
        }
    }
}

// for each read, check if it is time to break and pair up the reads
void Analysis (
    Options const& opts,
    BreakDancerData& bdancer,
    string lib,
    bam1_t *b,
    vector<vector<string> > &reg_seq,
    map<int, vector<int> > &reg_name,
    map<string, vector<int> > &read,
    map<int, vector<vector<string> > > &regs,
    int *idx_buff,
    int *nnormal_reads,
    int *normal_switch,
    int *reg_idx,
    map<uint32_t, map<string, int> > &x_readcounts,
    uint32_t reference_len,
    map<string, string> &ReadsOut,
    map<string, float> &mean_insertsize,
    map<string, string> &SVtype,
    map<string, int> &mapQual,
    map<string, float> &uppercutoff,
    map<string, float> &lowercutoff,
    int d,
    int *max_readlen,
    string ori,
    bam_header_t* bam_header,
    uint32_t *ntotal_nucleotides,
    map<string, float> &read_density,
    map<string, uint32_t> &possible_fake_data,
    map<string, string> libmaps,
    vector<string> maps
    )
{
    // for now, we can just set up references to the data struct so we
    // don't have to modify too much code
    int buffer_size = opts.buffer_size;
    int min_len = opts.min_len;
    bool transchr_rearrange = opts.transchr_rearrange;
    int min_map_qual = opts.min_map_qual;
    int Illumina_long_insert = opts.Illumina_long_insert;
    string const& prefix_fastq = opts.prefix_fastq;
    bool fisher = opts.fisher;
    int max_sd = opts.max_sd;
    int min_read_pair = opts.min_read_pair;
    string const& dump_BED = opts.dump_BED;
    int seq_coverage_lim = opts.seq_coverage_lim;
    bool CN_lib = opts.CN_lib;
    bool print_AF = opts.print_AF;
    int score_threshold = opts.score_threshold;

    int& begins = bdancer.begins;
    int& beginc = bdancer.beginc;
    int& lasts = bdancer.lasts;
    int& lastc = bdancer.lastc;
    map<string, uint32_t>& nread_ROI = bdancer.nread_ROI;
    map<int, map<string, uint32_t> >& read_count_ROI_map = bdancer.read_count_ROI_map;
    map<string, uint32_t>& nread_FR = bdancer.nread_FR;
    map<int, map<string, uint32_t> >& read_count_FR_map = bdancer.read_count_FR_map;


    string bam_name = libmaps[lib];
    //main analysis code

    // region between last and next begin
    // Store readdepth in nread_ROI by bam name (no per library calc) or by library
    // I believe this only counts normally mapped reads
    if(b->core.qual > min_map_qual && b->core.flag < 32 && b->core.flag >= 18){
        if(CN_lib == 1){
            if(nread_ROI.find(lib) == nread_ROI.end())
                nread_ROI[lib] = 1;
            else
                nread_ROI[lib] ++;
        }
        else{
            if(nread_ROI.find(bam_name) == nread_ROI.end())
                nread_ROI[bam_name] = 1;
            else
                nread_ROI[bam_name] ++;
        }
    // This stores the same counts. Not sure why.
    if(CN_lib == 1){
            if(possible_fake_data.find(lib) == possible_fake_data.end())
                possible_fake_data[lib] = 1;
            else
                possible_fake_data[lib] ++;
        }
        else{
            if(possible_fake_data.find(bam_name) == possible_fake_data.end())
                possible_fake_data[bam_name] = 1;
            else
                possible_fake_data[bam_name] ++;
        }

        // region between begin and last
        // This stores the exact same counts yet again. Still not sure why.
        if(CN_lib == 1){
            if(nread_FR.find(lib) == nread_FR.end())
                nread_FR[lib] = 1;
            else
                nread_FR[lib] ++;
        }
        else{
            if(nread_FR.find(bam_name) == nread_FR.end())
                nread_FR[bam_name] = 1;
            else
                nread_FR[bam_name] ++;
        }
    }
    // mapQual is part of the bam2cfg input. I infer it is a perlibrary mapping quality cutoff
  if(mapQual.find(lib) != mapQual.end()){
      if(b->core.qual <= mapQual[lib])
          return;
  }
  else{
      //here filter out if mapping quality is less than or equal to the min_map_qual.
      //Note that this doesn't make sense as a cutoff of 0 would still exclude reads with qual 0
      if(b->core.qual <= min_map_qual)
          return;
  }

    if(strcmp(bam_header->target_name[b->core.tid],"*")==0) // need to figure out how to compare a char and int //#ignore reads that failed to associate with a reference
        return;

    if(b->core.flag == 0)
        return; // return fragment reads

    //This is quite possibly the goofiest way to check core.flag == 32 that I know of.
    //Perhaps there are other flags I'm not considering?
  if((transchr_rearrange && b->core.flag < 32) || b->core.flag >=64) // only care flag 32 for CTX
        return;
    // for long insert
    // Mate pair libraries have different expected orientations so adjust
    // Also, aligner COULD have marked (if it was maq) that reads had abnormally large or small insert sizes
    // Remark based on BD options
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
    // This makes FF and RR the same thing
    if(b->core.flag == 8)
        b->core.flag = 1;
    // for normally mapped pairs, ignore if TOO far apart. Interesting... Isn't this exactly what we're trying to find?
    // For some reason doesn't apply to SmithWaterman MAQ reads. Don't think we have that info with BWA
    if(b->core.flag < 32 && abs(b->core.isize)>max_sd) // skip read pairs mapped too distantly on the same chromosome
        return;

    //count reads mapped by SW, FR and RF reads, but only if normal_switch is true
    //normal_switch is set to 1 as soon as reads are accumulated for dumping to fastq??? Not sure on this. Happens later in this function
    //I suspect this is to include those reads in the fastq dump for assembly!
    if(b->core.flag == 18 || b->core.flag == 20 || b->core.flag == 130){
        if(*normal_switch == 1 && b->core.isize > 0){
            (*nnormal_reads)++;
        }
        return;
    }

    if(*normal_switch == 1){
        *ntotal_nucleotides += b->core.l_qseq;
        *max_readlen = (*max_readlen < b->core.l_qseq) ? b->core.l_qseq : *max_readlen;
    }
    //This appears to test that you've exited a window after your first abnormal read by either reading off the chromosome or exiting the the window
    // d appears to be 1e8 at max (seems big), 50 at minimum or the smallest mean - readlen*2 for a given library
    int do_break = (int(b->core.tid) != lasts || int(b->core.pos) - lastc > d)?1:0;

    if(do_break){ // breakpoint in the assembly
        float seq_coverage = *ntotal_nucleotides/float(lastc - beginc + 1 + *max_readlen);
        if(lastc - beginc > min_len && seq_coverage < seq_coverage_lim){ // skip short/unreliable flnaking supporting regions
            // register reliable region and supporting reads across gaps
            int k = (*reg_idx) ++;  //assign an id to this region
            reg_name[k].push_back(begins); //chromosome
            reg_name[k].push_back(beginc); //starting coord
            reg_name[k].push_back(lastc);   //last coordinate
            reg_name[k].push_back(*nnormal_reads); //number of normal reads

            // never been to possible_fake in this turn, record ROI; or else the possible fake is not the fake, but the true one, doesn't need to record it in ROI, previous regions were recorded already
            // record nread_ROI
            // track the number of reads from each library for the region
            for(map<string, uint32_t>::const_iterator nread_ROI_it = nread_ROI.begin(); nread_ROI_it != nread_ROI.end(); nread_ROI_it ++){
                string lib_ = (*nread_ROI_it).first;
                if(read_count_ROI_map.find(k) == read_count_ROI_map.end()){
                    read_count_ROI_map[k][lib_] = nread_ROI[lib_];
                }
                else if(read_count_ROI_map[k].find(lib_) == read_count_ROI_map[k].end()){
                    read_count_ROI_map[k][lib_] = nread_ROI[lib_];
                }
                else{
                    read_count_ROI_map[k][lib_] += nread_ROI[lib_];
                }
            }

            // compute nread_FR and record it
            // track number of FR reads from the region.
            // From earlier, these numbers seem like they shoudl be the same unless they are being added to in multiple places
            for(map<string, uint32_t>::const_iterator nread_FR_it = nread_FR.begin(); nread_FR_it != nread_FR.end(); nread_FR_it ++){
                string lib_ = (*nread_FR_it).first;
                if(read_count_ROI_map[k].find(lib_) == read_count_ROI_map[k].end())
                    read_count_FR_map[k][lib_] = (*nread_FR_it).second;
                else{
                    // based on the variable names the ROI reads should contain ARPs as well as standard FR reads
                    // not sure where that would happen still
                    uint32_t diff = (*nread_FR_it).second - read_count_ROI_map[k][lib_];
                    if(diff < 0)
                        cout << "wrong, the subtraction is negative";
                    else if(diff >= 0){
                        read_count_FR_map[k][lib_] = diff;//(*nread_FR_it).second - read_count_ROI_map[k][lib_];
                    }
                    //else
                    //    cout << "lib only exist in ROI rather than FR";
                }

            }

            // reg_seq is an array of arrays of information about the reads. The first value is the read name.
            // This adds the region id to an array of region ids
            vector<vector<string> > p;
            for(vector<vector<string> >::const_iterator it_reg_seq = reg_seq.begin(); it_reg_seq != reg_seq.end(); it_reg_seq ++){
                p.push_back(*it_reg_seq);
                string s = (*it_reg_seq)[0];
                read[s].push_back(k);
            }

            regs[k] = p;    //store the array of reads at the region
            (*idx_buff)++; //increment tracking of number of regions in buffer??? Not quite sure if this is what idx_buff is
            if(*idx_buff > buffer_size){
                //flush buffer by building connection
                buildConnection(
                    bdancer,
                    read,
                    reg_name,
                    regs,
                    x_readcounts,
                    reference_len,
                    fisher,
                    min_read_pair,
                    dump_BED,
                    *max_readlen,
                    prefix_fastq,
                    ReadsOut,
                    SVtype,
                    mean_insertsize,
                    bam_header,
                    read_density,
                    CN_lib,
                    maps,
                    print_AF,
                    score_threshold,
                    libmaps
                );
                *idx_buff = 0;
            }
        }
        else{
          // possible fake
            // restore ROI for copy number since lastc will be cleared, but the new node has not been registered
            // possible fake is off, never gone to possible fake before, save the ROI
            // I don't understand exactly what this is doing. It is only hitting here to store the info if flanking region is too short or the coverage is too high
            // It appears to be used to pull in nearby neighboring regions to the last region identified if the distance between them is too short
            if(/* *possible_fake == 1 &&*/ *reg_idx >= 1){
                 for(map<string, uint32_t>::const_iterator possible_fake_data_it = possible_fake_data.begin(); possible_fake_data_it != possible_fake_data.end(); possible_fake_data_it ++){
                     string lib_ = (*possible_fake_data_it).first;
                     if(read_count_ROI_map.find(*reg_idx-1) == read_count_ROI_map.end()){
                         read_count_ROI_map[*reg_idx-1][lib_] = possible_fake_data[lib_];
                         //read_count_ROI_debug[*reg_idx-1][lib_] = possible_fake_data_debug[lib_];
                     }
                     else if(read_count_ROI_map[*reg_idx-1].find(lib_) == read_count_ROI_map[*reg_idx-1].end()){
                         read_count_ROI_map[*reg_idx-1][lib_] = possible_fake_data[lib_];
                         //read_count_ROI_debug[*reg_idx-1][lib_] = possible_fake_data_debug[lib_];
                     }
                     else{
                         read_count_ROI_map[*reg_idx-1][lib_] += possible_fake_data[lib_];
                     }
                 }
            }

            // remove any reads that are linking the last region with this new, merged in region
            if(reg_seq.size()>0){
                for(vector<vector<string> >::const_iterator it_reg_seq = reg_seq.begin(); it_reg_seq != reg_seq.end(); it_reg_seq ++){
                    string s= (*it_reg_seq)[0];
                    if(read.find(s) != read.end())
                        read.erase(read.find(s));
                }
            }
        }
        // clear out this node
        begins = int(b->core.tid);
        beginc = int(b->core.pos);
        reg_seq.clear();
        *normal_switch = 0;
        *nnormal_reads = 0;
        *max_readlen = 0;
        *ntotal_nucleotides = 0;

        // clear possible fake data
        possible_fake_data.clear();
        nread_ROI.clear();
        // clear FR
        nread_FR.clear();
    }

    // deal with the name string
    // this drops any trailing /1 or /2 strings on the read. Probably not necessary if the BAM is well formatted (but I'm not certain).
    string qname_tmp = bam1_qname(b);
    size_t found1 = qname_tmp.rfind("/1");
    size_t found2 = qname_tmp.rfind("/2");
    if(found1 != string::npos || found2 != string::npos){
        size_t found = (found1 == string::npos) ? found2 : found1;
        qname_tmp.replace(found,2,"");
    }
    //bam1_qname(b) = qname_tmp; // this might be quite hard to implement in samtools; ///////////////////////

    string seq = get_string(bam1_seq(b), b->core.l_qseq);
    string basequal = get_string_qual(bam1_qual(b), b->core.l_qseq);
    if(! prefix_fastq.empty() && ! seq.empty() && ! basequal.empty()){
        //This stores all the read info as a vector of strings. WTF that can't be efficient.
        vector<string> tmp_reg_seq;
        tmp_reg_seq.push_back(qname_tmp);
        tmp_reg_seq.push_back(itos(int(b->core.tid)));
        tmp_reg_seq.push_back(itos(int(b->core.pos)));
        tmp_reg_seq.push_back(ori);
        tmp_reg_seq.push_back(itos(int(b->core.isize)));
        tmp_reg_seq.push_back(itos(int(b->core.flag)));
        tmp_reg_seq.push_back(itos(int(b->core.qual)));
        tmp_reg_seq.push_back(itos(int(b->core.l_qseq)));
        tmp_reg_seq.push_back(lib);
        tmp_reg_seq.push_back(seq);
        tmp_reg_seq.push_back(basequal);
        reg_seq.push_back(tmp_reg_seq); // store each read in the region_sequence buffer
    }
    else{
        vector<string> tmp_reg_seq;
        tmp_reg_seq.push_back(qname_tmp);
        tmp_reg_seq.push_back(itos(int(b->core.tid)));
        tmp_reg_seq.push_back(itos(int(b->core.pos)));
        tmp_reg_seq.push_back(ori);
        tmp_reg_seq.push_back(itos(int(b->core.isize)));
        tmp_reg_seq.push_back(itos(int(b->core.flag)));
        tmp_reg_seq.push_back(itos(int(b->core.qual)));
        tmp_reg_seq.push_back(itos(int(b->core.l_qseq)));
        tmp_reg_seq.push_back(lib);
        reg_seq.push_back(tmp_reg_seq);
    }
    //If we just added the first read, flip the flag that lets us collect all reads
    if(reg_seq.size() == 1)
        *normal_switch = 1;
    lasts = int(b->core.tid);
    lastc = int(b->core.pos);
    nread_ROI.clear(); //Not sure why this is cleared here...
}

// pair up reads and print out results (SV estimation)
void buildConnection(
    BreakDancerData& bdancer,
    map<string, vector<int> > &read,
    map<int, vector<int> > &reg_name,
    map<int, vector<vector<string> > > &regs,
    map<uint32_t, map<string,int> > &x_readcounts,
    uint32_t reference_len,
    int fisher,
    int min_read_pair,
    string dump_BED,
    int max_readlen,
    string prefix_fastq,
    map<string, string> &ReadsOut,
    map<string, string> &SVtype,
    map<string, float> &mean_insertsize,
    bam_header_t* bam_header,
    map<string, float> &read_density,
    int CN_lib,
    vector<string> maps, // FIXME: should be constref
    int print_AF,
    int score_threshold,
    map<string, string> libmaps // FIXME: should be constref
    )
{
    map<int, map<string, uint32_t> >& read_count_ROI_map = bdancer.read_count_ROI_map;
    map<int, map<string, uint32_t> >& read_count_FR_map = bdancer.read_count_FR_map;

    // build connections
    // find paired regions that are supported by paired reads
    //warn("-- link regions\n");
    map<int, map<int, int> > clink;
    map<string,vector<int> >::iterator ii_read;
    //read is a map of readnames, each is associated with a vector of region ids
    // wtf is this using a vector? How would we ever have more than two regions? Multi-mapping?
    for(ii_read = read.begin(); ii_read != read.end(); ii_read++){
        // test
        //string tmp_str = (*ii_read).first;
        vector<int> p = (*ii_read).second;
        assert( p.size() < 3);
        if(p.size() != 2) // skip singleton read (non read pairs)
            continue;
        //cout << tmp_str << "\t" << p[0] << "\t" << p[1] << endl;
        //track the number of links between two nodes
        if(clink.find(p[0]) != clink.end() && clink[p[0]].find(p[1]) != clink[p[0]].end()){
            ++clink[p[0]][p[1]];
            ++clink[p[1]][p[0]];
        }
        else{
            clink[p[0]][p[1]] = 1;
            clink[p[1]][p[0]] = 1;
        }
        //cout << tmp_str << endl;
        //cout << p[0] << "\t" << p[1] << "\t" << link[p[0]][p[1]] << endl;
    }
    // segregate graph, find nodes that have connections
    map<int,int> free_nodes;
    map<int, map<int, int> >::iterator ii_clink;
    //  int tmp_size = clink.size();
    vector<int> s0_vec;
    //  int tmp_read_size = read.size();
    //cout << tmp_read_size << endl;

    // Grab the first region for every link
    for(ii_clink = clink.begin(); ii_clink != clink.end(); ii_clink++){
        s0_vec.push_back((*ii_clink).first);
        //cout << ",,,,," << (*ii_clink).first << endl;
        //      map<int,int> tmp_clink = (*ii_clink).second;
        /*int s1, value;
         for(map<int,int>::iterator ii_tmp_clink = tmp_clink.begin(); ii_tmp_clink != tmp_clink.end(); ii_tmp_clink++){
         s1 = (*ii_tmp_clink).first;
         //cout << ",,,," << s1 << endl;
         value = (*ii_tmp_clink).second;
         }*/
    }
    for(vector<int>::iterator ii_s0_vec = s0_vec.begin(); ii_s0_vec != s0_vec.end(); ii_s0_vec ++){
        int s0 = *ii_s0_vec;
        //cout << ",,,,," << s0 << endl;
        // assert( clink.find(s0) != clink.end() ); THIS ASSERT TRIPS
        if(clink.find(s0) == clink.end())
            continue;
        // construct a subgraph
        vector<int> tails;
        tails.push_back(s0);
        while(tails.size() > 0){
            vector<int> newtails;
            vector<int>::iterator it_tails;
            for(it_tails = tails.begin(); it_tails != tails.end(); it_tails ++){
                int tail = *it_tails;
                //cout << ",,,," << tail << endl;
                //assert(clink.find(*it_tails) != clink.end()); THIS ASSERT TRIPS
                if(clink.find(*it_tails) == clink.end())
                    continue;
                assert(reg_name.find(*it_tails) != reg_name.end());
                if(reg_name.find(*it_tails) == reg_name.end())
                    continue;
                vector<int> s1s; //accumulate all linked nodes for a  single node
                for(map<int, int>::iterator ii_clink_ = clink[tail].begin(); ii_clink_ != clink[tail].end(); ii_clink_++){
                    s1s.push_back((*ii_clink_).first);
                    //cout << ",,," << (*ii_clink_).first << endl;
                }
                for(vector<int>::iterator ii_s1s = s1s.begin(); ii_s1s != s1s.end(); ii_s1s++){
                    int s1 = *ii_s1s;
                    //cout << ",," << s1 << endl;
                    vector<string> free_reads;
                    map<int,map<int,int> > nodepair;
                    int nlinks = clink[tail][s1];
                    if(nlinks<min_read_pair) // require sufficient number of pairs
                        continue;
                    assert(nodepair.find(s1) == nodepair.end()); // a node only appear once in a pair
                    if(nodepair.find(s1) != nodepair.end()) // a node only appear once in a pair
                        continue;
                    assert(reg_name.find(s1) != reg_name.end()); // a node must be defined
                    if(reg_name.find(s1) == reg_name.end()) // a node must be defined
                        continue;
                    nodepair[tail][s1] = clink[tail][s1];
                    nodepair[s1][tail] = clink[s1][tail];
                    if(clink[tail].find(s1)!=clink[tail].end()){
                        clink[tail].erase(clink[tail].find(s1));    // use a link only once
                    }
                    //NOTE it is entirely possible that tail and s1 are the same.
                    if(tail != s1){
                        if(clink[s1].find(tail)!=clink[s1].end()){
                            clink[s1].erase(clink[s1].find(tail));
                        }
                    }
                    newtails.push_back(s1);

                    // analysis a nodepair
                    vector<int> snodes;
                    for(map<int,map<int,int> >::iterator ii_nodepair = nodepair.begin(); ii_nodepair != nodepair.end(); ii_nodepair ++){
                        snodes.push_back((*ii_nodepair).first);
                        //cout << "," << (*ii_nodepair).first << endl;
                    }
                    assert(snodes.size() < 3);
                    // track node1 and node2 as nodes that could potentially be freed
                    int node1 = snodes[0];
                    int node2;
                    if(snodes.size() == 1)
                        node2 = snodes[0];
                    else
                        node2 = snodes[1];
                    if(nodepair.find(node1) == nodepair.end() || nodepair[node1].find(node2) == nodepair[node1].end())
                        continue;

                    int nread_pairs = 0;
                    map<string,vector<string> > read_pair; //unpaired reads
                    map<string,int> type; // number of readpairs per each type/flag
                    map<string,map<string,int> > type_library_readcount; // number of readpairs per each type/flag (first key) then library (second key)
                    vector<map<string,int> > type_orient_counts; //vector of readcounts for each type/flag
                    map<string,map<string,int> > type_library_meanspan; //average ISIZE from BAM records
                    vector<vector<string> > support_reads; //reads supporting the SV
                    for(vector<int>::iterator ii_snodes = snodes.begin(); ii_snodes < snodes.end(); ii_snodes++){
                        int node = *ii_snodes;
                        //cout << node << endl;
                        map<string,int> orient_count; // number of reads per each orientation ('+' or '-')
                        vector<vector<string> > nonsupportives; // reads not supporting this SV
                        //debug
                        //int regs_size = regs[node].size();
                        //NOTE regs contains an array of information about the reads supporting the region (info is stored as a string array)
                        for(vector<vector<string> >::iterator ii_regs = regs[node].begin(); ii_regs != regs[node].end(); ii_regs++){
                            vector<string> y = *ii_regs;
                            //cout << y[3] << "\t" << y[0] << "\t" << y[2] << "\t" << orient_count[y[3]] << endl;
                            //skip things where the read name is no longer in our list of read names
                            //WHY ARE THESE CHECKS EVERYWHERE
                            if(read.find(y[0]) == read.end())
                                continue;
                            // initialize orient_count
                            // y[3] is the orientation. This is stored as a string value or - or +
                            if(orient_count.find(y[3]) == orient_count.end())
                                orient_count[y[3]] = 1;
                            else
                                orient_count[y[3]]++;

                            //START HERE
                            if(read_pair.find(y[0]) == read_pair.end()){
                                read_pair[y[0]] = y;
                                nonsupportives.push_back(y);
                                //cout << y[3] << "\t" << y[0] << "\t" << y[2] << "\t" << orient_count[y[3]] << endl;
                            }
                            else{
                                // see if initialized 'type' or not
                                // y[5] is our bastardized "flag" describing the pair orientation
                                if(type.find(y[5]) != type.end())
                                    type[y[5]]++;
                                else
                                    type[y[5]] = 1;
                                // see if initialized 'type_library_readcount' or not
                                // as you might guess, y[8] is the lib name
                                if(type_library_readcount.find(y[5]) != type_library_readcount.end() && type_library_readcount[y[5]].find(y[8]) != type_library_readcount[y[5]].end())
                                    type_library_readcount[y[5]][y[8]]++;
                                else
                                    type_library_readcount[y[5]][y[8]] = 1;

                                // y[4] is the insert size (as a string of course)
                                int y4_tmp = atoi(y[4].c_str());
                                if(type_library_meanspan.find(y[5]) != type_library_meanspan.end() && type_library_meanspan[y[5]].find(y[8]) != type_library_meanspan[y[5]].end()){
                                    type_library_meanspan[y[5]][y[8]]+=abs(y4_tmp);
                                }
                                else
                                    type_library_meanspan[y[5]][y[8]] = abs(y4_tmp);
                                nread_pairs++;
                                free_reads.push_back(y[0]);
                                //cout << y[0] << endl;
                                support_reads.push_back(y);
                                support_reads.push_back(read_pair[y[0]]);
                                if(read_pair.find(y[0])!=read_pair.end()){
                                    read_pair.erase(read_pair.find(y[0]));
                                }
                            }
                        }
                        regs[node] = nonsupportives;
                        type_orient_counts.push_back(orient_count);
                    }

                    //clean out supportive reads since the first read of every pair is stored in nonsupportives
                    //seems like it would be smarter to only populate nonsupportives after we determine the supportives...
                    //I think this must not be done because you don't know if the read pairs will occur on the other node or not
                    //so you build your list for each node and then go back and clean it up after the fact
                    //if you tracked which node each read came from then you could just reassign after the fact
                    for(vector<int>::iterator ii_snodes = snodes.begin(); ii_snodes != snodes.end(); ii_snodes++){
                        int node = *ii_snodes;
                        vector<vector<string> > nonsupportives;
                        for(vector<vector<string> >::iterator ii_regs = regs[node].begin(); ii_regs != regs[node].end(); ii_regs++){
                            vector<string> y = *ii_regs;
                            if(read_pair.find(y[0]) == read_pair.end())
                                continue;
                            nonsupportives.push_back(y);
                        }
                        regs[node] = nonsupportives;
                    }

                    //float score;//don't know if float; no usage actually
                    //int bestIndelSize;//don't know if int; no usage actually
                    // START HERE
                    if(nread_pairs >= min_read_pair){
                        map<string, int> diffspans;
                        map<string, string> sptypes;
                        string flag = choose_sv_flag(nread_pairs, type);
                        if(type[flag] >= min_read_pair) {
                            // print out result
                            int sv_chr1 = -1, sv_pos1 = 0, sv_chr2 = -1, sv_pos2 = 0;
                            string sv_ori1, sv_ori2;
                            int normal_rp;

                            int first_node = 0;
                            map<string, uint32_t> read_count;
                            // find inner most positions
                            for(vector<int>::iterator ii_snodes = snodes.begin(); ii_snodes != snodes.end(); ii_snodes ++){
                                int node = *ii_snodes;
                                //cout << node << "\t";
                                int chr = reg_name[node][0];
                                int start = reg_name[node][1];
                                int end = reg_name[node][2];
                                int nrp = reg_name[node][3];

                                //cout << " " << node << "\t" << start << "\t" << end << endl;
                                map<string,int> ori_readcount = type_orient_counts.front();
                                if(type_orient_counts.size()!=0){
                                    type_orient_counts.erase(type_orient_counts.begin());
                                }
                                if(sv_chr1 != -1 && sv_chr2 != -1){
                                    if(flag.compare("4") == 0)
                                        sv_pos2 = end + max_readlen - 5;
                                    else if(flag.compare("1") == 0){
                                        sv_pos1 = sv_pos2;
                                        sv_pos2 = end + max_readlen - 5;
                                    }
                                    else if(flag.compare("8") == 0)
                                        sv_pos2 = start;
                                    else{
                                        sv_pos1 = sv_pos2;
                                        sv_pos2 = start;
                                    }
                                    sv_chr1 = sv_chr2;
                                    //sv_pos1 = sv_pos2;
                                    sv_chr2 = chr;
                                    //sv_pos2 = start;
                                    string sv_ori2_tmp1 = "0";
                                    string sv_ori2_tmp2 = "0";
                                    if(ori_readcount.find("+") != ori_readcount.end())
                                        //sprintf(sv_ori2_tmp1, "%s", ori_readcount["+"]);
                                        sv_ori2_tmp1 = itos(ori_readcount["+"]);
                                    if(ori_readcount.find("-") != ori_readcount.end())
                                        //sprintf(sv_ori2_tmp2, "%s", ori_readcount["-"]);
                                        sv_ori2_tmp2 = itos(ori_readcount["-"]);
                                    sv_ori2 = sv_ori2_tmp1.append("+").append(sv_ori2_tmp2).append("-");

                                    // add up the read number
                                    for(int i_node = first_node; i_node < node; i_node++){
                                        map<string, uint32_t> read_count_ROI_map_second = read_count_ROI_map[i_node];
                                        map<string, uint32_t> read_count_FR_map_second = read_count_FR_map[i_node];
                                        //cout << "\nROI:" << reg_name[i_node][2] << "\t" << reg_name[i_node + 1][1] << "\n";
                                        for(map<string, uint32_t>::iterator read_count_ROI_map_second_it = read_count_ROI_map_second.begin(); read_count_ROI_map_second_it != read_count_ROI_map_second.end(); read_count_ROI_map_second_it ++){
                                            string lib = (*read_count_ROI_map_second_it).first;
                                            if(read_count.find(lib) == read_count.end())
                                                read_count[lib] = read_count_ROI_map[i_node][lib];
                                            else
                                                read_count[lib] += read_count_ROI_map[i_node][lib];

                                            //for(map<string, int>::iterator i_debug = read_count_ROI_debug[i_node][lib].begin(); i_debug != read_count_ROI_debug[i_node][lib].end(); i_debug++){
                                            //      cout << (*i_debug).first << "\t" << (*i_debug).second << "\n";
                                            //}
                                        }

                                        // flanking region doesn't contain the first node
                                        if(i_node == first_node)
                                            continue;
                                        //cout << "\nFR:" << reg_name[i_node][1] << "\t" << reg_name[i_node][2] << "\n";
                                        for(map<string, uint32_t>::iterator read_count_FR_map_second_it = read_count_FR_map_second.begin(); read_count_FR_map_second_it != read_count_FR_map_second.end(); read_count_FR_map_second_it ++){
                                            string lib = (*read_count_FR_map_second_it).first;
                                            if(read_count.find(lib) == read_count.end())
                                                read_count[lib] = read_count_FR_map[i_node][lib];
                                            else
                                                read_count[lib] += read_count_FR_map[i_node][lib];

                                            //  for(map<string, int>::iterator i_debug = read_count_FR_debug[i_node][lib].begin(); i_debug != read_count_FR_debug[i_node][lib].end(); i_debug++){
                                            //      cout << (*i_debug).first << "\t" << (*i_debug).second << "\n";
                                            //  }
                                        }
                                    }
                                }
                                else {
                                    first_node = node;
                                    sv_chr1 = chr;
                                    sv_chr2 = chr;
                                    sv_pos1 = start;
                                    sv_pos2 = end;
                                    string sv_ori2_tmp1 = "0";
                                    string sv_ori2_tmp2 = "0";
                                    if(ori_readcount.find("+") != ori_readcount.end())
                                        sv_ori2_tmp1 = itos(ori_readcount["+"]);
                                    if(ori_readcount.find("-") != ori_readcount.end())
                                        sv_ori2_tmp2 = itos(ori_readcount["-"]);
                                    sv_ori1 = sv_ori2_tmp1.append("+").append(sv_ori2_tmp2).append("-");
                                    sv_ori2 = sv_ori1;
                                    normal_rp = nrp;

                                    // get the read number for this region: which is not very accurate now. Need to count all of the flags.
                                    //read_count_ROI = normal_rp + read_count_intra;
                                }
                            }
                            // cout << "\n";
                            //cout << sv_pos1 + 1 << endl;

                            // get the copy_number from read_count
                            map<string, float> copy_number;
                            float copy_number_sum = 0;
                            for(map<string, uint32_t>:: iterator read_count_it = read_count.begin(); read_count_it != read_count.end(); read_count_it ++){
                                string lib = (*read_count_it).first;
                                copy_number[lib] = (float)((*read_count_it).second)/((float)read_density[lib] * float(sv_pos2 - sv_pos1))*2;
                                copy_number_sum += copy_number[lib];
                                //cout << lib << "\t" << (*read_count_it).second << "\t" << read_density[lib] << "\t" << sv_pos2-sv_pos1 << "\t" << copy_number[lib] << endl;

                            }
                            copy_number_sum /= (2.0*(float)read_count.size());

                            if(flag.compare("4") && flag.compare("8") && sv_pos1 + max_readlen - 5 < sv_pos2)
                                sv_pos1 += max_readlen - 5; // apply extra padding to the start coordinates

                            // deal with directly flag, rather than for each 'fl', since flag is already known, and diffspans and sptypes are only used for flag;
                            string sptype;
                            float diffspan = 0;
                            //debug
                            //int tmp_size_tlr = type_library_readcount[fl].size();
                            if(CN_lib == 1){
                                for(map<string,int>::iterator ii_type_lib_rc = type_library_readcount[flag].begin(); ii_type_lib_rc != type_library_readcount[flag].end(); ii_type_lib_rc ++){
                                    string sp = (*ii_type_lib_rc).first;
                                    // intialize to be zero, in case of no library, or DEL, or ITX.

                                    string copy_number_str = "NA";
                                    if(flag.compare("32")!=0){
                                        float copy_number_ = 0;

                                        if(copy_number.find(sp) != copy_number.end()){
                                            copy_number_ = copy_number[sp];
                                            stringstream sstr;
                                            sstr << fixed;
                                            sstr << setprecision(2) << copy_number_;
                                            copy_number_str = sstr.str();
                                        }
                                    }
                                    //string str_num_tmp;
                                    //sprintf(str_num_tmp, "%s", (*ii_type_lib_rc).second);
                                    if(!sptype.empty())
                                        sptype += ":" +  sp + "|" + itos((*ii_type_lib_rc).second) + "," + copy_number_str;
                                    else
                                        sptype = sp + "|" + itos((*ii_type_lib_rc).second) + "," + copy_number_str;
                                    diffspan += float(type_library_meanspan[flag][sp]) - float(type_library_readcount[flag][sp])*mean_insertsize[sp];
                                }
                            } // do lib for copy number and support reads
                            else{
                                map<string, int> type_bam_readcount;
                                for(map<string, int>::iterator ii_type_lib_rc = type_library_readcount[flag].begin(); ii_type_lib_rc != type_library_readcount[flag].end(); ii_type_lib_rc ++){
                                    string sp = (*ii_type_lib_rc).first;
                                    if(libmaps.find(sp) != libmaps.end()){
                                        string sp_bam = libmaps[sp];
                                        if(type_bam_readcount.find(sp_bam)!= type_bam_readcount.end()){
                                            type_bam_readcount[sp_bam] += (*ii_type_lib_rc).second;
                                        }
                                        else{
                                            type_bam_readcount[sp_bam] = (*ii_type_lib_rc).second;
                                        }
                                    }
                                    diffspan += float(type_library_meanspan[flag][sp]) - float(type_library_readcount[flag][sp])*mean_insertsize[sp];
                                }
                                for(map<string, int>::iterator ii_type_bam_rc = type_bam_readcount.begin(); ii_type_bam_rc != type_bam_readcount.end(); ii_type_bam_rc ++){
                                    string sp = (*ii_type_bam_rc).first;
                                    if(!sptype.empty())
                                        sptype += ":" + sp + "|" + itos((*ii_type_bam_rc).second);
                                    else
                                        sptype = sp + "|" + itos((*ii_type_bam_rc).second);
                                }
                                if(sptype.length() == 0){
                                    sptype = "NA";
                                }
                            } // do bam for support reads; copy number will be done later

                            //debug
                            //int tmp_tlm = type_library_meanspan[fl][sp];
                            //int tmp_tlr = type_library_readcount[fl][sp];
                            //int tmp_mi = mean_insertsize[sp];
                            diffspans[flag] = int(diffspan/float(type[flag]) + 0.5);
                            sptypes[flag] = sptype;


                            int flag_int = atoi(flag.c_str());
                            real_type LogPvalue = ComputeProbScore(snodes, type_library_readcount[flag], uint32_t(flag_int), x_readcounts, reference_len, fisher, reg_name);
                            real_type PhredQ_tmp = -10*LogPvalue/log(10);
                            int PhredQ = PhredQ_tmp>99 ? 99:int(PhredQ_tmp+0.5);
                            //float AF = float(type[flag])/float(type[flag]+normal_rp);
                            float AF = 1 - copy_number_sum;


                            string SVT = SVtype.find(flag)==SVtype.end()?"UN":SVtype[flag]; // UN stands for unknown
                            // make the coordinates with base 1
                            sv_pos1 = sv_pos1 + 1;
                            sv_pos2 = sv_pos2 + 1;
                            /*if(sv_pos1 == 17535916){
                             int k = 0;
                             k++;
                             }*/
                            //cout << bam_header->target_name[sv_chr1] << "\t" << sv_pos1 << "\t"  << sv_ori1 << "\t" << bam_header->target_name[sv_chr2] << "\t" << sv_pos2 << "\t" << sv_ori2 << "\t" << SVT << "\t" << diffspans[flag] << "\t" << PhredQ << "\t" << type[flag] << "\t" << sptypes[flag] << "\t" << AF << "\t" << version << "\t" << options << endl;
                            if(PhredQ > score_threshold){
                                cout << bam_header->target_name[sv_chr1] << "\t" << sv_pos1 << "\t"  << sv_ori1 << "\t" << bam_header->target_name[sv_chr2] << "\t" << sv_pos2 << "\t" << sv_ori2 << "\t" << SVT << "\t" << diffspans[flag] << "\t" << PhredQ << "\t" << type[flag] << "\t" << sptypes[flag];// << endl;
                                //printf("%d\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%s\t%.2f\t%s\t%s\n",sv_chr1,sv_pos1,sv_ori1,sv_chr2,sv_pos2,sv_ori2,SVT,diffspans[flag],PhredQ,type[flag],sptypes[flag],AF,version,options);// version and options should be figured out. Should do it later.
                                if(print_AF == 1)
                                    cout <<  "\t" << AF;
                                if(CN_lib == 0 && flag.compare("32")!=0){
                                    for(vector<string>::const_iterator iter = maps.begin(); iter != maps.end(); ++iter) {
                                        if(copy_number.find(*iter) == copy_number.end())
                                            cout << "\tNA";
                                        else {
                                            cout << "\t";
                                            cout << fixed;
                                            cout << setprecision(2) << copy_number[*iter];
                                        }
                                    }
                                }
                                cout << endl;



                                if(!prefix_fastq.empty()){ // print out supporting read pairs
                                    write_fastq_for_flag(flag, support_reads, ReadsOut);
                                }

                                if(!dump_BED.empty()){  // print out SV and supporting reads in BED format
                                    ofstream fh_BED;
                                    fh_BED.open(dump_BED.c_str(), ofstream::app);

                                    string trackname(bam_header->target_name[sv_chr1]);
                                    trackname = trackname.append("_").append(itos(sv_pos1)).append("_").append(SVT).append("_").append(itos(diffspans[flag]));
                                    //string fh_BED_tmp = "track name=".append(trackname).append("\tdescription=\"BreakDancer ").append(itoa(sv_chr1)).append(" ").append(itoa(sv_pos1)).append(" ").append(SVT).append(" ").append(itoa(diffspans[flag]));
                                    fh_BED << "track name=" << trackname << "\tdescription=\"BreakDancer" << " " << bam_header->target_name[sv_chr1] << " " << sv_pos1 << " " << SVT << " " << diffspans[flag] << "\"\tuseScore=0\n";// fh_BED is a file handle of BED
                                    for(vector<vector<string> >::iterator ii_support_reads = support_reads.begin(); ii_support_reads != support_reads.end(); ii_support_reads ++){
                                        vector<string> y = *ii_support_reads;
                                        if(y.size() < 8 || y[5].compare(flag))
                                            continue;
                                        if(y[1].find("chr")!=string::npos)
                                            y[1].erase(y[1].find("chr"),3);
                                        int y1_int = atoi(y[1].c_str());
                                        int y2_int, y7_int;
                                        y2_int = atoi(y[2].c_str()) + 1;
                                        y7_int = atoi(y[7].c_str());
                                        int aln_end = y2_int + y7_int - 1;
                                        string color = y[3].compare("+")?"0,0,255":"255,0,0";
                                        //fh_BED << "chr" << bam_header->target_name[y1_int] << "\t" << y2_int << "\t" << aln_end << "\t1\t" << y[0] << "\t" << y[3] << "\t" << y[4] << "\t" << y2_int << "\t" << aln_end << "\t" << color << "\n";//sprintf(fh_BED, "chr%s\t%s\t%s\t%s\t1\t%s\t%s\t%s\t%d\t%s\n",y[1],y[2],aln_end,y[0],y[3],y[4],y[2],aln_end,color);
                                        int aln_score = atoi(y[6].c_str()) * 10;
                                        fh_BED << "chr" << bam_header->target_name[y1_int] << "\t" << y2_int << "\t" << aln_end << "\t" << y[0] << "|" << y[8] << "\t" << aln_score << "\t" << y[3] << "\t" << y2_int << "\t" << aln_end << "\t" << color << "\n";
                                    }
                                    fh_BED.close();
                                }
                            }
                        }
                        // free reads
                        for(vector<string>::iterator ii_free_reads = free_reads.begin(); ii_free_reads != free_reads.end(); ii_free_reads ++){
                            if(read.find(*ii_free_reads) != read.end()){
                                read.erase(read.find(*ii_free_reads));
                            }
                        }
                        //free_reads.clear();
                        //record list of nodes that can be potentially freed
                        free_nodes[node1] = 1;
                        free_nodes[node2] = 1;
                    }
                }
                if(clink.find(tail)!=clink.end()){
                    clink.erase(clink.find(tail));
                }
            }
            tails = newtails;
        }
    }

    // free nodes
    for(map<int,int>::iterator ii_free_nodes = free_nodes.begin(); ii_free_nodes != free_nodes.end(); ii_free_nodes++){
        // remove reads in the regions
        int node = (*ii_free_nodes).first;
        vector<vector<string> > reads = regs[node];
        if(reads.size() < unsigned(min_read_pair)){
            for(vector<vector<string> >::iterator ii_reads = reads.begin(); ii_reads != reads.end(); ii_reads++){
                vector<string> y = *ii_reads;
                string readname = y[0];
                //cout << readname << endl;
                if(read.find(readname)!=read.end()){
                    read.erase(read.find(readname));
                }
            }
            // remove regions
            if(regs.find(node) != regs.end()){
                regs.erase(regs.find(node));
            }
            if(reg_name.find(node) != reg_name.end()){
                reg_name.erase(reg_name.find(node));
            }
        }
    }
}


// estimate prior parameters (not being used now)
void EstimatePriorParameters(
    Options const& opts,
    map<string,string> &fmaps,
    map<string,string> &readgroup_library,
    map<string, float> &mean_insertsize,
    map<string, float> &std_insertsize,
    map<string,float> &uppercutoff,
    map<string,float> &lowercutoff,
    map<string,float> &readlens,
    map<string, string> &readgroup_platform
    )
{
    map<string,float> es_means;
    map<string,float> es_stds;
    map<string,float> es_readlens;
    map<string,float> es_uppercutoff;
    map<string,float> es_lowercutoff;
    map<string,vector<int> > insert_stat;
    map<string,vector<int> > readlen_stat;

    bam1_t *b = bam_init1();
    for(map<string,string>::const_iterator ii=fmaps.begin(); ii!=fmaps.end(); ++ii){
        string bam_name = (*ii).first;
        RegionLimitedBamReader reader(bam_name, opts.chr.c_str());
        bam_header_t* header = reader.header();
        while (reader.next(b) > 0) {
            string format = "sam";
            string alt = "";
            int same_tid = strcmp(
                header->target_name[b->core.tid],
                header->target_name[b->core.mtid]) == 0 ? 1:0;

            breakdancer::Read aln2(b, format, readgroup_platform, readgroup_library);
            vector<string> aln_return = AlnParser(b, format, alt, readgroup_platform, same_tid, opts.platform);
            assert(b->core.flag == aln2.bdflag);
            string ori = aln_return[1];
            string readgroup = aln_return[0];

            // analyze the bam file line by line
            string lib = readgroup.empty()?(*ii).second:readgroup_library[readgroup];// when multiple libraries are in a BAM file
            if(lib.empty())
                continue;
            // analyze 1 chromosome
            if(opts.chr != "0" && opts.chr != header->target_name[b->core.tid])
                continue;
            //if(readlen_stat.find(lib) == readlen_stat.end())    // don't need to issue a new stat
            //readlen_stat[lib] = ; // Statistics::Descriptive::Sparse->new() // don't need to issue a new stat
            readlen_stat[lib].push_back(b->core.isize);
            if(b->core.qual <= opts.min_map_qual)    // skip low quality mapped reads
                continue;
            if((b->core.flag != 18 && b->core.flag != 20) || b->core.isize <= 0)
                continue;
            //if(insert_stat.find(lib) == insert_stat.end())    // don't need to issue a new stat
            insert_stat[lib].push_back(b->core.isize);
        }
    }
    bam_destroy1(b);
    for(map<string,vector<int> >::const_iterator ii_readlen_stat = readlen_stat.begin(); ii_readlen_stat != readlen_stat.end(); ii_readlen_stat ++){
        string lib = (*ii_readlen_stat).first;
        //double res = accumulate(insert_stat[lib].begin(), insert_stat[lib].end(), 0);
        //double mean = res/insert_stat[lib].size();
        float mean_insert = mean(insert_stat[lib]);
        float std_insert = standard_deviation(insert_stat[lib],mean_insert);
        float uppercutoff = mean_insert + std_insert*opts.cut_sd;
        float lowercutoff = mean_insert - std_insert*opts.cut_sd;
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
    vector<int>::const_iterator ii_stat;
    for(ii_stat = stat.begin(); ii_stat < stat.end(); ii_stat++){
        all += *ii_stat;
    }
    return (float)all/(float)(stat.size());
}

// compute the standard deviation of a vector of int knowing the mean
float standard_deviation(vector<int> &stat, float mean){
    int all = 0;
    vector<int>::const_iterator ii_stat;
    for(ii_stat = stat.begin(); ii_stat < stat.end(); ii_stat ++){
        all += (*ii_stat)*(*ii_stat);
    }
    return sqrt((float)all/(float)(stat.size()) - mean*mean);
}

// compute the probability score
real_type ComputeProbScore(vector<int> &rnode, map<string,int> &rlibrary_readcount, uint32_t type, map<uint32_t, map<string,int> > &x_readcounts, uint32_t reference_len, int fisher, map<int, vector<int> > &reg_name) {
    // rnode, rlibrary_readcount, type
    int total_region_size = PutativeRegion(rnode, reg_name);

    real_type lambda;
    real_type logpvalue = 0.0;
    real_type err = 0.0;
    for(map<string,int>::const_iterator ii_rlibrary_readcount = rlibrary_readcount.begin(); ii_rlibrary_readcount != rlibrary_readcount.end(); ii_rlibrary_readcount ++){
        string lib = (*ii_rlibrary_readcount).first;
        // debug
        //int db_x_rc = x_readcounts[type][lib];
        lambda = real_type(total_region_size)* (real_type(x_readcounts[type][lib])/real_type(reference_len));
        lambda = max(real_type(1.0e-10), lambda);
        poisson_distribution<real_type> poisson(lambda);
        real_type tmp_a = log(cdf(complement(poisson, rlibrary_readcount[lib]))) - err;
        real_type tmp_b = logpvalue + tmp_a;
        err = (tmp_b - logpvalue) - tmp_a;
        _max_kahan_err = max(_max_kahan_err, err);
        logpvalue = tmp_b;
    }

    if(fisher && logpvalue < 0) {
        // Fisher's Method
        chi_squared chisq(2*rlibrary_readcount.size());
        try {
            real_type fisherP = cdf(complement(chisq, -2*logpvalue));
            logpvalue = fisherP > ZERO ? log(fisherP) : LZERO;
        } catch (const std::exception& e) {
            cerr << "chi squared problem: N=" << 2*rlibrary_readcount.size()
                << ", log(p)=" << logpvalue << ", -2*log(p) = " << -2*logpvalue << "\n";
        }
    }

    return logpvalue;
}

// putative region
int PutativeRegion(vector<int> &rnode, map<int,vector<int> > &reg_name){
    int total_region_size = 0;
    for(vector<int>::const_iterator ii_node = rnode.begin(); ii_node < rnode.end(); ii_node++){
        int node = *ii_node;
        int clust_start = reg_name[node][1];
        int clust_end = reg_name[node][2];
        total_region_size += clust_end - clust_start + 1;
    }
    return total_region_size;
}

// read bam files all together, and merge them
int MergeBams_prep(string *fn, int n, samfile_t **in, heap1_t *heap, uint64_t *idx){
    for(int i = 0; i!=n; ++i){
        heap1_t *h;

        if ((in[i] = samopen(fn[i].c_str(), "rb", 0)) == 0) {
            fprintf(stderr, "[main_samview] fail to open file %s for reading.\n",
                fn[i].c_str());
            continue;
        }
        if (in[i]->header == 0) {
            fprintf(stderr, "[main_samview] fail to read the header (%s).\n",
                fn[i].c_str());
            continue;
        }
        h = heap + i;
        h->i = i;
        h->b = (bam1_t*)calloc(1,sizeof(bam1_t));
        int r;
        if ((r = samread(in[i], h->b)) >= 0) {
            h->pos = ((uint64_t)h->b->core.tid <<32) | (uint32_t)h->b->core.pos << 1 | bam1_strand(h->b);
            h->idx = (*idx)++;
        }
        else {
            h->pos = HEAP_EMPTY;
        }
    }

    ks_heapmake(heap, n, heap);
    return 1;
}

// read config file from bam2cfg.pl
string get_from_line(string line,string search,int flag){
    size_t pos;
    if(flag == 1)
        pos = line.find(search + ":");
    else
        pos = line.find(search);
    if(pos != string::npos){
        if(flag == 1){
            pos = pos + search.length() + 1;
            size_t pos_end = line.find("\t",pos);
            if(pos_end == string::npos)
                pos_end = line.find("\0",pos);
            size_t n = pos_end - pos;
            return line.substr(pos,n);
        }
        else{
            pos = pos + search.length();
            size_t pos_begin = line.find(":", pos);
            if(pos_begin != string::npos){
                string substr_tmp = line.substr(pos,pos_begin-pos);
                if(substr_tmp.find("\t") != string::npos || substr_tmp.find(" ") != string::npos){
                    return search_more(line, search, pos_begin);
                }
                else{
                    size_t pos_end = line.find("\t",pos_begin);
                    if(pos_end == string::npos)
                        pos_end = line.find("\0", pos_begin);
                    size_t n = pos_end - pos_begin - 1;
                    return line.substr(pos_begin+1,n);
                }
            }
            else
                return "NA";
        }
    }
    else
        return "NA";
    return "NA";
}

// augmenting function for reading config file: apply specifically to flag = 0, but the string appeared before, so search the following ones
string search_more(string line, string search,     size_t pos_begin){
    size_t pos = line.find(search, pos_begin);
    string ret;
    if(pos != string::npos){
        size_t pos_begin = line.find(":", pos);
        if(pos_begin != string::npos){
            string substr_tmp = line.substr(pos, pos_begin-pos);
            if(substr_tmp.find("\t") != string::npos || substr_tmp.find(" ") != string::npos){
                ret = search_more(line, search, pos_begin);
                return ret;
            }
            else{
                size_t pos_end = line.find("\t", pos_begin);
                if(pos_end == string::npos)
                    pos_end = line.find("\0",pos_begin);
                size_t n = pos_end - pos_begin - 1;
                return line.substr(pos_begin+1,n);
            }
        }
        else
            return "NA";
    }
    else
        return "NA";
    return "NA";
}

// augmenting fucntion for reading config file
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
                if(pos_begin != string::npos){

                    string substr_tmp = line.substr(pos,pos_begin-pos);
                    if(substr_tmp.find("\t") != string::npos || substr_tmp.find(" ") != string::npos)
                        return "NA";

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

// augmenting function from int to string
string itos(int i){
    stringstream i_str_stream;
    i_str_stream << i;
    return i_str_stream.str();
}

// get the string of the quality of the sequence
string get_string_qual(uint8_t *pt, int32_t length){
    string seq;
    seq.reserve(length);
    for(int i = 0; i < length ; i++){
        seq += char(pt[i] + 33);
    }
    return seq;
}

// get the string of the sequence
string get_string(uint8_t *pt, int32_t length){
    string seq;
    seq.reserve(length);
    for(int i = 0; i < length ; i++)
        seq += bam_nt16_rev_table[bam1_seqi(pt, i)];
    return seq;
}

void write_fastq_for_flag(const string &flag, const vector< vector<string> > &support_reads, const map<string, string> &ReadsOut) {
    map<string,int> pairing;
    for( vector<vector<string> >::const_iterator ii_support_reads = support_reads.begin(); ii_support_reads != support_reads.end(); ii_support_reads ++){
        vector<string> y = *ii_support_reads;
        if(y.size() != 11 || y[5].compare(flag))
            continue;
        //Paradoxically, the first read seen is put in file 2 and the second in file 1
        string fh_tmp_str = (pairing.find(y[0]) != pairing.end()) ? ReadsOut.at(y[8].append("1")) : ReadsOut.at(y[8].append("2"));
        ofstream fh;
        fh.open(fh_tmp_str.c_str(), ofstream::app);
        pairing[y[0]] = 1;
        //Note that no transformation on read bases based on read orientation is done here
        string str_tmp = "@" + y[0] + "\n" + y[9] + "\n" + "+\n" + y[10] + "\n";
        fh << str_tmp;
        fh.close();
    }
}

string choose_sv_flag(const int num_readpairs, const map<string, int> reads_per_type) {
    string flag = "0";
    float max_type_pct = 0.0;
    for(map<string,int>::const_iterator type_iter = reads_per_type.begin(); type_iter != reads_per_type.end(); ++type_iter){
        string current_flag = (*type_iter).first;
        float type_pct = static_cast<float>((*type_iter).second) / static_cast<float>(num_readpairs);
        if(max_type_pct < type_pct) {
            max_type_pct = type_pct;
            flag = current_flag;
        }
    }
    return flag;
}


// The following correspondence of the data structure is between perl and c in samtools. Left column: perl. Right column: cpp
/*
 t.readname         :        bam1_qname(b)  (length: b.core.l_qname)
 t.flag            :        b->core.flag
 t.chr            :        b->core.tid (this transit from char to int)
 t.pos            :        b->core.pos
 t.qual            :        b->core.qual
 mchr            :        b->core.mtid (this transit from char to int)
 mpos            :        b->core.mpos
 t.dist            :        b->core.isize
 t.seq            :        bam1_seq(b)        (length: b.core.l_qseq)
 t.basequal        :        bam1_qual(b)    (length: b.core.l_qseq)
 t.ori            :        ori as a return in AlnParser    (in samtools, can be derived by bam1_strand(b), bam1_mstrand(b))
 t.readgroup        :        readgroup as a pointer in AlnParser input
 t.readlen        :        b->core.l_qseq
 */
