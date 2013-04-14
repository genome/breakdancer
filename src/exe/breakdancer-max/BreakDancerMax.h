#include "breakdancer/AlnParser.h"
#include "breakdancer/samtools.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <list>
#include <iomanip>
#include <cmath>
#include <math.h>
#include <time.h>
#include <map>
#include <assert.h>
#include <sstream>
#include "sam.h"
#include "bam.h"
#include "ksort.h"
#include "khash.h"

using namespace std;

#define LZERO -99
#define ZERO exp(LZERO)

struct Options {
    Options()
        : chr("0")
        , min_len(7)
        , cut_sd(3)
        , max_sd(1000000000)
        , min_map_qual(35)
        , min_read_pair(2)
        , seq_coverage_lim(1000)
        , buffer_size(100)
        , learn_par(false)
        , prior_prob(0.001)
        , transchr_rearrange(false)
        , fisher(false)
        , Illumina_long_insert(false)
        , Illumina_to_SOLiD(false)
        , CN_lib(false)
        , print_AF(false)
        , score_threshold(30)
    {
    }

    string chr;
    int min_len;
    int cut_sd;
    int max_sd;
    int min_map_qual;
    int min_read_pair;
    int seq_coverage_lim;
    int buffer_size;
    bool learn_par;
    float prior_prob;
    bool transchr_rearrange;
    bool fisher;
    bool Illumina_long_insert;
    bool Illumina_to_SOLiD;
    bool CN_lib;
    bool print_AF;
    int score_threshold;
    string bam_file;
    string prefix_fastq;
    string dump_BED;
    string platform;
};

struct BreakDancerData {
    BreakDancerData()
        : begins(-1)
        , beginc(-1)
        , lasts(-1)
        , lastc(-1)
    {
    }

    int begins; // global (chr)
    int beginc; // global
    int lasts; // global (chr, should be int in samtools)
    int lastc; // global
    map<string, uint32_t> nread_ROI; // global
    map<int, map<string, uint32_t> > read_count_ROI_map; // global
    map<string, uint32_t> nread_FR;    // global
    map<int, map<string, uint32_t> > read_count_FR_map; // global
};

struct AnalysisData {
    map<string, uint32_t > possible_fake_data;
    map<int, vector<vector<string> > > regs;//global in analysis
    map<string, vector<int> > read;// global in analysis
    map<int,vector<int> > reg_name;// global in analysis
    vector<vector<string> > reg_seq; // global need to see if it's the key or value of one of the above global. should be a string
    vector<vector<string> >::const_iterator it_reg_seq; // global
};
    

/*template <class T>
T from_string(const std::string& s,
                 std::ios_base& (*f)(std::ios_base&))*/
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
    samfile_t *in,
    int seq_coverage_lim,
    uint32_t *ntotal_nucleotides,
    map<string, float> &read_density,
    int CN_lib,
    map<string, string> libmaps,
    vector<string> maps,
    int print_AF,
    int score_threshold
    );


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
    samfile_t *in,
    uint32_t *ntotal_nucleotides,
    map<string, float> &read_density,
    map<string, uint32_t> &possible_fake_data,
    map<string, string> libmaps,
    vector<string> maps
    );

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
    samfile_t *in,
    map<string, float> &read_density,
    int CN_lib,
    vector<string> maps, // FIXME: should be constref
    int print_AF,
    int score_threshold,
    map<string, string> libmaps // FIXME: should be constref
);

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
);

float mean(vector<int> &stat);

float standard_deviation(vector<int> &stat, float mean);

int PutativeRegion(vector<int> &rnode, map<int,vector<int> > &reg_name);

pair64_t * ReadBamChr_prep(string chr_str, string bam_name, int *tid, int *beg, int *end, samfile_t *in, int *n_off);

int ReadBamChr(bam1_t *b, bamFile fp, int tid, int beg, int end, uint64_t *curr_off, int *i, int *n_seeks, pair64_t *off, int n_off);

int MergeBams_prep(string *fn, int n, samfile_t **in, heap1_t *heap, uint64_t *idx);

int MergeBamsChr_prep(string *fn, int n, bamFile *fp, heap1_t *heap, string chr_str, int *tid, int *beg, int *end, samfile_t **in, pair64_t **off, int *n_off, uint64_t *idx);

string get_from_line(string line,string search,int flag);

string get_from_line_two(string line,string search1,string search2,int flag2);

/*string itoa(int i);

int atoi(string str);

float atof(string str);*/

string itos(int i);

string get_string_qual(uint8_t *pt, int32_t length);

string get_string(uint8_t *pt, int32_t length);

string get_string_qual(uint8_t *pt, int32_t length);

string search_more(string line, string search,     size_t pos_begin);

string char2str(char *str_);

// refactoring functions
void write_fastq_for_flag(const string &flag, const vector< vector<string> > &support_reads, const map<string, string> &ReadsOut);


