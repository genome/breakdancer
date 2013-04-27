#include "breakdancer/BDConfig.hpp"
#include "breakdancer/LegacyConfig.hpp"
#include "breakdancer/Read.hpp"
#include "breakdancer/saminternals.h"

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

class Options;

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
    map<int, vector<breakdancer::Read> > regs;//global in analysis
    map<string, vector<int> > read;// global in analysis
    map<int,vector<int> > reg_name;// global in analysis
    vector<breakdancer::Read> reg_seq; // global need to see if it's the key or value of one of the above global. should be a string
    vector<breakdancer::Read>::const_iterator it_reg_seq; // global
};

/*template <class T>
T from_string(const std::string& s,
                 std::ios_base& (*f)(std::ios_base&))*/
void do_break_func(
    Options const& opts,
    BreakDancerData& bdancer,
    LegacyConfig const& cfg,
    vector<breakdancer::Read> const& reg_seq,
    map<int, vector<int> >& reg_name,
    map<string, vector<int> >& read,
    map<int, vector<breakdancer::Read> > &regs,
    int *idx_buff,
    int *nnormal_reads,
    int *reg_idx,
    map<uint32_t, map<string,int> > &x_readcounts,
    uint32_t reference_len,
    map<breakdancer::pair_orientation_flag, string> &SVtype,
    int *max_readlen,
    bam_header_t* bam_header,
    uint32_t *ntotal_nucleotides,
    map<string, float>& read_density,
    vector<string> maps
    );

float mean(vector<int> &stat);

float standard_deviation(vector<int> &stat, float mean);

int PutativeRegion(vector<int> &rnode, map<int,vector<int> > &reg_name);

/*string itoa(int i);

int atoi(string str);

float atof(string str);*/

string itos(int i);

string get_string_qual(uint8_t *pt, int32_t length);

string get_string(uint8_t *pt, int32_t length);

string get_string_qual(uint8_t *pt, int32_t length);

// refactoring functions
// write out reads to fastq files
void write_fastq_for_flag(breakdancer::pair_orientation_flag const& flag, const vector<breakdancer::Read> &support_reads, ConfigMap<string, string>::type const& ReadsOut);

// choose the predominant type of read in a region
breakdancer::pair_orientation_flag choose_sv_flag(const int num_readpairs, const map<breakdancer::pair_orientation_flag, int> reads_per_type);
