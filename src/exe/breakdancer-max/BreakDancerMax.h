#include "breakdancer/BDConfig.hpp"
#include "breakdancer/LegacyConfig.hpp"
#include "breakdancer/Read.hpp"
#include "breakdancer/saminternals.h"
#include "breakdancer/BreakDancer.hpp"

#include <boost/format.hpp>

#include <cstddef>
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

/*template <class T>
T from_string(const std::string& s,
                 std::ios_base& (*f)(std::ios_base&))*/
void do_break_func(
    Options const& opts,
    BreakDancer& bdancer,
    LegacyConfig const& cfg,
    vector<breakdancer::Read> const& reg_seq,
    map<string, vector<int> >& read,
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

string itos(int i);

// refactoring functions
// write out reads to fastq files
void write_fastq_for_flag(breakdancer::pair_orientation_flag const& flag, const vector<breakdancer::Read> &support_reads, ConfigMap<string, string>::type const& ReadsOut);

// choose the predominant type of read in a region
breakdancer::pair_orientation_flag choose_sv_flag(const int num_readpairs, const map<breakdancer::pair_orientation_flag, int> reads_per_type);
