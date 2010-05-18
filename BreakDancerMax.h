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
#include "sam.h"
#include "bam.h"

void Analysis (string lib, bam1_t *b, vector<vector<string>> &reg_seq, vector<int,vector<int>> &reg_name, map<string,vector<int>> &read, map<int, vector<vector<string>>> &regs, int &begins, int &beginc, int &lasts, int &lastc, int &idx_buff, int buffer_size, int &nnormal_reads, int min_len, int &normal_switch, int &reg_idx, int transchr_rearrange, int min_map_qual, int Illumina_long_insert, int prefix_fastq);

void buildConnection(map<string,vector<int>> &read, map<int,vector<int>> &reg_name, map<int,vector<vector<string>>> &regs);

int PutativeRegion(vector<int> rnode, map<int,vector<int>> &reg_name);

void EstimatePriorParameters(map<string,string> &fmaps, map<string,string> &readgroup_library, map<string, float> &mean_insertsize, map<string, float> &std_insertsize, map<string,float> &uppercutoff, map<string,float> &lowercutoff, map<string,float> &readlens, int chr);

float mean(vector<int> &stat);

float standard_deviation(vector<int> &stat, float mean);

int PutativeRegion(vector<int> rnode, map<int,vector<int>> &reg_name);

bamFile ReadBamChr_prep(string chr_str, string bam_name, bam_index_t *idx, int tid, int beg, int end, samfile_t *in);

int ReadBamChr(bam1_t *b, bamFile fp, int tid, int beg, int end, bam_index_t *idx, uint64_t curr_off, int i, int n_seeks, pair64_t *off, int n_off);

string get_from_line(string line,string search,int flag);
