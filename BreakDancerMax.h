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

void buildConnection(map<string,vector<int>> &read, map<int,vector<int>> &reg_name, map<int,vector<vector<string>>> &regs);

int PutativeRegion(vector<int> rnode, map<int,vector<int>> &reg_name);

void EstimatePriorParameters(map<string,string> &fmaps, map<string,string> &readgroup_library, map<string, float> &mean_insertsize, map<string, float> &std_insertsize, map<string,float> &uppercutoff, map<string,float> &lowercutoff, map<string,float> &readlens, int chr);

float mean(vector<int> &stat);

float standard_deviation(vector<int> &stat, float mean);

int PutativeRegion(vector<int> rnode, map<int,vector<int>> &reg_name);

bamFile ReadBamChr_prep(string chr_str, string bam_name, bam_index_t *idx, int tid, int beg, int end, samfile_t *in);

int ReadBamChr(bam1_t *b, bamFile fp, int tid, int beg, int end, bam_index_t *idx, uint64_t curr_off, int i, int n_seeks, pair64_t *off, int n_off);

string get_from_line(string line,string search,int flag);
