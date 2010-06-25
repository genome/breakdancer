#ifndef ALNPARSER_H
#define ALNPARSER_H
#endif

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
#include <sstream>
#include "sam.h"
#include "bam.h"
#include <map>

using namespace std;

vector<string> AlnParser(bam1_t *b, string format, string alt, map<string, string> &readgroup_platform, int same_tid, string platform);

/*struct t_buf {
	string readname;
	string chr;
	int pos;
	int qual;
	int dist;
	string seq;
	string basequal;
	char ori;
	int readlen;
	int flag;
};

// no need for the struct for a map(hash)
struct readgroup_platform_buf {
	int readgroup;
	string platform;	
};*/

/*typedef struct {
	bam1_core_t core;
	int l_aux, data_len, m_data;
	uint8_t *data;
//	int readlen; -> bam2_t.core.l_qseq
	char ori;
//	int flag; -> bam2_t.core.flag
	char *readgroup;
//	int qual; -> bam2_t.core.qual
//	int pos; -> bam2_t.core.pos
} bam2_t;

#define bam_init2() ((bam2_t*)calloc(1, sizeof(bam2_t)))*/
