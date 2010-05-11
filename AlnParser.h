#ifndef ALNPARSER_H
#define ALNPARSER_H

struct t_buf {
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

struct readgroup_platform_buf {
	int readgroup;
	string platform;	
};