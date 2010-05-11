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
#include "BreakDancerMax.h"
#include "AlnParser.h"
#include "sam.h"
#include "bam.h"

extern string platform;

void AlnParser(bam2_t *b, string format, map readgroup_platform, string alt){
	string platform = platform;
	t_buf t;
	string tmp1, tmp2, tmp3, tmp4, tmp5, tmp6;
	// strcmp = 0 -> two strings are the same
/*	if(! compare(format, "maq")){
 		while (line >> t.readname >> t.chr >> t.pos >> t.ori >> t.dist >> tmp1 >> tmp2 >> t.flag >> tmp3 >> tmp4 >> tmp5 >> tmp6 >> t.readlen >> t.seq >> t.basequal){
 			if(alt != ""){
 				t.qual = tmp1;
 			}
 		}
	}*/ // do not do maq now
	/*else */if(! strcmp(format, "sam")){
//		string flag;
//		string mchr;
//		int mpos;
//		while (line >> t.readname >> flag >> t.chr >> t.pos >> t.qual >> tmp1 >> mchr >> mpos >> t.dist >> t.seq >> t.basequal){
			
		uint32_t flag = b.core.flag;
		b.core.l_qseq = strlen(bam1_seq(b));
		b.ori = (flag&0x0010 || strstr(flag,"r"))?"-":"+";
		b.flag = 0;
		
		if(uint8_t *tmp = bam_aux_get(b, "RG")){
			b.readgroup = bam_aux2Z(tmp);
			platform = readgroup_platform[b.readgroup]?readgroup_platform[b.readgroup]:"illumina";
		}
		
		string string2_cmp = "AQ:i"
		
		if(strlen(alt) != 0){
			if(uint8_t *tmp = bam_aux_get(b, "AQ"))
				b.core.qual = bam_aux2i(tmp);
				else if(uint8_t *tmp = bam_aux_get(b, "AM"))
			b.core.qual = bam_aux2i(tmp);					
		}
		
		// convert to Maq flag
		if(uint8_t *tmp = bam_aux_get(b, "MF"))
			b.flag = bam_aux2i(tmp);
		else{
			if(flag & 0x0400 || strstr(flag,"d"))
				b.flag = 0;
			else if(flag & 0x0001 || strstr(flag,"p")){
				char ori2;
				ori2 = (flag & 0x0020 || strstr(flag,"R"))?"-":"+";
				if(flag & 0x0004 || strstr(flag,"u"))
					b.flag = 192;
				else if(flag & 0x0008 || strstr(flag,,"U"))
					b.flag = 64;
				else if(strcmp(b.core.mtid, "="))
					b.flag = 32;
				else if(flag & 0x0002 || strstr(flag,"P")){
					if(strstr(platform,"solid")) // need to work on insensitive case
						b.flag = 18;
					else{
						if(b.core.pos < b.core.mpos)
							b.flag = (strcmp(b.ori,"+"))?20:18;
						else
							b.flag = (strcmp(b.ori,"+"))?18:20;
					}
				}
				else{
					if(strstr(platform,"solid")) {// insensitive case 
						if(strcmp(b.ori, ori2))
							b.flag = strcmp(ori2, "+")?8:1;
						else if(!strcmp(b.ori, "+")){
							if(flag & 0x0040)
								b.flag = (b.core.pos < b.core.mpos)?2:4;
							else
								b.flag = (b.core.pos > b.core.mpos)?2:4;
						}
						else{
							if(flag & 0x0040)
								b.flag = (b.core.pos > b.core.mpos)?2:4;
							else
								b.flag = (b.core.pos < b.core.mpos)?2:4;
						}
						else
							b.flag = 2;
					}
					else{
						if(!strcmp(b.ori,ori2))
							b.flag = strcmp(ori2, "+")?8:1;
						else if(b.core.mpos > b.core.pos && !strcmp(b.ori, "-") || b.core.pos > b.core.mpos && !strcmp(b.ori, "+"))
							b.flag = 4;
						else
							b.flag = 2;
					}
				}
			}
		}
	}
	return;
}




























