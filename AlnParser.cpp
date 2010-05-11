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

extern string platform;

t_buf *in(string line, string format, readgroup_platform_buf *readgroup_platform, string alt){
	string platform = platform;
	t_buf t;
	string tmp1, tmp2, tmp3, tmp4, tmp5, tmp6;
	// strcmp = 0 -> two strings are the same
	if(! compare(format, "maq")){
 		while (line >> t.readname >> t.chr >> t.pos >> t.ori >> t.dist >> tmp1 >> tmp2 >> t.flag >> tmp3 >> tmp4 >> tmp5 >> tmp6 >> t.readlen >> t.seq >> t.basequal){
 			if(alt != ""){
 				t.qual = tmp1;
 			}
 		}
	}
	else if(! compare(format, "sam")){
		string flag;
		string mchr;
		int mpos;
		while (line >> t.readname >> flag >> t.chr >> t.pos >> t.qual >> tmp1 >> mchr >> mpos >> t.dist >> t.seq >> t.basequal){
			t.readlen = t.seq.length();
			t.ori = (t.ori = flag&0x0010 || flag.contains("r"))?"-":"+"s
			t.flag = 0;
			string string1_cmp = "RG:Z:";
			string string1_tmp = string1_cmp + RXalpha;
			if(line.contains(string1_tmp)){
				t.readgroup = // how to extract the number? $1
				if(readgroup_platform.(t.readgroup))
					platform = readgroup_platform.(t.readgroup);
				else
					platform = "illumina";
			}
			string string2_cmp = "AQ:i"
			
			
			if(alt != ""){
				string string2_cmp = "AQ:i:";
				string string2_tmp = string2_cmp + RXint;
				string string3_cmp = "AM:i:";
				string string3_tmp = string3_cmp + RXint;
				if(line.contains(string2_tmp))
					t.qual = // same problem
				else if(line.contains(string3_tmp))
					t.qual = //same problem
			}
			
			// convert to Maq flag
			string string4_tmp = "MF:i:" + RXint;
			if(line.contains(string4_tmp))
				t.flag = //
			else{
				if(flag & 0x0400 || flag.contains("d"))
					t.flag = 0;
				else if(flag & 0x0001 || flag.contains("p")){
					char ori2;
					ori2 = (flag & 0x0020 || flag.contains("R"))?"-":"+";
					if(flag & 0x0004 || flag.contains("u"))
						t.flag = 192;
					else if(flag & 0x0008 || flag.contains("U"))
						t.flag = 64;
					else if(compare(mchr, "="))
						t.flag = 32;
					else if(flag & 0x0002 || flag.contains("P")){
						if(platform.contains("solid")) // need to work on insensitive case
							t.flag = 18;
						else{
							if(t.pos < mpos)
								t.flag = (compare(t.ori,"+"))?20:18;
							else
								t.flag = (compare(t.ori,"+"))?18:20;
						}
					}
					else{
						if(platform.contains("solid")) {// insensitive case 
							if(compare(t.ori, ori2))
								t.flag = compare(ori2, "+")?8:1;
							else if(!compare(t.ori, "+")){
								if(flag & 0x0040)
									t.flag = (t.pos < mpos)?2:4;
								else
									t.flag = (t.pos > mpos)?2:4;
							}
							else{
								if(flag & 0x0040)
									t.flag = (t.pos > mpos)?2:4;
								else
									t.flag = (t.pos < mpos)?2:4;
							}
							else
								t.flag = 2;
						}
						else{
							if(!compare(t.ori,ori2))
								t.flag = compare(ori2, "+")?8:1;
							else if(mpos > t.pos && !compare(t.ori, "-") || t.pos > mpos && !compare(t.ori, "+"))
								t.flag = 4;
							else
								t.flag = 2;
						}
					}
				}
			}
		}
	}
	else{}
	
	return t;
	
}




























