#include "AlnParser.h"

using namespace std;
extern string platform;
//extern std::map<char *, std::string> readgroup_platform;

char AlnParser(bam1_t *b, string format, string alt, char *readgroup, map<string, string> &readgroup_platform){
	char ori;
	string platform = platform;
	// strcmp = 0 -> two strings are the same
/*	if(! compare(format, "maq")){
 		while (line >> t.readname >> t.chr >> t.pos >> t.ori >> t.dist >> tmp1 >> tmp2 >> t.flag >> tmp3 >> tmp4 >> tmp5 >> tmp6 >> t.readlen >> t.seq >> t.basequal){
 			if(alt != ""){
 				t.qual = tmp1;
 			}
 		}
	}*/ // do not do maq now
	/*else */if(! format.compare("sam")){
//		string flag;
//		string mchr;
//		int mpos;
//		while (line >> t.readname >> flag >> t.chr >> t.pos >> t.qual >> tmp1 >> mchr >> mpos >> t.dist >> t.seq >> t.basequal){
			
		uint32_t flag = b->core.flag;
		// convert to string
		std::string str_flag;
		std::stringstream ss_flag;
		ss_flag << flag;
		ss_flag >> str_flag;
		// convert b->core.mtid to string
		std::string str_mtid;
		std::stringstream ss_mtid;
		ss_mtid << b->core.mtid;
		ss_mtid >> str_mtid;
		
		ori = (flag&0x0010 || str_flag.find("r") != string::npos)?'-':'+';
		b->core.flag = 0;
		
		if(uint8_t *tmp = bam_aux_get(b, "RG")){
			readgroup = bam_aux2Z(tmp);
			//platform = readgroup_platform[readgroup];
			platform = (readgroup_platform.count(readgroup)>0)?readgroup_platform[readgroup]:"illumina";
		}
		
		string string2_cmp = "AQ:i";
		
		if(alt.length() != 0){
			if(uint8_t *tmp = bam_aux_get(b, "AQ"))
				b->core.qual = bam_aux2i(tmp);
				else if(uint8_t *tmp = bam_aux_get(b, "AM"))
			b->core.qual = bam_aux2i(tmp);					
		}
		
		// convert to Maq flag
		if(uint8_t *tmp = bam_aux_get(b, "MF"))
			b->core.flag = bam_aux2i(tmp);
		else{
			if(flag & 0x0400 || str_flag.find("d")!=string::npos)
				b->core.flag = 0;
			else if(flag & 0x0001 || str_flag.find("p")){
				char ori2;
				ori2 = (flag & 0x0020 || str_flag.find("R"))?'-':'+';
				if(flag & 0x0004 || str_flag.find("u"))
					b->core.flag = 192;
				else if(flag & 0x0008 || str_flag.find("U"))
					b->core.flag = 64;
				else if(str_mtid.find("=")!=string::npos)
					b->core.flag = 32;
				else if(flag & 0x0002 || str_flag.find("P")){
					if(platform.compare("solid")) // need to work on insensitive case
						b->core.flag = 18;
					else{
						if(b->core.pos < b->core.mpos)
							b->core.flag = (ori != '+')?20:18;
						else
							b->core.flag = (ori != '+')?18:20;
					}
				}
				else{
					if(platform.find("solid")!=string::npos) {// insensitive case 
						if(ori != ori2)
							b->core.flag = (ori2 != '+')?8:1;
						else if(ori == '+'){
							if(flag & 0x0040)
								b->core.flag = (b->core.pos < b->core.mpos)?2:4;
							else
								b->core.flag = (b->core.pos > b->core.mpos)?2:4;
						}
						else{
							if(flag & 0x0040)
								b->core.flag = (b->core.pos > b->core.mpos)?2:4;
							else
								b->core.flag = (b->core.pos < b->core.mpos)?2:4;
						}
						//else
							//b->core.flag = 2;
					}
					else{
						if(ori != ori2)
							b->core.flag = (ori2 != '+')?8:1;
						else if(b->core.mpos > b->core.pos && (ori == '-') || b->core.pos > b->core.mpos && (ori == '+'))
							b->core.flag = 4;
						else
							b->core.flag = 2;
					}
				}
			}
		}
	}
	return ori;
}




























