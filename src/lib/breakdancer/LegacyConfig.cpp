#include "LegacyConfig.hpp"
#include "Options.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <cstdlib>

using namespace std;

namespace {
    // augmenting function for reading config file: apply specifically to flag = 0, but the string appeared before, so search the following ones
    string search_more(string line, string search, size_t pos_begin){
        size_t pos = line.find(search, pos_begin);
        string ret;
        if(pos != string::npos){
            size_t pos_begin = line.find(":", pos);
            if(pos_begin != string::npos){
                string substr_tmp = line.substr(pos, pos_begin-pos);
                if(substr_tmp.find("\t") != string::npos || substr_tmp.find(" ") != string::npos){
                    ret = search_more(line, search, pos_begin);
                    return ret;
                }
                else{
                    size_t pos_end = line.find("\t", pos_begin);
                    if(pos_end == string::npos)
                        pos_end = line.find("\0",pos_begin);
                    size_t n = pos_end - pos_begin - 1;
                    return line.substr(pos_begin+1,n);
                }
            }
            else
                return "NA";
        }
        else
            return "NA";
        return "NA";
    }


    string get_from_line(string line,string search,int flag){
        size_t pos;
        if(flag == 1)
            pos = line.find(search + ":");
        else
            pos = line.find(search);
        if(pos != string::npos){
            if(flag == 1){
                pos = pos + search.length() + 1;
                size_t pos_end = line.find("\t",pos);
                if(pos_end == string::npos)
                    pos_end = line.find("\0",pos);
                size_t n = pos_end - pos;
                return line.substr(pos,n);
            }
            else{
                pos = pos + search.length();
                size_t pos_begin = line.find(":", pos);
                if(pos_begin != string::npos){
                    string substr_tmp = line.substr(pos,pos_begin-pos);
                    if(substr_tmp.find("\t") != string::npos || substr_tmp.find(" ") != string::npos){
                        return search_more(line, search, pos_begin);
                    }
                    else{
                        size_t pos_end = line.find("\t",pos_begin);
                        if(pos_end == string::npos)
                            pos_end = line.find("\0", pos_begin);
                        size_t n = pos_end - pos_begin - 1;
                        return line.substr(pos_begin+1,n);
                    }
                }
                else
                    return "NA";
            }
        }
        else
            return "NA";
        return "NA";
    }

    // augmenting fucntion for reading config file
    string get_from_line_two(string line,string search1,string search2,int flag2){
        size_t pos = line.find(search1);
        size_t pos2;
        if(pos != string::npos){
            pos = pos + search1.length();
            pos2 = line.find(search2,pos);
            if(pos2 != string::npos){
                if(flag2 == 1){
                    pos2 = pos2 + search2.length();
                    if(line.substr(pos2,1).compare(":") == 0){
                        size_t pos_end = line.find("\t",pos2);
                        if(pos_end == string::npos)
                            pos_end = line.find("\0",pos2);
                        size_t n = pos_end - pos2 - 1;
                        return line.substr(pos2+1, n);
                    }
                    else
                        return "NA";
                }
                else{
                    pos2 = pos2 + search2.length();
                    size_t pos_begin = line.find(":", pos2);
                    if(pos_begin != string::npos){

                        string substr_tmp = line.substr(pos,pos_begin-pos);
                        if(substr_tmp.find("\t") != string::npos || substr_tmp.find(" ") != string::npos)
                            return "NA";

                        size_t pos_end = line.find("\t",pos2);
                        if(pos_end == string::npos)
                            pos_end = line.find("\0",pos2);
                        size_t n = pos_end - pos_begin - 1;
                        return line.substr(pos_begin, n);
                    }
                    else
                        return "NA";
                }
            }
            else
                return "NA";
        }
        else
            return "NA";
        return "NA";
    }
}


LegacyConfig::LegacyConfig(std::istream& in, Options const& opts)
    : d(1e8)
    , max_readlen(0)
{
    string line;
    while(getline(in, line)) {
            if(line.empty())
                break;

            // analyze the line
            string fmap = get_from_line(line,"map",1);
            string mean_ = get_from_line(line,"mean",0);
            string std_ = get_from_line(line,"std",0);
            string readlen_ = get_from_line(line,"readlen",0);
            string upper_ = get_from_line(line,"upp",0);
            string lower_ = get_from_line(line,"low",0);
            string mqual_ = get_from_line_two(line,"map","qual",0);
            string lib = get_from_line(line,"lib",0);
            float mean=0.0f,std=0.0f,readlen=0.0f,upper=0.0f,lower=0.0f;
            int mqual;
            if(lib.compare("NA")==0)
                lib = get_from_line(line,"samp",0);


            string readgroup = get_from_line(line,"group",1);
            if(readgroup.compare("NA")==0)
                readgroup = lib;
            readgroup_library[readgroup] = lib;

            string platform = get_from_line(line,"platform",1);
            if(opts.Illumina_to_SOLiD)
                readgroup_platform[readgroup] = "solid";
            else
                readgroup_platform[readgroup] = "illumina";
            readgroup_platform[readgroup] = platform;

            string exe = get_from_line(line,"exe",0);
            if(!opts.prefix_fastq.empty()) {
                //ofstream ReadsOut[lib.append("1")](prefix_fastq.append(lib).append(".1.fastq"), ios::out | ios::app | ios::binary);
                ReadsOut[lib + "1"] = opts.prefix_fastq + "." + lib + ".1.fastq";
                ofstream ReadsOutTmp;
                ReadsOutTmp.open(ReadsOut[lib+"1"].c_str());
                if(!ReadsOutTmp.is_open())
                    cout << "unable to open " << opts.prefix_fastq << "." << lib << ".1.fastq, check write permission\n";
                ReadsOutTmp.close();
                ReadsOut[lib + "2"] = opts.prefix_fastq + "." + lib + ".2.fastq";
                ReadsOutTmp.open(ReadsOut[lib+"2"].c_str());
                if(!ReadsOutTmp.is_open())
                    cout << "unable to open " << opts.prefix_fastq << "." << lib << ".2.fastq, check write permission\n";
                ReadsOutTmp.close();
            }

            libmaps[lib] = fmap;
            if(mqual_.compare("NA")){
                mqual = atoi(mqual_.c_str());
                mapQual[lib] = mqual;
            }
            fmaps[fmap] = lib;

            if(mean_.compare("NA") && std_.compare("NA")){
                mean = atof(mean_.c_str());
                std = atof(std_.c_str());
                if(!upper_.compare("NA") || !lower_.compare("NA")){
                    upper = mean + std*float(opts.cut_sd);
                    lower = mean - std*float(opts.cut_sd);
                    lower = lower > 0 ? lower : 0;//atof(lower_.c_str()) > 0 ? atof(lower_.c_str()):0; this is not related with lower_
                }
                else{
                    upper = atof(upper_.c_str());
                    lower = atof(lower_.c_str());
                }
            }

            if(readlen_.compare("NA")){
                readlen = atof(readlen_.c_str());
            }
            max_readlen = max_readlen < readlen ? readlen:max_readlen;

            mean_insertsize[lib] = mean;
            std_insertsize[lib] = std;
            uppercutoff[lib] = upper;
            lowercutoff[lib] = lower;
            readlens[lib] = readlen;
            if(exes.find(fmap) == exes.end())
                exes[fmap] = exe.compare("NA")?exe:"cat";
            else if(exes[fmap].compare(exe) != 0){
                throw runtime_error(
                    "Please use identical exe commands to open the same input file.\n"
                    );
            }

            int tmp = mean - readlen*2;    // this determines the mean of the max of the SV flanking region
            d = d<tmp ? d:tmp;
        }

    d = std::max(d, 50);
}
