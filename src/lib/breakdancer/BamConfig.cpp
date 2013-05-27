#include "BamConfig.hpp"

#include "BamIO.hpp"
#include "Options.hpp"
#include "Read.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <cstdlib>
#include <memory>

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


BamConfig::BamConfig()
    : max_read_window_size(1e8)
    , max_readlen(0)
    , _covered_ref_len(0)
{
}

BamConfig::BamConfig(std::istream& in, Options const& opts)
    : max_read_window_size(1e8)
    , max_readlen(0)
    , _covered_ref_len(0)
{
    map<string, LibraryInfo> temp_lib_info;
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
        _readgroup_library[readgroup] = lib;

        readgroup_platform[readgroup] = get_from_line(line,"platform",1);

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

        LibraryInfo lib_info;
        lib_info.name = lib;
        lib_info.bam_file = fmap;

        if(mqual_.compare("NA")){
            mqual = atoi(mqual_.c_str());
            lib_info.min_mapping_quality = mqual;
        } else {
            lib_info.min_mapping_quality = -1;
        }
        fmaps[fmap] = lib;

        if(mean_ != "NA" && std_ != "NA") {
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

        lib_info.mean_insertsize = mean;
        lib_info.std_insertsize = std;
        lib_info.uppercutoff = upper;
        lib_info.lowercutoff = lower;
        lib_info.readlens = readlen;

        // This is silly, library info can be repeated in the config file
        // hopefully whatever we are clobbering is equivalent!
        temp_lib_info[lib] = lib_info;

        if(exes.find(fmap) == exes.end())
            exes[fmap] = exe.compare("NA")?exe:"cat";
        else if(exes[fmap].compare(exe) != 0){
            throw runtime_error(
                "Please use identical exe commands to open the same input file.\n"
                );
        }

        int tmp = mean - readlen*2;    // this determines the mean of the max of the SV flanking region
        max_read_window_size = std::min(max_read_window_size, tmp);
    }


    for (map<string, LibraryInfo>::iterator i = temp_lib_info.begin(); i != temp_lib_info.end(); ++i) {
        i->second.index = _library_info.size();
        _lib_names_to_indices[i->first] = i->second.index;
        _library_info.push_back(i->second);
    }


    typedef ConfigMap<string, string>::type::const_iterator IterType;
    for (IterType iter = fmaps.begin(); iter != fmaps.end(); ++iter) {
        _bam_files.push_back(iter->first);
        auto_ptr<IBamReader> reader(openBam(iter->first, opts));
        _analyze_bam(*reader, opts);
    }

    max_read_window_size = std::max(max_read_window_size, 50);
}

void BamConfig::_analyze_bam(IBamReader& reader, Options const& opts) {
    int last_pos = 0;
    int last_tid = -1;

    size_t ref_len = 0;
    uint32_t read_count = 0;

    bam1_t* b = bam_init1();
    while (reader.next(b) > 0) {
        breakdancer::Read aln(b, false);

        string const& lib = readgroup_library(aln.readgroup());
        if (lib.empty())
            continue;

        size_t lib_idx = _lib_names_to_indices[lib];
        LibraryInfo& lib_info = _library_info[lib_idx];

        if (last_tid >= 0 && last_tid == aln.tid())
            ref_len += aln.pos() - last_pos;

        last_pos = aln.pos();
        last_tid = aln.tid();

        if(aln.bdqual() > opts.min_map_qual && (aln.bdflag() == breakdancer::NORMAL_FR || aln.bdflag() == breakdancer::NORMAL_RF)) {
            ++lib_info.read_count; // FIXME: this is a bad side effect modifying another object
            ++read_count;
        }

        //XXX This seems weird to me as well. Why are we checking the
        //mapping quality in two different places and performing some
        //calculations before this is applied?
        //-dlarson

        int min_mapq = lib_info.min_mapping_quality < 0 ?
                opts.min_map_qual : lib_info.min_mapping_quality;

        if (aln.bdqual() <= min_mapq)
            continue;

        if(aln.bdflag() == breakdancer::NA)
            continue;

        if((opts.transchr_rearrange && aln.bdflag() != breakdancer::ARP_CTX) || aln.bdflag() == breakdancer::MATE_UNMAPPED || aln.bdflag() == breakdancer::UNMAPPED)
            continue;

        //It would be nice if this was pulled into the Read class as well
        //for now, let's just set the bdflag directly here since it is public
        if(opts.Illumina_long_insert) {
            if(aln.abs_isize() > lib_info.uppercutoff && aln.bdflag() == breakdancer::NORMAL_RF) {
                aln.set_bdflag(breakdancer::ARP_RF);
            }
            if(aln.abs_isize() < lib_info.uppercutoff && aln.bdflag() == breakdancer::ARP_RF) {
                aln.set_bdflag(breakdancer::NORMAL_RF);
            }
            if(aln.abs_isize() < lib_info.lowercutoff && aln.bdflag() == breakdancer::NORMAL_RF) {
                aln.set_bdflag(breakdancer::ARP_FR_small_insert); //FIXME this name doesn't make a whole lot of sense here
            }
        }
        else{
            if(aln.abs_isize() > lib_info.uppercutoff && aln.bdflag() == breakdancer::NORMAL_FR) {
                aln.set_bdflag(breakdancer::ARP_FR_big_insert);
            }
            if(aln.abs_isize() < lib_info.uppercutoff && aln.bdflag() == breakdancer::ARP_FR_big_insert) {
                aln.set_bdflag(breakdancer::NORMAL_FR);
            }
            if(aln.abs_isize() < lib_info.lowercutoff && aln.bdflag() == breakdancer::NORMAL_FR) {
                aln.set_bdflag(breakdancer::ARP_FR_small_insert);
            }
        }

        if(aln.bdflag() == breakdancer::NORMAL_FR || aln.bdflag() == breakdancer::NORMAL_RF) {
            continue;
        }

        ++lib_info.read_counts_by_flag[aln.bdflag()];
    }
    bam_destroy1(b);

    if(ref_len == 0) {
        cerr << "Input file " << reader.path() <<
            " does not contain legitimate paired end alignment. "
            "Please check that you have the correct paths and the "
            "map/bam files are properly formated and indexed.";
    }

    _read_count_per_bam[reader.path()] = read_count;

    if (_covered_ref_len < ref_len)
        _covered_ref_len = ref_len;
}
