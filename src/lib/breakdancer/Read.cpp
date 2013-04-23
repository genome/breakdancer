#include "Read.hpp"
#include <cstdlib>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace breakdancer;

Read::Read(bam1_t const* record, string const& format, map<string, string> const& readgroup_platform, map<string, string> const& readgroup_library) {
    _record = bam_init1(); 
    bam_copy1(_record, record);
    
    bdflag = _determine_bdflag();
    bdqual = _determine_bdqual();

    platform = _platform(readgroup_platform);
    library = _library(readgroup_library);
    
    _string_record.push_back(queryname());
    _string_record.push_back(boost::lexical_cast<string>(_record->core.tid));
    _string_record.push_back(boost::lexical_cast<string>(_record->core.pos));
    _string_record.push_back(ori());
    _string_record.push_back(boost::lexical_cast<string>(_record->core.isize));
    _string_record.push_back(boost::lexical_cast<string>(bdflag));
    _string_record.push_back(boost::lexical_cast<string>(bdqual));
    _string_record.push_back(boost::lexical_cast<string>(_record->core.l_qseq));
    _string_record.push_back(library);

    if(query_sequence() != "*") {
        _string_record.push_back(query_sequence());
    }
    if(quality_string() != "*") {
        _string_record.push_back(quality_string());
    }
}

Read::~Read() {
    bam_destroy1(_record);
}

string Read::operator[](std::vector<std::string>::size_type idx) {
    return _string_record[idx];
}

string Read::readgroup() {
    if(uint8_t* tmp = bam_aux_get(_record, "RG")) {
        return string(bam_aux2Z(tmp));
    }
    else {
        return "";
    }
}

string Read::queryname() {
    return string(bam1_qname(_record));
}

string Read::query_sequence() {
    uint8_t* seq_ptr = bam1_seq(_record);
    if(_record->core.l_qseq) {
        string seq;
        seq.reserve(_record->core.l_qseq);
        for(int i = 0; i < _record->core.l_qseq; ++i) {
            seq += bam_nt16_rev_table[bam1_seqi(seq_ptr, i)];
        }
        return seq;
    }
    else {
        return "*"; //or maybe throw? I dunno
    }
}

string Read::quality_string() {
    uint8_t* qual_ptr = bam1_qual(_record);
    if(*qual_ptr != 0xff) { //fi 0xff then there is no qual string in BAM
        string qual;
        qual.reserve(_record->core.l_qseq);
        for(int i = 0; i < _record->core.l_qseq; ++i) {
            qual += char(qual_ptr[i] + 33);
        }
        return qual;
    }
    else {
        return "*";
    }
}

string Read::ori() {
    return  _record->core.flag & 0x0010 ? "-" : "+";
}

string Read::_platform(map<string, string> const& readgroup_platform) {
    map<string, string>::const_iterator platform_it = readgroup_platform.find(readgroup());
    if(platform_it != readgroup_platform.end()) {
        return platform_it->second;
    }
    else {
        return "illumina";
    }
}

string Read::_library(map<string, string> const& readgroup_library) {
    //contrary to the regular code, let's assume we've fixed things so that the map returns the filename if the readgroup is empty.
    map<string, string>::const_iterator lib = readgroup_library.find(readgroup());
    if(lib != readgroup_library.end()) {
        return lib->second;
    }
    else {
        return "";
    }
}

int Read::_determine_bdqual() {
    //Breakdancer always takes the alternative mapping quality, if available
    //it originally contained support for AQ, but the newer tag appears to be AM. Dropping support for AQ.
    if(uint8_t* alt_qual = bam_aux_get(_record, "AM")) {
         return bam_aux2i(alt_qual);					
    }
    else {
        return _record->core.qual;  //if no alternative mapping quality, use core quality
    }
}

//FIXME This is practically illegible and there are many flag definitions that need to be added
//In addition, there is some caching that could happen to make this whole thing shorter
//and less computationally expensive
pair_orientation_flag Read::_determine_bdflag() {
    pair_orientation_flag flag = NA;
    if(!(_record->core.flag & BAM_FDUP)) { //this should probably include QCfail as well
        if(_record->core.flag & BAM_FPAIRED) {
            if(_record->core.flag & BAM_FUNMAP) {
                flag = UNMAPPED;
            }
            else if(_record->core.flag & BAM_FMUNMAP) {
                flag = MATE_UNMAPPED;
            }
            else if(_record->core.tid != _record->core.mtid) {
                flag = ARP_CTX;
            }
            else if(_record->core.flag & BAM_FPROPER_PAIR) {
                if(platform == "solid") { //assuming config parser has normalized this for me
                    flag = NORMAL_FR; //normal insert size
                }
                else {
                    if(_record->core.pos < _record->core.mpos) {
                        flag = (_record->core.flag & BAM_FREVERSE) ? NORMAL_RF : NORMAL_FR;
                    }
                    else {
                        flag = (_record->core.flag & BAM_FREVERSE) ? NORMAL_FR : NORMAL_RF;
                    }
                }
            }
            else {
                if(platform == "solid") {
                    if( ((_record->core.flag & BAM_FREVERSE) && !(_record->core.flag & BAM_FMREVERSE)) ||
                        (!(_record->core.flag & BAM_FREVERSE) && (_record->core.flag & BAM_FMREVERSE))) { //do the mates have different orientation?
                        flag = (_record->core.flag & BAM_FREVERSE) ? ARP_RR : ARP_FF;
                    }
                    else if( !(_record->core.flag & BAM_FREVERSE)) {
                        if(_record->core.flag & BAM_FREAD1) {
                            flag = (_record->core.pos < _record->core.mpos) ? ARP_FR_big_insert : ARP_RF;
                        }
                        else {
                            flag = (_record->core.pos > _record->core.mpos) ? ARP_FR_big_insert : ARP_RF;
                        }
                    }
                    else {
                        if(_record->core.flag & BAM_FREAD1) {
                            flag = (_record->core.pos > _record->core.mpos) ? ARP_FR_big_insert : ARP_RF;
                        }
                        else {
                            flag = (_record->core.pos < _record->core.mpos) ? ARP_FR_big_insert : ARP_RF;
                        }
                    }
                }
                else {
                    if( ((_record->core.flag & BAM_FREVERSE) && (_record->core.flag & BAM_FMREVERSE)) ||
                        (!(_record->core.flag & BAM_FREVERSE) && !(_record->core.flag & BAM_FMREVERSE))) { //do the mates have the same orientation?
                    
                        flag = (_record->core.flag & BAM_FMREVERSE) ? ARP_RR : ARP_FF;
                    }
                    else if((_record->core.mpos > _record->core.pos && (_record->core.flag & BAM_FREVERSE)) || (_record->core.pos > _record->core.mpos && !(_record->core.flag & BAM_FREVERSE))) {
                        flag = ARP_RF;
                    }
                    else {
                        flag = ARP_FR_big_insert;
                    }
                }
            }
        }
    }
    return flag;
}
