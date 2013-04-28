#include "Read.hpp"
#include <cstdlib>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace breakdancer;

Read::Read(bam1_t const* record, string const& format, map<string, string> const& readgroup_platform, map<string, string> const& readgroup_library) :
    _record(bam_init1()), _query_name_cached(false), _query_seq_cached(false), _quality_string_cached(false) {
    //create copy of bam record
    bam_copy1(_record, record);
    
    //set flag and qualities
    _bdflag = _determine_bdflag();
    _bdqual = _determine_bdqual();

    //Be careful below. These must be in this order as _platform and _library depend on readgroup
    readgroup = _readgroup();
    platform = _platform(readgroup_platform);
    library = _library(readgroup_library);
    
    _string_record.push_back("NA");
    _string_record.push_back(boost::lexical_cast<string>(_record->core.tid));
    _string_record.push_back(boost::lexical_cast<string>(_record->core.pos));
    _string_record.push_back(ori());
    _string_record.push_back(boost::lexical_cast<string>(_record->core.isize));
    _string_record.push_back(boost::lexical_cast<string>(_bdflag));
    _string_record.push_back(boost::lexical_cast<string>(_bdqual));
    _string_record.push_back(boost::lexical_cast<string>(_record->core.l_qseq));
    _string_record.push_back(library);
}

Read::Read(const Read& other) : _query_name(other._query_name), _query_name_cached(other._query_name_cached), _query_sequence(other._query_sequence), _query_seq_cached(other._query_seq_cached), _quality_string(other._quality_string), _quality_string_cached(other._quality_string_cached), _bdflag(other._bdflag), _bdqual(other._bdqual), _string_record(other._string_record), readgroup(other.readgroup), platform(other.platform), library(other.library) {
    if(other._record) {
        _record = bam_init1();
        bam_copy1(_record, other._record);
    }
    else {
        _record = NULL;
    }
}

Read::~Read() {
    if(_record != NULL) {
        bam_destroy1(_record);
    }
}

string Read::operator[](std::vector<std::string>::size_type idx) const {
    if(idx == 0) {
        throw "Accessing query name through index now unsupported";
    }
    return _string_record[idx];
}

Read& Read::operator=(const Read& other) {
    if(this != &other) {
        if(_record) {
            bam_destroy1(_record);
        }
        if(other._record != NULL) {
            _record = bam_init1();
            bam_copy1(_record, other._record);
        }
        else {
            _record = NULL;
        }
        _bdflag = other._bdflag;
        _bdqual = other._bdqual;
        platform = other.platform;
        library = other.library;
        readgroup = other.readgroup;
        _query_name = other._query_name;
        _query_name_cached = other._query_name_cached;
        _query_sequence = other._query_sequence;
        _query_seq_cached = other._query_seq_cached;
        _quality_string = other._quality_string;
        _quality_string_cached = other._quality_string_cached;
        _string_record = other._string_record;
    }
    return *this;
}

string Read::_readgroup() {
    if(uint8_t* tmp = bam_aux_get(_record, "RG")) {
        return string(bam_aux2Z(tmp));
    }
    else {
        return "";
    }
}

string Read::_library(map<string, string> const& readgroup_library) {
    //contrary to the regular code, let's assume we've fixed things so that the map returns the filename if the readgroup is empty.
    map<string, string>::const_iterator lib = readgroup_library.find(readgroup);
    if(lib != readgroup_library.end()) {
        return lib->second;
    }
    else {
        return "";
    }
}

string Read::_platform(map<string, string> const& readgroup_platform) {
    map<string, string>::const_iterator platform_it = readgroup_platform.find(readgroup);
    if(platform_it != readgroup_platform.end()) {
        return platform_it->second;
    }
    else {
        return "illumina";
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
//ideally the majority of this code could be pulled into another class.
pair_orientation_flag Read::_determine_bdflag() {
    pair_orientation_flag flag = NA;
    int read_reversed = _record->core.flag & BAM_FREVERSE;
    int mate_reversed = _record->core.flag & BAM_FMREVERSE;
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
                        flag = (read_reversed) ? NORMAL_RF : NORMAL_FR;
                    }
                    else {
                        flag = (read_reversed) ? NORMAL_FR : NORMAL_RF;
                    }
                }
            }
            else {
                if(platform == "solid") {
                    if( ((read_reversed) && !(mate_reversed)) ||
                        (!(read_reversed) && (mate_reversed))) { //do the mates have different orientation?
                        flag = (read_reversed) ? ARP_RR : ARP_FF;
                    }
                    else if( !(read_reversed)) {
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
                    if( ((read_reversed) && (mate_reversed)) ||
                        (!(read_reversed) && !(mate_reversed))) { //do the mates have the same orientation?
                    
                        flag = (mate_reversed) ? ARP_RR : ARP_FF;
                    }
                    else if((_record->core.mpos > _record->core.pos && (read_reversed)) || (_record->core.pos > _record->core.mpos && !(read_reversed))) {
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

vector<string>::size_type Read::size() {
    return _string_record.size();
}

string const& Read::query_name() {
    if(!_query_name_cached) {
        _query_name = string(bam1_qname(_record));
        _query_name_cached = true;
    }
    return _query_name;
}

string const& Read::query_sequence() {
    if(!_query_seq_cached) {
        uint8_t* seq_ptr = bam1_seq(_record);
        if(_record->core.l_qseq) {
            string seq;
            seq.reserve(_record->core.l_qseq);
            for(int i = 0; i < _record->core.l_qseq; ++i) {
                seq += bam_nt16_rev_table[bam1_seqi(seq_ptr, i)];
            }
            _query_sequence = seq;
        }
        else {
            _query_sequence = "*"; //or maybe throw? I dunno
        }
        _query_seq_cached = true;
    }
    return _query_sequence;
}

string const& Read::quality_string() {
    if(!_quality_string_cached) {
    uint8_t* qual_ptr = bam1_qual(_record);
    if(*qual_ptr != 0xff) { //fi 0xff then there is no qual string in BAM
        string qual;
        qual.reserve(_record->core.l_qseq);
        for(int i = 0; i < _record->core.l_qseq; ++i) {
            qual += char(qual_ptr[i] + 33);
        }
        _quality_string = qual;
    }
    else {
        _quality_string = "*";
    }
        _quality_string_cached = true;
    }
    return _quality_string;
}

void Read::set_bdflag(pair_orientation_flag const& new_flag) {
    _bdflag = new_flag;
    _string_record[5] = boost::lexical_cast<string>(new_flag); //this may need a check
}

pair_orientation_flag const& Read::bdflag() {
    return _bdflag;
}

int const& Read::bdqual() {
    return _bdqual;
}

string Read::ori() {
    return  _record->core.flag & 0x0010 ? "-" : "+";
}
