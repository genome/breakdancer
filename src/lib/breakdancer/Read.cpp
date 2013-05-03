#include "Read.hpp"
#include "LegacyConfig.hpp"

#include <cstdlib>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace breakdancer;

Read::Read(
        bam1_t const* record,
        string const& format,
        LegacyConfig const& cfg
        )
    : _record(bam_init1())
    , _query_name_cached(false)
    , _query_seq_cached(false)
    , _quality_string_cached(false)
    , _tid(record->core.tid)
    , _pos(record->core.pos)
    , _isize(record->core.isize)
    , _abs_isize_cached(false)
    , _query_length(record->core.l_qseq)
    , _ori(record->core.flag & BAM_FREVERSE ? '-' : '+')
    , _bdqual(determine_bdqual(record))
{
    //create copy of bam record
    bam_copy1(_record, record);

    //Be careful below. These must be in this order as _platform and _library depend on readgroup
    readgroup = _readgroup();
    library = _library(cfg.readgroup_library);
    _bdflag = determine_bdflag(record, cfg.platform_for_readgroup(readgroup));
}

Read::Read(const Read& other)
    : _query_name(other._query_name)
    , _query_name_cached(other._query_name_cached)
    , _query_sequence(other._query_sequence)
    , _query_seq_cached(other._query_seq_cached)
    , _quality_string(other._quality_string)
    , _quality_string_cached(other._quality_string_cached)
    , _bdflag(other._bdflag)
    , _ori(other._ori)
    , _bdqual(other._bdqual)
    , _abs_isize(other._abs_isize)
    , _abs_isize_cached(other._abs_isize_cached)
    , readgroup(other.readgroup)
    , library(other.library)
    , platform(other.platform)
{
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
        _ori = other._ori;
        _abs_isize = other._abs_isize;
        _abs_isize_cached = other._abs_isize_cached;
        platform = other.platform;
        library = other.library;
        readgroup = other.readgroup;
        _query_name = other._query_name;
        _query_name_cached = other._query_name_cached;
        _query_sequence = other._query_sequence;
        _query_seq_cached = other._query_seq_cached;
        _quality_string = other._quality_string;
        _quality_string_cached = other._quality_string_cached;
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

string Read::_library(ConfigMap<string, string>::type const& readgroup_library) {
    //contrary to the regular code, let's assume we've fixed things so that the map returns the filename if the readgroup is empty.
    ConfigMap<string, string>::type::const_iterator lib = readgroup_library.find(readgroup);
    if(lib != readgroup_library.end()) {
        return lib->second;
    }
    else {
        return "";
    }
}

//FIXME This is practically illegible and there are many flag definitions that need to be added
//In addition, there is some caching that could happen to make this whole thing shorter
//and less computationally expensive
//ideally the majority of this code could be pulled into another class.
string const& Read::query_name() const {
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
}

pair_orientation_flag const& Read::bdflag() {
    return _bdflag;
}

int const& Read::bdqual() {
    return _bdqual;
}

int const& Read::tid() {
    return _record->core.tid;
}

int const& Read::pos() {
    return _record->core.pos;
}

int const& Read::query_length() {
    return _record->core.l_qseq;
}

char const& Read::ori() {
    return  _ori;
}

int const& Read::isize() {
    return _record->core.isize;
}

int const& Read::abs_isize() {
    if(!_abs_isize_cached) {
        _abs_isize = abs(_record->core.isize);
        _abs_isize_cached = true;
    }
    return _abs_isize;
}
