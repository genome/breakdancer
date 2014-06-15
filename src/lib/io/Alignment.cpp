#include "Alignment.hpp"

#include <boost/lexical_cast.hpp>

#include <cstddef>
#include <cstdlib>
#include <iostream>

using namespace std;


uint8_t determine_bdqual(bam1_t const* record) {
    // Breakdancer always takes the alternative mapping quality, if available
    // it originally contained support for AQ, but the newer tag appears to be
    // AM. Dropping support for AQ.
    if(uint8_t* alt_qual = bam_aux_get(record, "AM")) {
         return bam_aux2i(alt_qual);
    }
    else {
        // if no alternative mapping quality, use core quality
        return record->core.qual;
    }
}

std::string determine_read_group(bam1_t const* record) {
    if(uint8_t* tmp = bam_aux_get(record, "RG"))
        return bam_aux2Z(tmp);
    return "";
}

Alignment::Alignment()
    : _tid(-1)
    , _pos(-1)
    , _query_length(0)
    , _mtid(-1)
    , _mpos(-1)
    , _abs_isize(-1)
    , _sam_flag(0)
    , _bdqual(0)
    , _lib_index(-1)
    , _bdflag(ReadFlag::NA)
{
}

Alignment::Alignment(bam1_t const* record, bool seq_data)
    : _tid(record->core.tid)
    , _pos(record->core.pos)
    , _query_length(record->core.l_qseq)
    , _mtid(record->core.mtid)
    , _mpos(record->core.mpos)
    , _abs_isize(abs(record->core.isize))
    , _sam_flag(record->core.flag)
    , _bdqual(determine_bdqual(record))
    , _query_name(bam1_qname(record))
    , _lib_index(~0)
    , _bdflag(ReadFlag::NA)
{
    if (seq_data) {
        uint8_t* end = bam1_qual(record);
        if (*end != 0xff)
            end += record->core.l_qseq;
        _bam_data.assign(bam1_seq(record), end);
    }
}

void Alignment::to_fastq(std::ostream& stream) const {
    assert(!_bam_data.empty());
    stream << "@" << query_name() << "\n";
    for (int i = 0; i < _query_length; ++i) {
        stream << char(bam_nt16_rev_table[bam1_seqi(_bam_data.data(), i)]);
    }

    stream << "\n+\n";
    uint8_t const* qdata = _bam_data.data() + ((_query_length+1) >> 1);
    if (qdata[0] != 0xff) {
        assert(qdata + _query_length <= &*_bam_data.end());
        for (int i = 0; i < _query_length; ++i)
            stream << char(qdata[i] + 33);
    }
    else {
        std::cerr << "Warning: no quality data for read " << query_name() << "\n";
    }
    stream << "\n";
}
