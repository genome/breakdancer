#include "Read.hpp"

#include <cstdlib>
#include <cstddef>
#include <boost/lexical_cast.hpp>

using namespace std;


int determine_bdqual(bam1_t const* record) {
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

//FIXME This is practically illegible and there are many flag definitions that need to be added
//In addition, there is some caching that could happen to make this whole thing shorter
//and less computationally expensive
//ideally the majority of this code could be pulled into another class.
ReadFlag determine_bdflag(bam1_t const* record) {

    // this should probably include QCfail as well
    if((record->core.flag & BAM_FDUP) || !(record->core.flag & BAM_FPAIRED))
        return NA;

    bool read_reversed = record->core.flag & BAM_FREVERSE;
    bool mate_reversed = record->core.flag & BAM_FMREVERSE;
    ReadFlag flag = NA;

    if(record->core.flag & BAM_FUNMAP) {
        flag = UNMAPPED;
    }
    else if(record->core.flag & BAM_FMUNMAP) {
        flag = MATE_UNMAPPED;
    }
    else if(record->core.tid != record->core.mtid) {
        flag = ARP_CTX;
    }
    else if(record->core.flag & BAM_FPROPER_PAIR) {
        if(record->core.pos < record->core.mpos) {
            flag = (read_reversed) ? NORMAL_RF : NORMAL_FR;
        }
        else {
            flag = (read_reversed) ? NORMAL_FR : NORMAL_RF;
        }
    }
    else if( ((read_reversed) && (mate_reversed)) ||
            (!(read_reversed) && !(mate_reversed))) { //do the mates have the same orientation?

            flag = (mate_reversed) ? ARP_RR : ARP_FF;
    }
    else if((record->core.mpos > record->core.pos && (read_reversed)) || (record->core.pos > record->core.mpos && !(read_reversed))) {
        flag = ARP_RF;
    }
    else {
        flag = ARP_LARGE_INSERT;
    }

    return flag;
}

Read::Read(bam1_t const* record, bool seq_data)
    : _sam_flag(record->core.flag)
    , _abs_isize(abs(record->core.isize))
    , _bdqual(determine_bdqual(record))
    , _pos(record->core.pos)
    , _mpos(record->core.mpos)
    , _query_length(record->core.l_qseq)
    , _tid(record->core.tid)
    , _mtid(record->core.mtid)
    , _query_name(bam1_qname(record))
    , _lib_index(~0)
    , _bdflag(determine_bdflag(record))
{
    if (seq_data) {
        uint8_t* end = bam1_qual(record);
        if (*end != 0xff)
            end += record->core.l_qseq;
        _bam_data.assign(bam1_seq(record), end);
    }

    if(uint8_t* tmp = bam_aux_get(record, "RG"))
        _readgroup = bam_aux2Z(tmp);
}

void Read::to_fastq(std::ostream& stream) const {
    stream << "@" << query_name() << "\n";
    for (int i = 0; i < _query_length; ++i) {
        stream << char(bam_nt16_rev_table[bam1_seqi(_bam_data.data(), i)]);
    }

    stream << "\n+\n";
    uint8_t const* qdata = _bam_data.data() + ((_query_length+1) >> 1);
    if (qdata[0] != 0xff) {
        for (int i = 0; i < _query_length; ++i)
            stream << char(qdata[i] + 33);
    }
    stream << "\n";
}
