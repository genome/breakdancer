#include "Read.hpp"

#include <cstdlib>
#include <cstddef>
#include <boost/lexical_cast.hpp>

using namespace std;

namespace {
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
        flag = ARP_FR_big_insert;
    }

    return flag;
}

Read::Read(bam1_t const* record, bool seq_data)
    : _sam_flag(record->core.flag)
    , _bdflag(determine_bdflag(record))
    , _ori(record->core.flag & BAM_FREVERSE ? REV : FWD)
    , _abs_isize(abs(record->core.isize))
    , _bdqual(determine_bdqual(record))
    , _isize(record->core.isize)
    , _pos(record->core.pos)
    , _mpos(record->core.mpos)
    , _endPos(bam_calend(&record->core, bam1_cigar(record)))
    , _query_length(record->core.l_qseq)
    , _tid(record->core.tid)
    , _mtid(record->core.mtid)
    , _query_name(bam1_qname(record))
    , _seq_converted(false)
    , _quality_converted(false)
    , _lib_index(~0)
{
    // FIXME: if *bam1_qual(record) = 0xff, there is no quality string?
    // we should test for that
    if (seq_data) {
        _bam_data.assign(bam1_seq(record),
                bam1_qual(record) + record->core.l_qseq);
    }

    if(uint8_t* tmp = bam_aux_get(record, "RG"))
        _readgroup = bam_aux2Z(tmp);
}

int Read::sam_flag() const {
    return _sam_flag;
}

bool Read::proper_pair() const {
    static int const mask = BAM_FPROPER_PAIR | BAM_FUNMAP | BAM_FMUNMAP | BAM_FPAIRED | BAM_FDUP;
    static int const want = BAM_FPROPER_PAIR | BAM_FPAIRED;
    return (_sam_flag & mask) == want;
}

bool Read::either_unmapped() const {
    return _sam_flag & (BAM_FUNMAP | BAM_FMUNMAP);
}

int Read::pair_overlap() const {
    if (leftmost() && _endPos > _mpos) {
        return _endPos - _mpos;
    }
    else {
        return 0;
    }
}

bool Read::interchrom_pair() const {
    return _tid != _mtid;
}
