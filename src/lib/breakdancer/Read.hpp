#pragma once

#include "LegacyConfig.hpp"

#include <string>
#include <vector>

extern "C" {
    #include <sam.h>
    #include <bam.h>
}

class LegacyConfig;

namespace breakdancer {

enum pair_orientation_flag {
    NA = 0, //NA means not applicable.
    ARP_FF = 1,
    ARP_FR_big_insert = 2,
    ARP_FR_small_insert = 3,
    ARP_RF = 4,
    ARP_RR = 8,
    NORMAL_FR = 18,
    NORMAL_RF = 20,
    ARP_CTX = 32,
    UNMAPPED = 192,
    MATE_UNMAPPED = 64
};

class Read {
public:
    std::string readgroup;
    std::string library;

    Read(bam1_t const* record,
        std::string const& format,
        LegacyConfig const& cfg);

    Read()
        : _record(NULL)
        , _bdflag(NA)
        , _ori(0)
        , _bdqual(0)
        , _abs_isize(0)
    {}

    Read(const Read& other);
    ~Read();

    Read& operator=(const Read& other);

    std::string const& query_name() const;
    std::string const& query_sequence();
    std::string const& quality_string();
    pair_orientation_flag const& bdflag();
    int const& bdqual();
    void set_bdflag(pair_orientation_flag const& new_flag);
    int const& tid();
    int const& pos();
    int const& query_length();
    char const& ori();
    int const& isize();
    int const& abs_isize();

private: // Functions
    std::string _readgroup();
    std::string _library(ConfigMap<std::string, std::string>::type const& readgroup_library);

private: // Data
    bam1_t* _record;
    mutable std::string _query_name;
    mutable bool _query_name_cached;

    std::string _query_sequence;
    bool _query_seq_cached;

    std::string _quality_string;
    bool _quality_string_cached;

    pair_orientation_flag _bdflag;

    int _tid;
    int _pos;
    int _isize;
    int _query_length;
    char _ori;
    int _bdqual;
    int _abs_isize;
    bool _abs_isize_cached;

    std::string platform;
};

inline
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

inline
pair_orientation_flag
determine_bdflag(bam1_t const* record, std::string const& platform) {
    pair_orientation_flag flag = NA;
    int read_reversed = record->core.flag & BAM_FREVERSE;
    int mate_reversed = record->core.flag & BAM_FMREVERSE;

    // this should probably include QCfail as well
    if(!(record->core.flag & BAM_FDUP)) {
        if(record->core.flag & BAM_FPAIRED) {
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
                if(platform == "solid") { //assuming config parser has normalized this for me
                    flag = NORMAL_FR; //normal insert size
                }
                else {
                    if(record->core.pos < record->core.mpos) {
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
                        if(record->core.flag & BAM_FREAD1) {
                            flag = (record->core.pos < record->core.mpos) ? ARP_FR_big_insert : ARP_RF;
                        }
                        else {
                            flag = (record->core.pos > record->core.mpos) ? ARP_FR_big_insert : ARP_RF;
                        }
                    }
                    else {
                        if(record->core.flag & BAM_FREAD1) {
                            flag = (record->core.pos > record->core.mpos) ? ARP_FR_big_insert : ARP_RF;
                        }
                        else {
                            flag = (record->core.pos < record->core.mpos) ? ARP_FR_big_insert : ARP_RF;
                        }
                    }
                }
                else {
                    if( ((read_reversed) && (mate_reversed)) ||
                        (!(read_reversed) && !(mate_reversed))) { //do the mates have the same orientation?

                        flag = (mate_reversed) ? ARP_RR : ARP_FF;
                    }
                    else if((record->core.mpos > record->core.pos && (read_reversed)) || (record->core.pos > record->core.mpos && !(read_reversed))) {
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

}
