#pragma once

#include "LegacyConfig.hpp"
#include "namespace.hpp"

#include <string>
#include <vector>

extern "C" {
    #include <sam.h>
    #include <bam.h>
}

class LegacyConfig;

BEGIN_NAMESPACE(breakdancer)

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

pair_orientation_flag determine_bdflag(bam1_t const* record, std::string const& platform);

class Read {
public:
    std::string readgroup;
    std::string library;

    Read(bam1_t const* record,
        std::string const& format,
        LegacyConfig const& cfg);

    Read()
        : _bdflag(NA)
        , _ori(0)
        , _abs_isize(0)
        , _bdqual(0)
    {}

    void set_bdflag(pair_orientation_flag const& new_flag);

    std::string const& query_name() const;
    std::string const& query_sequence() const;
    std::string const& quality_string() const;
    pair_orientation_flag const& bdflag() const;
    int const& bdqual() const;
    int const& tid() const;
    int const& pos() const;
    int const& query_length() const;
    char const& ori() const;
    int const& isize() const;
    int const& abs_isize() const;

private: // Data
    pair_orientation_flag _bdflag;

    char _ori;
    int _abs_isize;
    int _bdqual;
    int _isize;
    int _pos;
    int _query_length;
    int _tid;
    std::string _query_name;


    mutable std::string _query_sequence;
    mutable bool _seq_converted;

    mutable std::string _quality_string;
    mutable bool _quality_converted;

    std::string _bam_data;
};

inline
void Read::set_bdflag(pair_orientation_flag const& new_flag) {
    _bdflag = new_flag;
}

inline
pair_orientation_flag const& Read::bdflag() const {
    return _bdflag;
}

inline
int const& Read::bdqual() const {
    return _bdqual;
}

inline
int const& Read::tid() const {
    return _tid;
}

inline
int const& Read::pos() const {
    return _pos;
}

inline
int const& Read::query_length() const {
    return _query_length;
}

inline
char const& Read::ori() const {
    return  _ori;
}

inline
int const& Read::isize() const {
    return _isize;
}

inline
int const& Read::abs_isize() const {
    return _abs_isize;
}

inline
std::string const& Read::query_name() const {
    return _query_name;
}

inline
std::string const& Read::query_sequence() const {
    if (!_seq_converted) {
        _query_sequence.reserve(_query_length);
        for (int i = 0; i < _query_length; ++i) {
            _query_sequence += bam_nt16_rev_table[bam1_seqi(_bam_data.data(), i)];
        }
        _seq_converted = true;
    }
    return _query_sequence;
}

inline
std::string const& Read::quality_string() const {
    if (!_quality_converted) {
        char const* qdata = _bam_data.data() + ((_query_length+1) >> 1);
        if (qdata[0] != 0xff) {
            _quality_string.reserve(_query_length);
            for (int i = 0; i < _query_length; ++i)
                _quality_string += char(qdata[i] + 33);
        }
        _quality_converted = true;
    }

    return _quality_string;
}

END_NAMESPACE(breakdancer)
