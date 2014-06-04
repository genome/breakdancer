#pragma once

#include "common/ReadFlags.hpp"

#include <cassert>
#include <ostream>
#include <stdint.h>
#include <string>
#include <vector>

extern "C" {
    #include <sam.h>
    #include <bam.h>
}

struct LibraryInfo;

ReadFlag determine_bdflag(bam1_t const* record);

class Read {
public:
    Read(bam1_t const* record, bool seq_data = true);

    Read()
        : _sam_flag(0)
        , _bdflag(ReadFlag::NA)
        , _ori(FWD)
        , _abs_isize(0)
        , _bdqual(0)
        , _lib_index(-1)
    {}

    void set_bdflag(ReadFlag const& new_flag);
    bool proper_pair() const;
    bool either_unmapped() const;
    bool interchrom_pair() const;

    std::string const& query_name() const;
    std::string const& query_sequence() const;
    std::string const& quality_string() const;
    std::string const& readgroup() const;
    ReadFlag const& bdflag() const;
    int sam_flag() const;
    int const& bdqual() const;
    int const& tid() const;
    int const& pos() const;
    int const& query_length() const;
    strand_e const& ori() const;
    int const& isize() const;

    void set_lib_index(std::size_t const& index);
    std::size_t const& lib_index() const;
    int const& abs_isize() const;
    int pair_overlap() const;

    void to_fastq(std::ostream& stream) const;
    bool leftmost() const;


private: // Data
    int _sam_flag;
    ReadFlag _bdflag;
    strand_e _ori;
    int _abs_isize;
    int _bdqual;
    int _isize;
    int _pos;
    int _mpos;
    int _endPos;
    int _query_length;
    int _tid;
    int _mtid;
    std::string _query_name;
    std::string _readgroup;



    mutable std::string _query_sequence;
    mutable bool _seq_converted;

    mutable std::string _quality_string;
    mutable bool _quality_converted;

    std::vector<uint8_t> _bam_data;

    std::size_t _lib_index;
};

inline
bool Read::leftmost() const {
    return _pos < _mpos;
}

inline
void Read::set_bdflag(ReadFlag const& new_flag) {
    _bdflag = new_flag;
}

inline
ReadFlag const& Read::bdflag() const {
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
strand_e const& Read::ori() const {
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
        uint8_t const* qdata = _bam_data.data() + ((_query_length+1) >> 1);
        if (qdata[0] != 0xff) {
            _quality_string.reserve(_query_length);
            for (int i = 0; i < _query_length; ++i)
                _quality_string += char(qdata[i] + 33);
        }
        _quality_converted = true;
    }

    return _quality_string;
}

inline
void Read::set_lib_index(std::size_t const& index) {
    _lib_index = index;
}


inline
std::size_t const& Read::lib_index() const {
    assert(_lib_index != ~0ull);
    return _lib_index;
}

inline
void Read::to_fastq(std::ostream& stream) const {
    stream << "@" << query_name() << "\n"
        << query_sequence()
        << "\n+\n"
        << quality_string() << "\n";
}

inline
std::string const& Read::readgroup() const {
    return _readgroup;
}
