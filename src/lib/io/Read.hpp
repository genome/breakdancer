#pragma once

#include "common/ReadFlags.hpp"
#include "common/namespace.hpp"

#include <cassert>
#include <ostream>
#include <string>
#include <vector>

extern "C" {
    #include <sam.h>
    #include <bam.h>
}

class LibraryInfo;

BEGIN_NAMESPACE(breakdancer)

class Read {
public:
    Read(bam1_t const* record, bool seq_data = true);

    Read()
        : _bdflag(NA)
        , _ori(FWD)
        , _abs_isize(0)
        , _bdqual(0)
        , _lib_index(-1)
    {}

    void set_bdflag(pair_orientation_flag const& new_flag);

    std::string const& query_name() const;
    std::string const& query_sequence() const;
    std::string const& quality_string() const;
    std::string const& readgroup() const;
    pair_orientation_flag const& bdflag() const;
    int const& bdqual() const;
    int const& tid() const;
    int const& pos() const;
    int const& query_length() const;
    strand_e const& ori() const;
    int const& isize() const;

    void set_lib_index(std::size_t const& index);
    std::size_t const& lib_index() const;
    int const& abs_isize() const;

    void to_fastq(std::ostream& stream) const;


private: // Data
    pair_orientation_flag _bdflag;
    strand_e _ori;
    int _abs_isize;
    int _bdqual;
    int _isize;
    int _pos;
    int _query_length;
    int _tid;
    std::string _query_name;
    std::string _readgroup;



    mutable std::string _query_sequence;
    mutable bool _seq_converted;

    mutable std::string _quality_string;
    mutable bool _quality_converted;

    std::string _bam_data;

    std::size_t _lib_index;
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


END_NAMESPACE(breakdancer)
