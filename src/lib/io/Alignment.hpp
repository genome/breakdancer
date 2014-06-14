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
int determine_bdqual(bam1_t const* record);
std::string determine_read_group(bam1_t const* record);

class Alignment {
public:
    Alignment(bam1_t const* record, bool seq_data = true);

    Alignment()
        : _sam_flag(0)
        , _abs_isize(0)
        , _bdqual(0)
        , _lib_index(-1)
        , _bdflag(ReadFlag::NA)
    {}

    void set_bdflag(ReadFlag const& new_flag);
    bool proper_pair() const;
    bool either_unmapped() const;
    bool interchrom_pair() const;

    bool has_sequence() const;
    std::string const& query_name() const;
    std::string const& readgroup() const;
    ReadFlag const& bdflag() const;
    int sam_flag() const;
    int bdqual() const;
    int tid() const;
    int pos() const;
    int query_length() const;
    strand_e ori() const;

    void set_lib_index(std::size_t const& index);
    std::size_t const& lib_index() const;
    int abs_isize() const;
    int pair_overlap() const;

    void to_fastq(std::ostream& stream) const;
    bool leftmost() const;

private: // Data
    int _sam_flag;
    int _abs_isize;
    int _bdqual;
    int _pos;
    int _mpos;
    int _query_length;
    int _tid;
    int _mtid;
    std::string _query_name;
    std::string _readgroup;

    std::vector<uint8_t> _bam_data;
    std::size_t _lib_index;
    ReadFlag _bdflag;
};

inline
bool Alignment::leftmost() const {
    return _pos < _mpos;
}

inline
void Alignment::set_bdflag(ReadFlag const& new_flag) {
    _bdflag = new_flag;
}

inline
ReadFlag const& Alignment::bdflag() const {
    return _bdflag;
}

inline
int Alignment::bdqual() const {
    return _bdqual;
}

inline
int Alignment::tid() const {
    return _tid;
}

inline
int Alignment::pos() const {
    return _pos;
}

inline
int Alignment::query_length() const {
    return _query_length;
}

inline
strand_e Alignment::ori() const {
    return sam_flag() & BAM_FREVERSE ? REV : FWD;
}

inline
int Alignment::abs_isize() const {
    return _abs_isize;
}

inline
std::string const& Alignment::query_name() const {
    return _query_name;
}

inline
void Alignment::set_lib_index(std::size_t const& index) {
    _lib_index = index;
}


inline
std::size_t const& Alignment::lib_index() const {
    assert(_lib_index != ~0ull);
    return _lib_index;
}

inline
std::string const& Alignment::readgroup() const {
    return _readgroup;
}

inline
bool Alignment::has_sequence() const {
    return !_bam_data.empty() && query_length() > 0;
}

inline
int Alignment::sam_flag() const {
    return _sam_flag;
}

inline
bool Alignment::proper_pair() const {
    static int const mask = BAM_FPROPER_PAIR | BAM_FUNMAP | BAM_FMUNMAP | BAM_FPAIRED | BAM_FDUP;
    static int const want = BAM_FPROPER_PAIR | BAM_FPAIRED;
    return (sam_flag() & mask) == want;
}

inline
bool Alignment::either_unmapped() const {
    return sam_flag() & (BAM_FUNMAP | BAM_FMUNMAP);
}

inline
bool Alignment::interchrom_pair() const {
    return _tid != _mtid;
}
