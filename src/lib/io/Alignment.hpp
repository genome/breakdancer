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

uint8_t determine_bdqual(bam1_t const* record);
std::string determine_read_group(bam1_t const* record);

class Alignment {
public:
    Alignment();
    Alignment(bam1_t const* record, bool seq_data = true);

    void set_bdflag(ReadFlag const& new_flag);
    bool proper_pair() const;
    bool either_unmapped() const;
    bool interchrom_pair() const;

    bool has_sequence() const;

    std::string const& query_name() const;

    ReadFlag bdflag() const;
    uint16_t sam_flag() const;
    int32_t tid() const;
    int32_t pos() const;
    int32_t query_length() const;
    int32_t abs_isize() const;
    uint8_t bdqual() const;
    strand_e ori() const;

    void set_lib_index(std::size_t const& index);
    std::size_t lib_index() const;

    void to_fastq(std::ostream& stream) const;
    bool leftmost() const;

private: // Data
    int32_t _tid;
    int32_t _pos;
    int32_t _query_length;
    int32_t _mtid;
    int32_t _mpos;
    int32_t _abs_isize;
    uint16_t _sam_flag;
    uint8_t _bdqual;

    std::string _query_name;

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
ReadFlag Alignment::bdflag() const {
    return _bdflag;
}

inline
uint8_t Alignment::bdqual() const {
    return _bdqual;
}

inline
int32_t Alignment::tid() const {
    return _tid;
}

inline
int32_t Alignment::pos() const {
    return _pos;
}

inline
int32_t Alignment::query_length() const {
    return _query_length;
}

inline
strand_e Alignment::ori() const {
    return sam_flag() & BAM_FREVERSE ? REV : FWD;
}

inline
int32_t Alignment::abs_isize() const {
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
std::size_t Alignment::lib_index() const {
    assert(_lib_index != ~0ull);
    return _lib_index;
}

inline
bool Alignment::has_sequence() const {
    return !_bam_data.empty() && query_length() > 0;
}

inline
uint16_t Alignment::sam_flag() const {
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
