#pragma once

#include "common/namespace.hpp"

#include <boost/array.hpp>
#include <cassert>

enum strand_e {
    FWD = 0,
    REV = 1
};

enum ReadFlag {
    NA = 0, //NA means not applicable.
    ARP_FF,
    ARP_FR_big_insert,
    ARP_FR_small_insert,
    ARP_RF,
    ARP_RR,
    NORMAL_FR,
    NORMAL_RF,
    ARP_CTX,
    MATE_UNMAPPED,
    UNMAPPED,
    NUM_ORIENTATION_FLAGS
};

template<typename T>
struct PerFlagArray {
    typedef boost::array<T, NUM_ORIENTATION_FLAGS> type;
};

struct FlagValues {
    FlagValues();

    int operator[](ReadFlag idx) const;

private:
    PerFlagArray<int>::type _values;
};

inline
int FlagValues::operator[](ReadFlag idx) const {
    return _values[int(idx)];
}

extern const FlagValues FLAG_VALUES;
