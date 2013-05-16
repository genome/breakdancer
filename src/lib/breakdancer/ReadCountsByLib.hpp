#pragma once

#include "utility.hpp"

#include <stdint.h>

#include <functional>
#include <map>
#include <string>

class ReadCountsByLib {
public:
    typedef std::string LibId;
    typedef uint32_t IntType;
    typedef std::map<LibId, IntType> MapType;

    void increment(LibId const& lib, IntType n) {
        _counts[lib] += n;
    }

    IntType& operator[](LibId const& lib) {
        return _counts[lib];
    }

    IntType& at(LibId const& lib) {
        return _counts.at(lib);
    }

    IntType const& at(LibId const& lib) const {
        return _counts.at(lib);
    }

    ReadCountsByLib& operator+=(ReadCountsByLib const& rhs) {
        merge_maps(_counts, rhs._counts, std::plus<IntType>());
        return *this;
    }

    size_t size() const {
        return _counts.size();
    }

    bool empty() const {
        return _counts.empty();
    }

private:
    MapType _counts;
};

inline
ReadCountsByLib operator+(ReadCountsByLib lhs, ReadCountsByLib const& rhs) {
    lhs += rhs;
    return lhs;
}
