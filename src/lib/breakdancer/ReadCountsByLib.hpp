#pragma once

#include "common/utility.hpp"

#include <stdint.h>

#include <iostream>
#include <functional>
#include <map>
#include <string>

class ReadCountsByLib {
public:
    typedef std::string LibId;
    typedef uint32_t IntType;
    typedef std::map<LibId, IntType> MapType;
    typedef MapType::iterator iterator;
    typedef MapType::const_iterator const_iterator;

    void clear() {
        _counts.clear();
    }

    const_iterator find(LibId const& lib) const {
        return _counts.find(lib);
    }

    iterator find(LibId const& lib) {
        return _counts.find(lib);
    }

    const_iterator begin() const {
        return _counts.begin();
    }

    iterator begin() {
        return _counts.begin();
    }

    const_iterator end() const {
        return _counts.end();
    }

    iterator end() {
        return _counts.end();
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

    size_t size() const {
        return _counts.size();
    }

    bool empty() const {
        return _counts.empty();
    }

    // Operators
    ReadCountsByLib& operator+=(ReadCountsByLib const& rhs) {
        merge_maps(_counts, rhs._counts, std::plus<IntType>());
        return *this;
    }

    bool operator==(ReadCountsByLib const& rhs) const {
        return _counts == rhs._counts;
    }

    ReadCountsByLib& operator-=(ReadCountsByLib const& rhs) {
        using namespace std;
        for (const_iterator i = rhs.begin(); i != rhs.end(); ++i) {
            pair<iterator, bool> inserted = _counts.insert(make_pair(i->first, -i->second));
            if (!inserted.second) {
                inserted.first->second -= i->second;
                if (inserted.first->second == 0)
                    _counts.erase(inserted.first);
            }
        }
        return *this;
    }

    std::ostream& toStream(std::ostream& s) const {
        s << "(";
        for (const_iterator i = begin(); i != end(); ++i) {
            if (i != begin())
                s << ", ";
            s << i->first << ": " << i->second;
        }
        s << ")";
        return s;
    }

private:
    MapType _counts;
};

inline
ReadCountsByLib operator+(ReadCountsByLib lhs, ReadCountsByLib const& rhs) {
    lhs += rhs;
    return lhs;
}


inline
ReadCountsByLib operator-(ReadCountsByLib lhs, ReadCountsByLib const& rhs) {
    lhs -= rhs;
    return lhs;
}

inline
std::ostream& operator<<(std::ostream& s, ReadCountsByLib const& counts) {
    return counts.toStream(s);
}

