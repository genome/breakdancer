#pragma once

#include "Read.hpp"

#include <boost/function.hpp>
#include <boost/iterator/filter_iterator.hpp>

#include <vector>

class BasicRegion {
public:
    typedef breakdancer::Read ReadType;
    typedef std::vector<ReadType> ReadVector;
    typedef boost::filter_iterator<
            boost::function<bool(ReadType const&)>,
            ReadVector::const_iterator
            > const_read_iterator;


    BasicRegion() {}
    BasicRegion(int idx, int chr, int start, int end, int normal_read_pairs)
        : index(idx)
        , chr(chr)
        , start(start)
        , end(end)
        , normal_read_pairs(normal_read_pairs)
    {
    }

    int index;
    int chr;
    int start;
    int end;
    int normal_read_pairs;

    void swap_reads(ReadVector& reads) {
        _reads.swap(reads);
    }

    int size() const {
        return end - start + 1;
    }

    ReadVector const& reads() const {
        return _reads;
    }

    template<typename PredType>
    const_read_iterator reads_begin(PredType pred) const {
        return const_read_iterator(pred, _reads.begin(), _reads.end());
    }

    template<typename PredType>
    const_read_iterator reads_end(PredType pred) const {
        return const_read_iterator(pred, _reads.end(), _reads.end());
    }

private:
    ReadVector _reads;
};
