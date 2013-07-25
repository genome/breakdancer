#pragma once

#include "io/Read.hpp"

#include <boost/function.hpp>
#include <boost/iterator/filter_iterator.hpp>
#include <boost/range/iterator_range.hpp>

#include <vector>

class BasicRegion {
public:
    typedef breakdancer::Read ReadType;
    typedef std::vector<ReadType> ReadVector;
    typedef boost::filter_iterator<
            boost::function<bool(ReadType const&)>,
            ReadVector::const_iterator
            > const_read_iterator;

    typedef boost::iterator_range<const_read_iterator> iterator_range;


    BasicRegion() {}
    BasicRegion(int idx, int chr, int start, int end, int normal_read_pairs)
        : index(idx)
        , chr(chr)
        , start(start)
        , end(end)
        , normal_read_pairs(normal_read_pairs)
        , fwd_read_count(0)
        , rev_read_count(0)
        , times_accessed(0)
        , times_collapsed(0)
    {
    }

    int index;
    int chr;
    int start;
    int end;
    int normal_read_pairs;
    int fwd_read_count;
    int rev_read_count;
    int times_accessed;
    int times_collapsed;

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
    iterator_range reads_range(PredType pred) const {
        return make_iterator_range(reads_begin(pred), reads_end(pred));
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
