#pragma once

#include <bam.h>

// This struct exists to enable RAII for bam1_t* objects.
// It should be interchangeable for most purposes.
// Note that much (except cast operator?) can be done by C++11's
// std::unique_ptr, but we're not using C++11 yet.
// We definitely don't want to use shared_ptr as that imposes additional
// overhead.
class RawBamEntry {
public:
    RawBamEntry();
    ~RawBamEntry();

    // cast conversion operator allows passing RawBamEntry where bam1_t* is
    // expected.
    operator bam1_t const*() const;
    operator bam1_t*();

    // const and non-const pointer operators allow transparent usage like
    //  b->core.tid
    //  (*b).core.tid
    bam1_t const* operator->() const;
    bam1_t* operator->();
    bam1_t const& operator*() const;
    bam1_t& operator*();

private:
    bam1_t* entry;
};


inline
RawBamEntry::RawBamEntry()
    : entry(bam_init1())
{
}

inline
RawBamEntry::~RawBamEntry() {
    bam_destroy1(entry);
    entry = 0;
}

inline
RawBamEntry::operator bam1_t const*() const {
    return entry;
}

inline
RawBamEntry::operator bam1_t*() {
    return entry;
}

inline
bam1_t* RawBamEntry::operator->() {
    return entry;
}

inline
bam1_t const* RawBamEntry::operator->() const{
    return entry;
}

inline
bam1_t& RawBamEntry::operator*() {
    return *entry;
}

inline
bam1_t const& RawBamEntry::operator*() const{
    return *entry;
}
