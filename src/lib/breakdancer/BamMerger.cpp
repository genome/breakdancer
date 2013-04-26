#include "BamMerger.hpp"

#include <cassert>
#include <stdexcept>

using namespace std;

BamMerger::Stream::Stream()
    : bam(0)
    , entry(0)
    , status(-1)
{
}

BamMerger::Stream::Stream(IBamReader* bam)
    : bam(bam)
    , entry(bam_init1())
{
    advance();
}

BamMerger::Stream::~Stream() {
    if (entry)
        bam_destroy1(entry);
    entry = 0;
}

bool BamMerger::Stream::operator>(Stream const& rhs) const {
    assert(valid() && rhs.valid());
    bam1_t const* x = entry;
    bam1_t const* y = rhs.entry;

    // TODO: determine if it is better to do this samtools style
    // (by cramming tid, pos, and strand into a uint64_t like:
    // (tid << 32) | (pos << 1) | strand)
    if (x->core.tid > y->core.tid)
        return true;

    if (y->core.tid > x->core.tid)
        return false;

    if (x->core.pos > y->core.pos)
        return true;

    if (y->core.pos > x->core.pos)
        return false;

    return bam1_strand(x) > bam1_strand(y);
}

bool BamMerger::Stream::valid() const {
    return status > 0;
}

bool BamMerger::Stream::advance() {
    assert(bam && entry);
    status = bam->next(entry);
    return status > 0;
}

BamMerger::BamMerger(std::vector<IBamReader*> const& streams) {
    if (streams.empty())
        throw runtime_error("BamMerger created with no input streams!");

    // TODO: validate that headers are all "compatible".
    _header = streams[0]->header();

    typedef std::vector<IBamReader*>::const_iterator IterType;
    for (IterType i = streams.begin(); i != streams.end(); ++i) {
        Stream* s = new Stream(*i);
        if (s->valid()) {
            _streams.push(s);
        } else {
            delete s;
        }
    }
}

BamMerger::~BamMerger() {
    while (!_streams.empty()) {
        delete _streams.top();
        _streams.pop();
    }
}

bam_header_t* BamMerger::header() const {
    return _header;
}

int BamMerger::next(bam1_t* entry) {
    if (_streams.empty())
        return -1;

    Stream* s = _streams.top();
    _streams.pop();

    assert(s && s->entry && s->status > 0);

    int rv = s->status;
    bam_copy1(entry, s->entry);

    if (s->advance()) {
        _streams.push(s);
    } else {
        delete s;
    }

    return rv;
}
