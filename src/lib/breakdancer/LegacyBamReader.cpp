#include "LegacyBamReader.hpp"

#include <iostream>

#define BAM_LIDX_SHIFT    14
#define MAX_BIN 37450 // =(8^6-1)/7+1
#define pair64_lt(a,b) ((a).u < (b).u)

using namespace std;

KHASH_MAP_INIT_INT(i, bam_binlist_t)
struct __bam_index_t {
    int32_t n;
    khash_t(i) **index;
    bam_lidx_t *index2;
};

namespace {
    KSORT_INIT(off, pair64_t, pair64_lt)

    int is_overlap(uint32_t beg, uint32_t end, const bam1_t *b)
    {
        uint32_t rbeg = b->core.pos;
        uint32_t rend = b->core.n_cigar ?
            bam_calend(&b->core, bam1_cigar(b))
            : b->core.pos + 1;
        return (rend > beg && rbeg < end);
    }

    int reg2bins(uint32_t beg, uint32_t end, uint16_t list[MAX_BIN])
    {
        uint32_t i = 0, k;
        --end;
        list[i++] = 0;
        for (k =    1 + (beg>>26); k <=    1 + (end>>26); ++k) list[i++] = k;
        for (k =    9 + (beg>>23); k <=    9 + (end>>23); ++k) list[i++] = k;
        for (k =   73 + (beg>>20); k <=   73 + (end>>20); ++k) list[i++] = k;
        for (k =  585 + (beg>>17); k <=  585 + (end>>17); ++k) list[i++] = k;
        for (k = 4681 + (beg>>14); k <= 4681 + (end>>14); ++k) list[i++] = k;
        return i;
    }
}

// bam_fetch helper function retrieves 
pair64_t * get_chunk_coordinates(const bam_index_t *idx, int tid, int beg, int end, int* cnt_off)
{
    uint16_t *bins;
    int i, n_bins, n_off;
    pair64_t *off;
    khint_t k;
    khash_t(i) *index;
    uint64_t min_off;

    bins = (uint16_t*)calloc(MAX_BIN, 2);
    n_bins = reg2bins(beg, end, bins);
    index = idx->index[tid];
    min_off = (beg>>BAM_LIDX_SHIFT >= idx->index2[tid].n)? 0 : idx->index2[tid].offset[beg>>BAM_LIDX_SHIFT];
    for (i = n_off = 0; i < n_bins; ++i) {
        if ((k = kh_get(i, index, bins[i])) != kh_end(index))
            n_off += kh_value(index, k).n;
    }
    if (n_off == 0) {
        free(bins); return 0;
    }
    off = (pair64_t*)calloc(n_off, 16);
    for (i = n_off = 0; i < n_bins; ++i) {
        if ((k = kh_get(i, index, bins[i])) != kh_end(index)) {
            uint32_t j;
            bam_binlist_t *p = &kh_value(index, k);
            for (j = 0; j < p->n; ++j)
                if (p->list[j].v > min_off) off[n_off++] = p->list[j];
        }
    }
    free(bins);
    {
        bam1_t *b = (bam1_t*)calloc(1, sizeof(bam1_t));
        int l;
        ks_introsort(off, n_off, off);
        // resolve completely contained adjacent blocks
        for (i = 1, l = 0; i < n_off; ++i)
            if (off[l].v < off[i].v)
                off[++l] = off[i];
        n_off = l + 1;
        // resolve overlaps between adjacent blocks; this may happen due to the merge in indexing
        for (i = 1; i < n_off; ++i)
            if (off[i-1].v >= off[i].u) off[i-1].v = off[i].u;
        { // merge adjacent blocks
#if defined(BAM_TRUE_OFFSET) || defined(BAM_VIRTUAL_OFFSET16)
            for (i = 1, l = 0; i < n_off; ++i) {
#ifdef BAM_TRUE_OFFSET
                if (off[l].v + BAM_MIN_CHUNK_GAP > off[i].u) off[l].v = off[i].v;
#else
                if (off[l].v>>16 == off[i].u>>16) off[l].v = off[i].v;
#endif
                else off[++l] = off[i];
            }
            n_off = l + 1;
#endif
        }
        bam_destroy1(b);
    }
    *cnt_off = n_off;
    return off;
}



// prepare: read bam file by a chromosome
pair64_t* ReadBamChr_prep(
    string const& chr_str,
    string const& bam_name,
    int *tid,
    int *beg,
    int *end,
    samfile_t *in,
    int *n_off)
{
    pair64_t *off;
    bam_index_t *idx;
    idx = bam_index_load(bam_name.c_str());// index
    if(idx == 0){
        off = (pair64_t*)calloc(1, 16);
        off[0].u = uint64_t(-1);
        off[0].v = uint64_t(-1);
        cout << "Error: should do sort and index first if specifying chromosome!" << endl;
        return off;
    }
    char *chr_str_;
    chr_str_ = new char[chr_str.length()+1];
    strcpy(chr_str_, chr_str.c_str());
    bam_parse_region(in->header, chr_str_, tid, beg, end);// parse
    delete []chr_str_;

    // return the file handle for handle
    //*fp = in->x.bam;
    //bamFile fp = in->x.bam;
    off = get_chunk_coordinates(idx, *tid, *beg, *end, n_off);
    bam_index_destroy(idx);
    return off;
}

// read bam file by a chromosome by one line; fp will track where we are
int ReadBamChr(
    bam1_t *b,
    bamFile fp,
    int tid,
    int beg,
    int end,
    uint64_t *curr_off,
    int *i,
    int *n_seeks,
    pair64_t *off,
    int n_off)
{
    if (off == 0) return 0;

    if (*curr_off == 0 || (*i>=0 && *curr_off >= off[*i].v)) { // then jump to the next chunk
        if (*i == n_off - 1) return 0; // no more chunks
        if (*i >= 0) return 0;//assert(*curr_off == off[*i].v); // otherwise bug
        if (*i < 0 || (*i>=0 && off[*i].v != off[*i+1].u)) { // not adjacent chunks; then seek
            if(*i + 1 >= 0){
                bam_seek(fp, off[*i+1].u, SEEK_SET);
                *curr_off = bam_tell(fp);
                ++(*n_seeks);
            }
        }
        ++(*i);
    }
    if (bam_read1(fp, b) > 0) {
        *curr_off = bam_tell(fp);
        if (b->core.tid != tid || b->core.pos >= end) return 0; // no need to proceed
        else if (is_overlap(beg, end, b)) return 1;
    }
    else
        return 0;
    return 1;
}
