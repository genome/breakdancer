#include "sam.h"
#include "bam.h"
#include "ksort.h"
#include "khash.h"

#define HEAP_EMPTY 0xffffffffffffffffull
#define __pos_cmp(a, b) ((a).pos > (b).pos || ((a).pos == (b).pos && ((a).i > (b).i || ((a).i == (b).i && (a).idx > (b).idx))))
// define not by query name anyway
#define g_is_by_qname 0	

typedef struct {
	int i;
	uint64_t pos, idx;
	bam1_t *b;
} heap1_t;

static inline int strnum_cmp(const char *a, const char *b)
{
        char *pa, *pb;
        pa = (char*)a; pb = (char*)b;
        while (*pa && *pb) {
	        if (isdigit(*pa) && isdigit(*pb)) {
			long ai, bi;
     			ai = strtol(pa, &pa, 10);
     			bi = strtol(pb, &pb, 10);
     			if (ai != bi) return ai<bi? -1 : ai>bi? 1 : 0;
     		} else {
     			if (*pa != *pb) break;
     		        ++pa; ++pb;
     		}
     	}
     	if (*pa == *pb)
     		return (pa-a) < (pb-b)? -1 : (pa-a) > (pb-b)? 1 : 0;
	return *pa<*pb? -1 : *pa>*pb? 1 : 0;
}

static inline int heap_lt(const heap1_t a, const heap1_t b)
{
        if (g_is_by_qname) {
        int t;
        if (a.b == 0 || b.b == 0) return a.b == 0? 1 : 0;
        	t = strnum_cmp(bam1_qname(a.b), bam1_qname(b.b));
        	return (t > 0 || (t == 0 && __pos_cmp(a, b)));
	} else return __pos_cmp(a, b);
}

KSORT_INIT(heap, heap1_t, heap_lt)

typedef struct {
        uint64_t u, v;
} pair64_t;

#define pair64_lt(a,b) ((a).u < (b).u)
KSORT_INIT(off, pair64_t, pair64_lt)

typedef struct {
	uint32_t m, n;
	pair64_t *list;
} bam_binlist_t;

KHASH_MAP_INIT_INT(i, bam_binlist_t)

typedef struct {
	int32_t n, m;
	uint64_t *offset;
} bam_lidx_t;

//KHASH_MAP_INIT_INT(i, bam_binlist_t)

struct __bam_index_t {
	int32_t n;
	khash_t(i) **index;
	bam_lidx_t *index2;
};

//pair64_t * get_chunk_coordinates(const bam_index_t *idx, int tid, int beg, int end, int* cnt_off);

#define BAM_LIDX_SHIFT    14
#define MAX_BIN 37450 // =(8^6-1)/7+1
static inline int reg2bins(uint32_t beg, uint32_t end, uint16_t list[MAX_BIN])
{
	int i = 0, k;
	--end;
	list[i++] = 0;
	for (k =    1 + (beg>>26); k <=    1 + (end>>26); ++k) list[i++] = k;
	for (k =    9 + (beg>>23); k <=    9 + (end>>23); ++k) list[i++] = k;
	for (k =   73 + (beg>>20); k <=   73 + (end>>20); ++k) list[i++] = k;
	for (k =  585 + (beg>>17); k <=  585 + (end>>17); ++k) list[i++] = k;
	for (k = 4681 + (beg>>14); k <= 4681 + (end>>14); ++k) list[i++] = k;
	return i;
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
			int j;
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


/*static inline int heap_lt(const heap1_t a, const heap1_t b)
{
        //if (g_is_by_qname) {
        //        int t;
        //        if (a.b == 0 || b.b == 0) return a.b == 0? 1 : 0;
        //        t = strnum_cmp(bam1_qname(a.b), bam1_qname(b.b));
        //        return (t > 0 || (t == 0 && __pos_cmp(a, b)));
        //} else 
	return __pos_cmp(a, b);
}*/

static inline int is_overlap(uint32_t beg, uint32_t end, const bam1_t *b)
{
        uint32_t rbeg = b->core.pos;
        uint32_t rend = b->core.n_cigar? bam_calend(&b->core, bam1_cigar(b)) : b->core.pos + 1;
        return (rend > beg && rbeg < end);
}



/*void ks_heapadjust_heap(int i, int n, heap1_t l[])
{
	int k = i;
	heap1_t tmp = l[i];
	while ((k = (k << 1) + 1) < n) {
		if (k != n - 1 && heap_lt(l[k], l[k+1])) ++k;
		if (heap_lt(l[k], tmp)) break;
		l[i] = l[k]; i = k;
	}
	l[i] = tmp;
}
void ks_heapmake_heap(int lsize, heap1_t l[])
{
	int i;
	for (i = (lsize >> 1) - 1; i != (int)(-1); --i)
		ks_heapadjust_heap(i, lsize, l);
}*/
