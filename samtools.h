#include "sam.h"
#include "bam.h"
#include "ksort.h"
#include "khash.h"

#define HEAP_EMPTY 0xffffffffffffffffull
//#define __pos_cmp(a, b) ((a).pos > (b).pos || ((a).pos == (b).pos && ((a).i > (b).i || ((a).i == (b).i && (a).idx > (b).idx))))

typedef struct {
	int i;
	uint64_t pos, idx;
	bam1_t *b;
} heap1_t;

typedef struct {
        uint64_t u, v;
} pair64_t;

pair64_t * get_chunk_coordinates(const bam_index_t *idx, int tid, int beg, int end, int* cnt_off);

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
