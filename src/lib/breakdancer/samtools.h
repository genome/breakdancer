#include "saminternals.h"

#include "sam.h"
#include "bam.h"
#include "ksort.h"
#include "khash.h"
#include <cctype>

#define HEAP_EMPTY 0xffffffffffffffffull
#define __pos_cmp(a, b) ((a).pos > (b).pos || ((a).pos == (b).pos && ((a).i > (b).i || ((a).i == (b).i && (a).idx > (b).idx))))
// define not by query name anyway
#define g_is_by_qname 0    

namespace {
    int strnum_cmp(const char *a, const char *b)
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

    int heap_lt(const heap1_t a, const heap1_t b)
    {
            if (g_is_by_qname) {
            int t;
            if (a.b == 0 || b.b == 0) return a.b == 0? 1 : 0;
                t = strnum_cmp(bam1_qname(a.b), bam1_qname(b.b));
                return (t > 0 || (t == 0 && __pos_cmp(a, b)));
        } else return __pos_cmp(a, b);
    }

}
