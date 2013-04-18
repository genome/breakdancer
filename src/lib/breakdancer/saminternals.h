#pragma once

#include <sam.h>
#include <stdint.h>

typedef struct {
    int i;
    uint64_t pos, idx;
    bam1_t *b;
} heap1_t;

typedef struct {
        uint64_t u, v;
} pair64_t;

typedef struct {
    uint32_t m, n;
    pair64_t *list;
} bam_binlist_t;

typedef struct {
    int32_t n, m;
    uint64_t *offset;
} bam_lidx_t;


