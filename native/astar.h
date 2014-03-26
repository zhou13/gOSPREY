#ifndef __ASTAR_H
#define __ASTAR_H

#define _WIN32_WINNT 0x0600

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <limits.h>
#include <float.h>
#include <assert.h>
#include <math.h>
#include <inttypes.h>

#include "utility.h"

#define nUSE_OPENCL

#ifndef MAX_ROTAMER
#  define MAX_ROTAMER 700
#endif
#ifndef MAX_LEVEL
#  define MAX_LEVEL 20
#endif
#ifndef MAX_ROUND
#  define MAX_ROUND 120
#endif

extern int tree_level;
extern int rot_cnt;

extern bool do_pert;

extern int node_offset[MAX_LEVEL];
extern int rot_per_level[MAX_LEVEL];

extern bool do_pair_prune;
extern bool **pair_pruned;

extern bool do_triple_prune;
extern bool ***triple_pruned;

extern float reduce_energy[MAX_ROTAMER][MAX_ROTAMER];
extern float pm_energy[MAX_ROTAMER][MAX_LEVEL+1];

extern bool enable_cpu;
extern bool enable_gpu;
extern size_t cpu_memory;
extern size_t gpu_memory;
extern int num_gpu_group;
extern int num_gpu_item;
extern int num_gpu_item2;
extern double shrink_ratio;

typedef unsigned int uint;
typedef unsigned char uchar;
typedef unsigned short ushort;

static inline void error(const char *msg)
{
    fprintf(stderr, "\t**error**: %s\n", msg);
    fflush(stderr);
    exit(EXIT_FAILURE);
}

static inline void warning(const char *msg)
{
    fprintf(stderr, "\t**warning**: %s\n", msg);
}

#ifdef CRASH_DEBUG
#  define PRINT_LINE() { fprintf(stderr, "Here is line %d\n", __LINE__); fflush(stderr); }
#else
#  define PRINT_LINE() ((void)0)
#endif

#ifdef __WIN32
#include <windows.h>

static uint64_t wall_time = 0;

__attribute__ ((unused))
static void wall_time_begin(void)
{
    wall_time = timeGetTime();
}
__attribute__ ((unused))
static uint64_t wall_time_elapsed(void)
{
    return timeGetTime() - wall_time;
}
#endif

#ifdef __unix
#include <sys/time.h>
static uint64_t wall_time = 0;
__attribute__ ((unused))
static void wall_time_begin(void)
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    wall_time = (uint64_t)tv.tv_usec/1000 + (uint64_t)tv.tv_sec*1000;
}
__attribute__ ((unused))
static uint64_t wall_time_elapsed(void)
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    uint64_t now = (uint64_t)tv.tv_usec/1000 + (uint64_t)tv.tv_sec*1000;
    return now - wall_time;
}
#endif

#define UNUSED(x) (void)(x)
#define STATIC_ASSERT(e) typedef char _STATIC_ASSERT[(e)?1:-1]

#define swap(x, y) ({ \
    (void) (&x == &y); \
    typeof(x) z; \
    z = x; \
    x = y; \
    y = z; \
})

#define min(x, y) ({ \
        typeof(x) _min1 = (x); \
        typeof(y) _min2 = (y); \
        (void) (&_min1 == &_min2); \
        _min1 < _min2 ? _min1 : _min2; })

#define max(x, y) ({ \
        typeof(x) _max1 = (x); \
        typeof(y) _max2 = (y); \
        (void) (&_max1 == &_max2); \
        _max1 > _max2 ? _max1 : _max2; })

#define pair(x, y) struct { x a; y b; }

#define __plus(x, y) ((x)+(y))
#define __minus(x, y) ((x)+(y))
#define __multiply(x, y) ((x)*(y))
#define __divide(x, y) ((x)/(y))

#define func_map(begin, end, op) \
    for (typeof(&(*begin)) __func_p = begin; __func_p != end; ++__func_p) \
        *__func_p = op(*__func_p);

#define func_reduce(result, begin, end, op) \
    for (typeof(&(*begin)) __func_p = begin; __func_p != end; ++__func_p) \
        result = op(result, *__func_p);

#endif
