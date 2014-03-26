#pragma once

#include "array.h"
#include <assert.h>

#define priority_queue array

#define pq_init(heap) \
({ \
    arr_init(heap); \
    arr_resize(heap, 11); \
})

#define pq_push(heap, val, cmp) \
({ \
    int __now = ++heap.size; \
    typeof(val) __val = val; \
    arr_resize(heap, heap.size+1); \
    while (__now > 1) { \
        int __next = __now / 2; \
        if (cmp(__val, heap.v[__next])) { \
            heap.v[__now] = heap.v[__next]; \
            __now = __next; \
        } else \
            break; \
    } \
    heap.v[__now] = __val; \
})

#define pq_empty(heap) (arr_empty(heap))

#define pq_top(heap) (heap.v[1])

#define pq_pop(heap, cmp) \
({ \
    assert(heap.size); \
    typeof(heap.v[0]) __ret = heap.v[1]; \
    typeof(heap.v[0]) __val = heap.v[heap.size--]; \
    int __now = 1; \
 \
    int __next; \
    while ((__next = __now*2) <= heap.size) { \
        if (__next+1 <= heap.size && cmp(heap.v[__next+1], heap.v[__next])) \
            ++__next; \
        if (cmp(heap.v[__next], __val)) { \
            heap.v[__now] = heap.v[__next]; \
            __now = __next; \
        } else \
            break; \
    } \
    heap.v[__now] = __val; \
    __ret; \
})

#define pq_for(i, x) for (typeof((x).v) i = (x).v+1; i <= (x).v + (x).size; ++i)
