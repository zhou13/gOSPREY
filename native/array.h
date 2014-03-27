/* The MIT License (MIT)
 *
 * Copyright (c) 2012, Yichao Zhou
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

/* This is a macro library to implement a dynamic-sizegth array.  Notice that
 * arr_pop will *NOT* shrink the actual memory it uses. */

#pragma once

#define array(type) struct { \
    type *v; \
    int size; \
    int mem_size; \
}

#define arr_init(x) ({ \
    (x).v = NULL; \
    (x).size = 0; \
    (x).mem_size = 0; \
})

#define arr_push(x, y) \
({ \
    if ((x).size == (x).mem_size) { \
        (x).mem_size = (x).mem_size * 2 + 2; \
        (x).v = (typeof((x).v))realloc((x).v, sizeof((x).v[0]) * (x).mem_size); \
    } \
    (x).v[(x).size++] = (y); \
})

#define arr_empty(x) (x.size == 0)

#define arr_resize(x, new_size) \
({ \
    if ((x).mem_size < (new_size)) { \
        (x).mem_size = max((x).mem_size*2, (new_size)); \
        (x).v = (typeof((x).v))realloc((x).v, sizeof((x).v[0]) * (x).mem_size); \
    } \
})

#define arr_pop(x)   ({(x).size--; (x).v[(x).size];})
#define arr_free(x)  (free((x).v))
#define arr_back(x)  ((x).v[(x).size-1])
#define arr_front(x) ((x).v[0])

#define arr_for(i, x) for (typeof((x).v) i = (x).v; i < (x).v + (x).size; ++i)

