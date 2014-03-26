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

#pragma once

#ifdef __GNUC__
#  define NUNUSED __attribute__ ((unused))
#else
#  define NUNUSED
#endif


#define QSORT_CMP_INT(x, y) (*(x) < *(y))
#define QSORT_CMP_STR(x, y) (strcmp(*(x), *(y)) < 0)
#define QSORT_CMP_DBL(x, y) (*(x) < *(y))
#define QSORT_INIT_INT QSORT_INIT(int, QSORT_CMP_INT, int);
#define QSORT_INIT_STR QSORT_INIT(char *, QSORT_CMP_STR, str);
#define QSORT_INIT_DBL QSORT_INIT(double, QSORT_CMP_DBL, dbl);

#define QSORT_SWAP(x, y) { typeof(*x) z = *x; *x = *y; *y = z; }
#define QSORT_MAGIC_NUMBER 16

#define QSORT_INIT(type, cmp, suffix) \
static void __swap_median_first_##suffix(type *m1, type *m2, type *m3) \
{ \
    if (cmp(*m1, *m2)) { \
        if (cmp(*m2, *m3)) { \
            QSORT_SWAP(m1, m2); \
        } else { \
            if (cmp(*m1, *m3)) { \
                QSORT_SWAP(m1, m3); \
            } \
        } \
    } else { \
        if (cmp(*m3, *m1)) { \
            if (cmp(*m2, *m3)) { \
                QSORT_SWAP(m1, m3); \
            } else { \
                QSORT_SWAP(m1, m2); \
            } \
        } \
    } \
} \
 \
static void __quick_sort_##suffix(type *l, type *r) \
{ \
    while (r - l > QSORT_MAGIC_NUMBER) { \
        __swap_median_first_##suffix(l, l+(r-l)/2, r-1); \
 \
        type *i = l+1; \
        type *j = r; \
        for (;;) { \
            while (cmp(*i, *l)) \
                ++i; \
            --j; \
            while (cmp(*l, *j)) \
                --j; \
            if (i >= j) \
                break; \
            QSORT_SWAP(i, j); \
            ++i; \
        } \
        __quick_sort_##suffix(i, r); \
        r = i; \
    } \
} \
 \
static void __insertion_sort_guarded_##suffix(type *l, type *r) \
{ \
    for (type *i = l+1; i < r; ++i) { \
        type t = *i; \
        type *j = i-1; \
        while (j >= l && cmp(t, *j)) { \
            *(j+1) = *j; \
            --j; \
        } \
        *(j+1) = t; \
    } \
} \
 \
static void __insertion_sort_ungarded_##suffix(type *l, type *r) \
{ \
    for (type *i = l; i < r; ++i) { \
        type t = *i; \
        type *j = i-1; \
        while (cmp(t, *j)) { \
            *(j+1) = *j; \
            --j; \
        } \
        *(j+1) = t; \
    } \
} \
 \
static void __insertion_sort_##suffix(type *l, type *r) \
{ \
    if (r - l > QSORT_MAGIC_NUMBER) { \
        __insertion_sort_guarded_##suffix(l, l + QSORT_MAGIC_NUMBER); \
        __insertion_sort_ungarded_##suffix(l + QSORT_MAGIC_NUMBER, r); \
    } else { \
        __insertion_sort_guarded_##suffix(l, r); \
    } \
} \
 \
static void NUNUSED qsort_##suffix(type *l, type *r) \
{ \
    __quick_sort_##suffix(l, r); \
    __insertion_sort_##suffix(l, r); \
}
