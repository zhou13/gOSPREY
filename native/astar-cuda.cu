/* 
 * If you think this code is a mess... Well, this is originally written in C
 * but then I found that nvcc does not support C99.  So now it is a frank which
 * combines custom either c and cpp.
 *
 * To successor:  Please write detailed comment on function header and global
 * variable!
 */
#define __STDC_LIMIT_MACROS
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_profiler_api.h>

#include <thrust/device_ptr.h>
#include <thrust/host_vector.h>
#include <thrust/sort.h>
#include <thrust/scan.h>

#include <algorithm>
#include <stdint.h>

#include "pq.h"

#include "astar-cuda.h"

#define nVERBOSE_DEBUG

typedef struct data_t {
    float f_score;
    float g_score;
    int level;
    uchar conf[MAX_LEVEL];
} data_t;

typedef struct heap_t {
    float v;
    int idx;
} heap_t;

bool operator< (const data_t &a, const data_t &b)
{
    return a.f_score < b.f_score;
}

// number of data_t s on GPU
static int data_capacity;
// number of data_t s reserved on GPU before for ``shrink''
static int data_reserved;

// GPU parameter
static int num_block;
static int num_local;
static int num_local2;
static int num_global; // =num_block*num_local
static int heap_capacity;  // =data_capacity/num_global
static __device__ __constant__ int d_heap_capacity;

// priority_queue on CPU to store the result return by GPU
#define heap_cmp(a, b) ((a)->f_score < (b)->f_score)
static priority_queue(data_t *) heap;

// GPU heap, contain point to d_data
// Heap of heap: d_heap[heap_capacity], d_heap[2*heap_capacity], ...
static heap_t *d_heap;
// Memory pool
static data_t *d_data;
// Used by shrink function
static data_t *d_data2;
// Answer is put in this array so that CPU can copy.
static data_t *d_output;

// Used by shrink function for sorting
static float *d_data_val;
// Used by shrink function
static int   *d_data_used;

// d_heap_size[i]: heap size of ith heap
static int *d_heap_size;
// d_parent[i]: the node preparing to extend
static int *d_parent;
// radix sort by level of node to decrease branch divergence
static int *d_radix;
// node we need to calculate energy
static int2 *d_input;
// number of child for d_parent[i].
static int *d_node_cnt;

// these point to a number, not a array
static int *d_data_size;
static int *d_output_size;
static int *d_begin_index;
static int *d_begin_index2;
static uint *d_optimal;

// necessary to compute energy
static int *d_node_offset;
static int *d_rot_per_level;
static float *d_self_energy;
static float *d_reduce_energy;
static float *d_pm_energy;

static int rounds = 0;
static uint optimal = 0;
static int curr_conf[MAX_LEVEL];
static int h_num_child = 0;
static int h_output_size = 0;
static int h_data_size = 0;

static int max_data;
static float throw_min;

// map a float number to a unsigned such that:
//     if a < b, flip_float(a) < flip_float(b)
    __host__ __device__
uint flip_float(float fl)
{
    union {
        float fl;
        int  u;
    } un;
    un.fl = fl;
    return un.u ^ ((un.u >> 31) | 0x80000000);
}

    __attribute__((unused))
float reverse_flip_float(int u)
{
    union {
        float f;
        int u;
    } un;
    un.u = u ^ ((~u >> 31) | 0x80000000);
    return un.f;
}

#define CHK_CUDA(exp) \
    if (1)  { \
        cudaError_t v = (exp); \
        if (v != cudaSuccess) { \
            fprintf(stderr, "CUDA ERR [File %s, Line: %d]: %s\n", \
                    __FILE__, __LINE__, cudaGetErrorString(v)); \
            exit(EXIT_FAILURE); \
        } \
    } else


#define copy_from_device(d_ptr0) ({ \
                                  typeof(d_ptr0) d_ptr = d_ptr0; \
                                  typeof(*d_ptr) h; \
                                  CHK_CUDA(cudaMemcpy(&h, d_ptr, sizeof(h), cudaMemcpyDeviceToHost)); \
                                  h; });

    template <typename T>
inline void copy_to_device(T *d_ptr, const T &val)
{
    CHK_CUDA(cudaMemcpy(d_ptr, &val, sizeof(T), cudaMemcpyHostToDevice));
}

    template <typename T>
inline void cuda_free(T *&d_ptr)
{
    if (d_ptr)
        CHK_CUDA(cudaFree(d_ptr));
    d_ptr = NULL;
}

// compute minimal energy for ``data_t *data''
    __device__
void __compute_score(
        int tree_level,
        int rot_cnt,
        data_t *data,
        float old_g_score,
        const int node_offset[],
        const int rot_per_level[],
        const float self_energy[],
        const float reduce_energy[],
        const float pm_energy[])
{
    int level = data->level;
    float g_score = old_g_score; /* use previous g score */
    float h_score = 0.f;

    /* compute g delta */
    {
        int idx = node_offset[level] + data->conf[level];
        for (int j = 0; j < level; ++j) {
            int cord = idx*rot_cnt + node_offset[j]+data->conf[j];
            g_score += reduce_energy[cord];
        }
        g_score += self_energy[idx];
    }

    /* compute h */
    {
        for (int i = level+1; i < tree_level; ++i) {
            float min_energy = FLT_MAX;
            for (int j = 0; j < rot_per_level[i]; ++j) {
                int idx = node_offset[i]+j;
                float cur_energy = self_energy[idx];

                for (int k = 0; k <= level; ++k) {
                    int cord = (idx*rot_cnt) + (node_offset[k]+data->conf[k]);
                    cur_energy += reduce_energy[cord];
                }
                /*
                   int k;
                   for (k = 0; k+1 <= level; k += 2) {
                   int cord0 = (idx*rot_cnt) + (node_offset[k]+data->conf[k]);
                   int cord1 = (idx*rot_cnt) + (node_offset[k+1]+data->conf[k+1]);
                   t0 += reduce_energy[cord0];
                   t1 += reduce_energy[cord1];
                   }
                   if (k <= level) {
                   int cord0 = (idx*rot_cnt) + (node_offset[k]+data->conf[k]);
                   t0 += reduce_energy[cord0];
                   }
                 */

                int cord = (idx)*(tree_level+1) + i+1;
                cur_energy += pm_energy[cord];
                min_energy = min(min_energy, cur_energy);
            }
            h_score += min_energy;
        }
    }
    data->f_score = g_score + h_score;
    data->g_score = g_score;
}

// initialize everything
__global__
void d_initialize(
        int tree_level,
        int rot_cnt,

        heap_t *g_heap,
        data_t *g_data,
        int  g_heap_size[],
        int *g_data_size,

        uint *g_optimal,

        int *g_output_size,
        int *g_begin_index,
        int *g_begin_index2,

        int g_radix[],

        const int   g_node_offset[],
        const int   g_rot_per_level[],
        const float g_self_energy[],
        const float g_reduce_energy[],
        const float g_pm_energy[])
{
    int id  = blockDim.x*blockIdx.x + threadIdx.x;
    int lid = threadIdx.x;

    __shared__ int   node_offset[MAX_LEVEL];
    __shared__ int   rot_per_level[MAX_LEVEL];
    extern __shared__ float self_energy[];

    if (lid < tree_level) {
        node_offset[lid] = g_node_offset[lid];
        rot_per_level[lid] = g_rot_per_level[lid];
    }
    for (int i = id; i < rot_cnt; i += blockDim.x)
        self_energy[i] = g_self_energy[i];

    __syncthreads();

    g_heap_size[id] = 0;

    if (id < rot_per_level[0]) {
        data_t data;
        data.level = 0;
        data.conf[0] = id;
        __compute_score(tree_level,
                        rot_cnt,
                        &data,
                        0.f,
                        node_offset,
                        rot_per_level,
                        self_energy,
                        g_reduce_energy,
                        g_pm_energy);
        g_data[id] = data;

        int index = id*d_heap_capacity + 1;
        g_heap[index].v = data.f_score;
        g_heap[index].idx = id;
        g_heap_size[id] = 1;
        atomicMin(g_optimal, flip_float(data.f_score));

#ifdef VERBOSE_DEBUG
        printf("[%d]init: fscore: %.3f gscore: %.3f\n",
               id, data.f_score, data.g_score);
#endif
    }

    if (id < tree_level)
        g_radix[id] = 0;

    if (id == 0) {
        *g_data_size = rot_per_level[0];
        *g_output_size = *g_begin_index = *g_begin_index2 = 0;
    }
}

// extract the minimal element from heap
__global__ void d_delete_min(
        int tree_level,
        int num_child,

        heap_t g_heap[],
        data_t g_data[],
        data_t g_output[],

        int  g_heap_size[],
        int *g_data_size,
        int *g_output_size,

        int *g_begin_index,
        int *g_begin_index2,

        uint *g_optimal,

        int  g_parent[],
        int  g_radix[],
        int  g_node_cnt[],

        const int g_rot_per_level[])
{
    int id  = blockDim.x*blockIdx.x + threadIdx.x;
    int lid = threadIdx.x;

    int heap_size = g_heap_size[id];
    bool in_work = (heap_size != 0);

    /* =========================== delete_min ============================ */
    heap_t *heap = g_heap + id*d_heap_capacity;

    float node_score;
    int node_index = -1;
    int cur_level;

    if (in_work) {
        node_score = heap[1].v;
        node_index = heap[1].idx;
        cur_level  = g_data[node_index].level + 1;

#ifdef VERBOSE_DEBUG
        printf("[%d]pop node: %d:%.3f\n", id, node_index, heap[1].v);
#endif

        heap_t now_val = heap[heap_size--];
        g_heap_size[id] = heap_size;

        /* pop from heap */
        int now = 1;
        int next;
        while ((next = now*2) <= heap_size) {
            heap_t next_val = heap[next];
            heap_t next_val2 = heap[next+1];
            bool inc = (next+1 <= heap_size) && (next_val2.v < next_val.v);
            if (inc) {
                next += 1;
                next_val = next_val2;
            }

            if (next_val.v < now_val.v) {
                heap[now] = next_val;
                now = next;
            } else
                break;
        }
        heap[now] = now_val;
    }

    __shared__ int  rot_per_level[MAX_LEVEL];
    __shared__ int  l_radix[MAX_LEVEL];
    __shared__ uint l_optimal;
    l_optimal = UINT_MAX;
    if (lid < tree_level) {
        l_radix[lid] = 0;
        rot_per_level[lid] = g_rot_per_level[lid];
    }
    __syncthreads();

    /* =========================== update answer ============================ */
    if (in_work) {
        if (cur_level == tree_level) {
            /* add to answer array if we are leaves */
            int idx = atomicAdd(g_output_size, 1);
#ifdef VERBOSE_DEBUG
            printf("   >>   puts %d(%.8f, %.8f) on %d\n", node_index, node_score, g_data[node_index].f_score, idx);
#endif
            g_output[idx] = g_data[node_index];
            if (heap_size > 0)
                atomicMin(&l_optimal, flip_float(heap[1].v));

            in_work = false;
            node_index = -1;
        } else
            atomicMin(&l_optimal, flip_float(node_score));
    }
    __syncthreads();
    if (lid == 0) {
        atomicMin(g_optimal, l_optimal);
    }
    /* ======================== radix sort first part ======================== */
    g_parent[id] = node_index;

    if (in_work)
        atomicAdd(&l_radix[cur_level], rot_per_level[cur_level]);
    __syncthreads();
    if (lid < tree_level) {
        g_node_cnt[id] = l_radix[lid];
    }

    for (int i = 1; i < tree_level; i *= 2) {
        if (lid - i >= 0 && lid < tree_level)
            l_radix[lid] += l_radix[lid - i];
        __syncthreads();
    }
    if (lid < tree_level) {
        atomicAdd(&g_radix[lid], l_radix[lid]);
    }

    if (id == 0) {
        *g_data_size += num_child;
        *g_begin_index = *g_begin_index2;
    }
}

// sort extracted node by their level to decreasing branch divergance
__global__ void d_radix_sort(
        int tree_level,
        data_t g_data[],

        int *g_data_size,

        int g_radix[],
        int g_parent[],
        int g_node_cnt[],

        int2 g_input[],

        const int g_rot_per_level[])
{
    __shared__ int rot_per_level[MAX_LEVEL];
    __shared__ int l_radix[MAX_LEVEL];

    int id  = blockDim.x*blockIdx.x + threadIdx.x;
    int lid = threadIdx.x;

    if (lid < tree_level) {
        rot_per_level[lid] = g_rot_per_level[lid];

        int t = g_node_cnt[id];
        l_radix[lid] = atomicSub(&g_radix[lid], t) - t;
    }
    __syncthreads();

    int node_index = g_parent[id];
#ifdef VERBOSE_DEBUG
    printf("[%d]radix sort: node_index: %d\n", id, node_index);
#endif

    if (node_index >= 0) {
        int cur_level = g_data[node_index].level + 1;
        int index = atomicAdd(&l_radix[cur_level], rot_per_level[cur_level]);
        for (int i = 0; i < rot_per_level[cur_level]; ++i) {
            g_input[index++] = make_int2(node_index, i);
        }
    }
}

__global__ void d_compute_score(
        int tree_level,
        int rot_cnt,
        int num_child,

        data_t g_data[],
        int2 g_input[],

        int *g_data_size,

        const int   g_node_offset[],
        const int   g_rot_per_level[],
        const float g_self_energy[],
        const float g_reduce_energy[],
        const float g_pm_energy[])
{
    int id  = blockDim.x*blockIdx.x + threadIdx.x;
    int lid = threadIdx.x;

    __shared__ int   node_offset[MAX_LEVEL];
    __shared__ int   rot_per_level[MAX_LEVEL];
    extern __shared__ float self_energy[];

    if (lid < tree_level) {
        node_offset[lid] = g_node_offset[lid];
        rot_per_level[lid] = g_rot_per_level[lid];
    }
    for (int i = lid; i < rot_cnt; i += blockDim.x)
        self_energy[i] = g_self_energy[i];

    __syncthreads();

    if (id >= num_child)
        return;

    int2   input = g_input[id];
    data_t data  = g_data[input.x];
    float  score = data.g_score;

    data.level += 1;
    data.conf[data.level] = input.y;

    __compute_score(tree_level,
                    rot_cnt,
                    &data,
                    score,
                    node_offset,
                    rot_per_level,
                    self_energy,
                    g_reduce_energy,
                    g_pm_energy);

    int index = *g_data_size + id;
    g_data[index] = data;

#ifdef VERBOSE_DEBUG
    printf("[%d]Compute Score: (%d %d) to %d\n",
           id, input.x, input.y, index);
    printf("[%d]Compute Score: fscore: %.3f gscore: %.3f\n",
           id, data.f_score, data.g_score);
#endif
}

// put the generated node back to heap
__global__ void d_push_back(
        int tree_level,
        int num_child,

        heap_t g_heap[],
        data_t g_data[],
        int    g_heap_size[],
        int   *g_data_size,

        uint *g_optimal,

        int g_radix[],

        int *g_begin_index,
        int *g_begin_index2,
        int *g_output_size)
{
    int global_size = gridDim.x * blockDim.x;
    int id  = blockDim.x*blockIdx.x + threadIdx.x;
    int lid = threadIdx.x;

    int data_size  = *g_data_size;
    int data_size2 = data_size + num_child;

    int index = id - *g_begin_index;
    index = (index < 0 ? index + global_size : index);
    index += data_size;

    heap_t *heap = g_heap + d_heap_capacity*id;
    int heap_size = g_heap_size[id];

    uint optimal = (uint)-1;
    while (index < data_size2) {
        heap_t val;
        val.v = g_data[index].f_score;
        val.idx = index;

        optimal = min(optimal, flip_float(val.v));
#ifdef VERBOSE_DEBUG
        printf("[%d]: assign node (%.3f, %d) to this heap\n", id, val.v, val.idx);
#endif

        int now = ++heap_size;
        while (now > 1) {
            int next = now / 2;
            heap_t next_val = heap[next];
            if (val.v < next_val.v) {
                heap[now] = next_val;
                now = next;
            } else
                break;
        }
        heap[now] = val;

        index += global_size;
    }
    g_heap_size[id] = heap_size;

    __shared__ uint l_optimal;
    l_optimal = (uint)-1;
    __syncthreads();
    atomicMin(&l_optimal, optimal);
    __syncthreads();
    if (lid == 0)
        atomicMin(g_optimal, l_optimal);

    if (index == data_size2)
        *g_begin_index2 = id;
    if (id < tree_level)
        g_radix[id] = 0;
    if (id == 0)
        *g_output_size = 0;
}

// Shrink operation: delete node which unlikely become the answer {{{
// This can speed up calculation and do GMAC in a bounded memroy
// But it may not give the optimal answer if we apply shrink operation
// See also http://en.wikipedia.org/wiki/SMA* for detailed explaination.
__global__ void d_tagging(
        heap_t g_heap[],
        int    g_heap_size[],
        int    g_data_used[])
{
    int id  = blockDim.x*blockIdx.x + threadIdx.x;

    heap_t *heap = g_heap + d_heap_capacity*id;
    int heap_size = g_heap_size[id];

    for (int i = 1; i <= heap_size; ++i)
        g_data_used[heap[i].idx] = 1;
}

__global__ void d_scatter_data(
        int    data_size,
        int    g_data_used[],
        data_t g_old_data[],
        data_t g_new_data[],
        float  g_data_val[])
{
    int id  = blockDim.x*blockIdx.x + threadIdx.x;

    // When id == 0, we assume that it should be free.
    if (id == 0 || id >= data_size)
        return;

    if (id > 0) {
        int pos = g_data_used[id-1];
        if (g_data_used[id] != pos) {
            data_t data = g_old_data[id];
            g_new_data[pos] = data;
            g_data_val[pos] = data.f_score;
        }
    }
}

__global__ void d_reassign(
        heap_t g_heap[],
        data_t g_data[],
        int    g_heap_size[],
        int   *g_data_size)
{
    int id  = blockDim.x*blockIdx.x + threadIdx.x;
    int global_size = gridDim.x * blockDim.x;

    heap_t *heap = g_heap + d_heap_capacity*id;
    int data_size = *g_data_size;

    int cnt = 0;
    for (int i = id; i < data_size; i += global_size) {
        ++cnt;
        heap[cnt].v = g_data[i].f_score;
        heap[cnt].idx = i;
    }
    g_heap_size[id] = cnt;
}

void __print_data(int data_size)
{
    puts(">>> Print heap data <<<");

    printf("data size: %d\n", data_size);
    for (int i = 0; i < data_size; ++i) {
        data_t t = copy_from_device(d_data+i);
        printf("[%d]fscore: %.3f, level: %d\n", i, t.f_score, t.level);
    }

    puts(">>> End heap data <<<");
}

void shrink(int *__data_size, int *__num_child)
{
    int data_size = *__data_size + *__num_child;
    if (data_size + data_reserved > data_capacity || data_size > max_data) {
        if (shrink_ratio == 1) {
            puts("  >>> GPU run out of memory! Running Failed");
            exit(EXIT_FAILURE);
        }
        puts("  >>> GPU run out of memory! Rescan to free memory");
        *__num_child = 0;

        if (!d_data_val)
            CHK_CUDA(cudaMalloc(&d_data_val, data_capacity * sizeof(float)));
        if (!d_data_used)
            CHK_CUDA(cudaMalloc(&d_data_used, data_capacity * sizeof(int)));
        if (!d_data2)
            CHK_CUDA(cudaMalloc(&d_data2, data_capacity * sizeof(data_t)));

        puts("OK0");fflush(stdout);

        CHK_CUDA(cudaMemset(d_data_used, 0, data_size * sizeof(int)));
        d_tagging<<<num_block, num_local>>>(
                d_heap,
                d_heap_size,
                d_data_used);

        puts("OK1");fflush(stdout);

        thrust::device_ptr<int> dev_used(d_data_used);
        thrust::inclusive_scan(dev_used, dev_used+data_size, dev_used);
        int new_data_size = dev_used[data_size-1];

        printf("data_size: %d\nnew_data_size: %d\n", data_size, new_data_size);
        fflush(stdout);

        d_scatter_data<<<(data_size-1) / num_local + 1, num_local>>>(
                data_size,
                d_data_used,
                d_data,
                d_data2,
                d_data_val);

        std::swap<data_t *>(d_data, d_data2);

        size_t free0;
        size_t total;
        cuMemGetInfo(&free0, &total);
        printf("free: %zuMB, total: %zuMB\n", free0/1024/1024, total/1024/1024);

        thrust::device_ptr<data_t> dev_data(d_data);
        thrust::device_ptr<float>  dev_val(d_data_val);
        thrust::sort_by_key(dev_val, dev_val+new_data_size, dev_data);
        // thrust::host_vector<data_t> h_data(dev_data, dev_data + new_data_size);
        // thrust::sort(h_data.begin(), h_data.end());
        // thrust::copy(h_data.begin(), h_data.end(), dev_data);


        new_data_size *= shrink_ratio;

        float curr_throw_min = dev_val[new_data_size];
        throw_min = min(throw_min, curr_throw_min);

        copy_to_device(d_data_size, new_data_size);
        *__data_size = new_data_size;
        printf("  >>> %d elements left\n", new_data_size);

        d_reassign<<<num_block, num_local>>>(
                d_heap,
                d_data,
                d_heap_size,
                d_data_size);
    }
}
// }}}

extern "C" void init_cuda()
{
    data_capacity = (gpu_memory - 100*1024*1024) / (sizeof(data_t)*(shrink_ratio == 1.f ? 1 : 3) + sizeof(heap_t));
    data_capacity = int(data_capacity / 1.2f);

    printf("data_capacity: %d\n", data_capacity);
    data_reserved = data_capacity / 12;
    printf("data_reserved: %d\n", data_reserved);
    num_block     = num_gpu_group;
    printf("num_block: %d\n", num_block);
    num_local     = num_gpu_item;
    printf("num_local: %d\n", num_local);
    num_local2    = num_gpu_item2;
    printf("num_local2: %d\n", num_local2);
    num_global    = num_block * num_local;
    printf("num_global: %d\n", num_global);
    heap_capacity = data_capacity / num_global;
    printf("heap_capacity: %d\n", heap_capacity);
    cudaMemcpyToSymbol(d_heap_capacity, &heap_capacity, sizeof(heap_capacity));

    int max_child = 0;
    func_reduce(max_child, rot_per_level, rot_per_level+tree_level, max);
    if (max_child >= 255) {
        printf("CUDA K* do not support more than 255 rotemar in single\n"
               "residue due to memory limit.\n\n"
               "You may want to change some data type char to short in\n"
               "source jni/astar-cuda.cu to avoid this limitation.\n"
               "But that would almost double the memory usage\n");
        exit(EXIT_FAILURE);
    }
    if (tree_level > MAX_LEVEL) {
        printf("CUDA K* do not support more than %d residues due to memory"
               "bound limit.\n\n"
               "You may want to change constant in jni/astar.h to avoid this"
               "limitation.  But that would almost double the memory usage\n",
               MAX_LEVEL);
        exit(EXIT_FAILURE);
    }

    /*cuda_free(d_heap);*/
    /*cuda_free(d_heap);*/
    /*cuda_free(d_data);*/
    /*cuda_free(d_heap_size);*/
    /*cuda_free(d_output);*/
    /*cuda_free(d_parent);*/
    /*cuda_free(d_node_cnt);*/
    /*cuda_free(d_radix);*/
    /*cuda_free(d_input);*/
    /*cuda_free(d_data_size);*/
    /*cuda_free(d_optimal);*/
    /*cuda_free(d_output_size);*/
    /*cuda_free(d_begin_index);*/
    /*cuda_free(d_begin_index);*/
    /*cuda_free(d_node_offset);*/
    /*cuda_free(d_rot_per_level);*/
    /*cuda_free(d_self_energy);*/
    /*cuda_free(d_reduce_energy);*/
    /*cuda_free(d_pm_energy);*/

    /*cuda_free(d_data_val);*/
    /*cuda_free(d_data_used);*/
    /*cuda_free(d_data2);*/

    if (!d_heap)
        CHK_CUDA(cudaMalloc(&d_heap, heap_capacity * num_global * sizeof(heap_t)));
    if (!d_data)
        CHK_CUDA(cudaMalloc(&d_data, data_capacity * sizeof(data_t)));
    if (!d_heap_size)
        CHK_CUDA(cudaMalloc(&d_heap_size, num_global * sizeof(int)));
    if (!d_output)
        CHK_CUDA(cudaMalloc(&d_output, num_global * sizeof(data_t)));
    if (!d_parent)
        CHK_CUDA(cudaMalloc(&d_parent, num_global * sizeof(int)));
    if (!d_node_cnt)
        CHK_CUDA(cudaMalloc(&d_node_cnt, num_global * sizeof(int)));
    if (!d_radix)
        CHK_CUDA(cudaMalloc(&d_radix, tree_level * sizeof(int)));
    if (!d_input)
        CHK_CUDA(cudaMalloc(&d_input, num_global * MAX_ROTAMER * sizeof(int2)));

    if (!d_data_size)
        CHK_CUDA(cudaMalloc(&d_data_size, sizeof(int)));
    if (!d_optimal)
        CHK_CUDA(cudaMalloc(&d_optimal, sizeof(int)));
    if (!d_output_size)
        CHK_CUDA(cudaMalloc(&d_output_size, sizeof(int)));
    if (!d_begin_index)
        CHK_CUDA(cudaMalloc(&d_begin_index, sizeof(int)));
    if (!d_begin_index2)
        CHK_CUDA(cudaMalloc(&d_begin_index2, sizeof(int)));

    if (!d_node_offset)
        CHK_CUDA(cudaMalloc(&d_node_offset, tree_level * sizeof(int)));
    if (!d_rot_per_level)
        CHK_CUDA(cudaMalloc(&d_rot_per_level, tree_level * sizeof(int)));
    if (!d_self_energy)
        CHK_CUDA(cudaMalloc(&d_self_energy, MAX_ROTAMER * sizeof(float)));
    if (!d_reduce_energy)
        CHK_CUDA(cudaMalloc(&d_reduce_energy, MAX_ROTAMER * MAX_ROTAMER * sizeof(float)));
    if (!d_pm_energy)
        CHK_CUDA(cudaMalloc(&d_pm_energy, MAX_ROTAMER * (tree_level+1) * sizeof(float)));

    static float h_self_energy[MAX_ROTAMER];
    static float h_reduce_energy[MAX_ROTAMER*MAX_ROTAMER];
    static float h_pm_energy[MAX_ROTAMER*(MAX_LEVEL+1)];

    for (int i = 0; i < rot_cnt; ++i)
        h_self_energy[i] = reduce_energy[rot_cnt][i] + reduce_energy[i][rot_cnt];
    float *ptr;
    ptr = h_reduce_energy;
    for (int i = 0; i < rot_cnt; ++i)
        for (int j = 0; j < rot_cnt; ++j)
            *ptr++ = reduce_energy[i][j];
    ptr = h_pm_energy;
    for (int i = 0; i < rot_cnt; ++i)
        for (int j = 0; j <= tree_level; ++j)
            *ptr++ = pm_energy[i][j];

    printf("node offset:\n");
    for (int i = 0; i < tree_level; ++i)
        printf("%d ", node_offset[i]);
    printf("\nrot_per_level:\n");
    for (int i = 0; i < tree_level; ++i)
        printf("%d ", rot_per_level[i]);
    printf("\nself_energy:\n");
    /*
       for (int i = 0; i < rot_cnt; ++i)
       printf("%.3f ", h_self_energy[i]);
       printf("\n");
     */

    CHK_CUDA(cudaMemcpy(d_node_offset,
                        node_offset,
                        tree_level * sizeof(int),
                        cudaMemcpyHostToDevice));
    CHK_CUDA(cudaMemcpy(d_rot_per_level,
                        rot_per_level,
                        tree_level * sizeof(int),
                        cudaMemcpyHostToDevice));
    CHK_CUDA(cudaMemcpy(d_self_energy,
                        h_self_energy,
                        rot_cnt * sizeof(float),
                        cudaMemcpyHostToDevice));
    CHK_CUDA(cudaMemcpy(d_reduce_energy,
                        h_reduce_energy,
                        rot_cnt * rot_cnt * sizeof(float),
                        cudaMemcpyHostToDevice));
    CHK_CUDA(cudaMemcpy(d_pm_energy,
                        h_pm_energy,
                        rot_cnt * (tree_level+1) * sizeof(float),
                        cudaMemcpyHostToDevice));

    rounds = 0;
    optimal = 0;
    h_num_child = 0;
    h_output_size = 0;
    h_data_size = 0;

    pq_for(node, heap) {
        /*printf("%p\n", node);*/
        /*fflush(stdout);*/
        free(*node);
    }
    pq_init(heap);

    max_data = INT32_MAX;
    char *quota_s = getenv("KSTAR_MAX_NODES");
    if (quota_s) {
        max_data = atoi(quota_s);
        printf("Environment KSTAR_MAX_NODES is setting to %d\n", max_data);
    }

    throw_min = +INFINITY;

    printf("CUDA init finish!\n");
}

extern "C" int *astar_cuda(bool first_run)
{
    puts("\n====== GPU A* start ======");
    fflush(stdout);

    cudaProfilerStart();
    if (first_run) {
        CHK_CUDA(cudaMemset(d_optimal, -1, sizeof(uint)));
        d_initialize<<<num_block, num_local, rot_cnt*sizeof(float)>>>(
                tree_level,
                rot_cnt,
                d_heap,
                d_data,
                d_heap_size,
                d_data_size,
                d_optimal,
                d_output_size,
                d_begin_index,
                d_begin_index2,
                d_radix,
                d_node_offset,
                d_rot_per_level,
                d_self_energy,
                d_reduce_energy,
                d_pm_energy);
        cudaDeviceSynchronize();
    }


    static data_t *h_data = NULL;
    if (h_data == NULL)
        CHK_CUDA(cudaMallocHost(&h_data, num_global * sizeof(data_t)));

#define check_return() \
    if (!pq_empty(heap) && flip_float(pq_top(heap)->f_score) <= optimal) { \
        printf("GPU best result: %.9f\n", pq_top(heap)->f_score); \
        if (pq_top(heap)->f_score > throw_min) \
            printf("!!! fscore is greater than the minimal throw element.\n" \
                   "!!! GMEC is not guaranteed:(\n"); \
        for (int i = 0; i < tree_level; ++i) { \
            curr_conf[i] = pq_top(heap)->conf[i]; \
            printf("%d ", curr_conf[i]); \
        } \
        printf("\n"); \
        free(pq_pop(heap, heap_cmp)); \
        printf("GPU native A* finish in %d ms and %d rounds\n", \
               (int)wall_time_elapsed(), rounds); \
        printf("GPU memory: %d out of %d\n", h_data_size, data_capacity); \
        puts("====== GPU A* finished ======"); \
        fflush(stdout); \
        cudaProfilerStop(); \
        return curr_conf; \
    } else

    wall_time_begin();
    check_return();

    for (;;) {
        rounds++;
        // printf("Current round: %d\n", rounds);

        d_delete_min<<<num_block, num_local>>>(
                tree_level,
                h_num_child,
                d_heap,
                d_data,
                d_output,
                d_heap_size,
                d_data_size,
                d_output_size,
                d_begin_index,
                d_begin_index2,
                d_optimal,
                d_parent,
                d_radix,
                d_node_cnt,
                d_rot_per_level);

        h_num_child   = copy_from_device(d_radix+tree_level-1);
        h_output_size = copy_from_device(d_output_size);
        // printf("num_child: %d\n", h_num_child);

        // No answer
        if (!first_run && h_num_child == 0) {
            memset(curr_conf, -1, sizeof(int)*(uint)tree_level);
            printf("GPU native A* finish in %d ms and %d rounds\n",
                   (int)wall_time_elapsed(), rounds);
            printf("GPU memory: %d out of %d\n", h_data_size, data_capacity);
            puts("====== GPU A* finished ======");
            return curr_conf;
        }
        if (h_output_size > 0) {
            CHK_CUDA(cudaMemcpy(h_data,
                                d_output,
                                h_output_size * sizeof(data_t),
                                cudaMemcpyDeviceToHost));
            for (int i = 0; i < h_output_size; ++i) {
                data_t *node = (data_t *)malloc(sizeof(data_t));
                *node = h_data[i];
#ifdef VERBOSE_DEBUG
                printf("I reveice %.8f\n", node->f_score);
#endif
                pq_push(heap, node, heap_cmp);
            }
        }

        d_radix_sort<<<num_block, num_local>>>(
                tree_level,
                d_data,
                d_data_size,
                d_radix,
                d_parent,
                d_node_cnt,
                d_input,
                d_rot_per_level);
        int num_block2 = (h_num_child-1) / num_local2 + 1;
        d_compute_score<<<num_block2, num_local2, rot_cnt*sizeof(float)>>>(
                tree_level,
                rot_cnt,
                h_num_child,
                d_data,
                d_input,
                d_data_size,
                d_node_offset,
                d_rot_per_level,
                d_self_energy,
                d_reduce_energy,
                d_pm_energy);
        d_push_back<<<num_block, num_local>>>(
                tree_level,
                h_num_child,
                d_heap,
                d_data,
                d_heap_size,
                d_data_size,
                d_optimal,
                d_radix,
                d_begin_index,
                d_begin_index2,
                d_output_size);

        optimal = copy_from_device(d_optimal);
        CHK_CUDA(cudaMemset(d_optimal, -1, sizeof(uint)));
        // printf("optimal: %.9f\n", reverse_flip_float(optimal));

        h_data_size = copy_from_device(d_data_size);
        shrink(&h_data_size, &h_num_child);
        // printf("h_data_size: %d\n", h_data_size);
        // printf("h_num_child: %d\n", h_num_child);

        check_return();
    }
}
