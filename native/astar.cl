#ifndef NUM_LEVEL
# define NUM_LEVEL 24
#endif
#ifndef HEAP_SIZE
# define HEAP_SIZE 997
#endif
#ifndef NUM_ROTAMER
# define NUM_ROTAMER 100
#endif
#ifndef NUM_MAX_CHILD
# define NUM_MAX_CHILD 35
#endif
#ifndef GROUP_SIZE
# define GROUP_SIZE 256
#endif

#define nVERBOSE_DEBUG

#define SMP (CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_NONE | CLK_FILTER_NEAREST)
#define PR  (__constant char *)

typedef struct data_t {
    float f_score;
    float g_score;
    int level;
    uchar conf[NUM_LEVEL];
} data_t;

typedef struct heap_t {
    float v;
    int idx;
} heap_t;

uint flip_float(float f)
{
    int u = as_int(f);
    return u ^ ((u >> 31) | 0x80000000);
}

void print_node(int level, float f_score, float g_score, int conf[])
{
    printf(PR"level: %d\t f_score: %.4f g_score: %.4f\n",
           level, f_score, g_score);
    for (int i = 0; i <= level; ++i)
        printf(PR"%d ", conf[i]);
    printf(PR"\n");
}

/* return a pair of int (fscore, gscore) */
float2 __compute_score(int *conf,
                       int level,
                       float old_g_score,
                       const __local ushort node_offset[],
                       const __local uchar  rot_per_level[],
                       const __local float  self_energy[],
                       __read_only image2d_t reduce_energy,
                       __read_only image2d_t pm_energy)
{
    sampler_t sampler = SMP;
    float g_score = old_g_score; /* use previous g score */
    float h_score = 0.f;

    /* compute g delta */
    {
        int idx = node_offset[level] + conf[level];
        for (int j = 0; j < level; ++j) {
            int2 cord = (int2)(idx, node_offset[j] + conf[j]);
            /*cord = select(cord, cord.yx, cord.x < cord.y); */
            g_score += read_imagef(reduce_energy, sampler, cord).x;
        }
        g_score += self_energy[idx];
    }

    /* compute h */
    {
        for (int i = level+1; i < NUM_LEVEL; ++i) {
            float min_energy = FLT_MAX;
            for (int j = 0; j < rot_per_level[i]; ++j) {
                int idx = node_offset[i]+j;
                float cur_energy = self_energy[idx];

                for (int k = 0; k <= level; ++k) {
                    int2 cord = (int2)(idx, node_offset[k] + conf[k]);
                    /*cord = select(cord, cord.yx, cord.x < cord.y);*/
                    cur_energy += read_imagef(reduce_energy, sampler, cord).x;
                }

                cur_energy += read_imagef(pm_energy, sampler, (int2)(i+1, idx)).x;
                min_energy = min(min_energy, cur_energy);
            }
            h_score += min_energy;
        }
    }

    return (float2)(g_score + h_score, g_score);
}

__kernel
void initialize(
        __global heap_t *g_heap,
        __global int *g_heap_size,
        __global data_t *g_data,
        __global int *g_data_size,     /* [0 .. data_size-1] was filled */

        __global uint *g_optimal,  /* this is a reinterperation of float */

        const __global ushort *g_node_offset,
        const __global uchar  *g_rot_per_level,
        const __global float  *g_self_energy,
        __read_only image2d_t reduce_energy,
        __read_only image2d_t pm_energy)
{
    int id = get_global_id(0);
#ifdef VERBOSE_DEBUG
    printf(PR"\n\n================initialize[%d]===================\n", id);
#endif

    __local ushort node_offset[NUM_LEVEL];
    __local uchar  rot_per_level[NUM_LEVEL];
    __local float  self_energy[NUM_ROTAMER];

    event_t event[3];
    event[0] = async_work_group_copy(node_offset,   g_node_offset,   NUM_LEVEL,   0);
    event[1] = async_work_group_copy(rot_per_level, g_rot_per_level, NUM_LEVEL,   0);
    event[2] = async_work_group_copy(self_energy,   g_self_energy,   NUM_ROTAMER, 0);

#ifdef VERBOSE_DEBUG
    printf(PR"size of data_t: %d\n", sizeof(data_t));
#endif

    int conf[NUM_LEVEL];
    float2 score;

    conf[0] = id;

    wait_group_events(3, event);
    score = __compute_score(conf,
                            0,
                            0.f,
                            node_offset,
                            rot_per_level,
                            self_energy,
                            reduce_energy,
                            pm_energy);
    g_heap[id*HEAP_SIZE + 1] = ((heap_t){ score.x, id });
    g_heap_size[id] = 1;
    atomic_min(g_optimal, flip_float(score.x));

    g_data[id].f_score = score.x;
    g_data[id].g_score = score.y;
    g_data[id].level   = 0;
    g_data[id].conf[0] = id;

#ifdef VERBOSE_DEBUG
    print_node(0, score.x, score.y, conf);
#endif
}

__kernel __attribute__((reqd_work_group_size(GROUP_SIZE, 1, 1)))
void delete_min(
        __global heap_t *g_heap,
        __global int    *g_heap_size,
        __global data_t *g_data,
        __global int    *g_data_size,
        __global uint   *g_optimal,  /* uint is a reinterperation of float */

        __global int    *g_father,
        __global int    *g_radix,
        __global int    *g_node_cnt,
        __global data_t *g_answer,
        __global int    *g_answer_size,

        const __global uchar  *g_rot_per_level)
{
    __local uchar rot_per_level[NUM_LEVEL];

    int id  = get_global_id(0);
    int lid = get_local_id(0);

#ifdef VERBOSE_DEBUG
    printf(PR"\n\n================delete_min[%d]=================\n", id);
#endif

    int heap_size = g_heap_size[id];
    bool in_work = (heap_size != 0);

    /* =========================== delete_min ============================ */
    __global heap_t *heap = g_heap + id * HEAP_SIZE;

    float node_score;
    int node_index = -1;
    int cur_level;

    if (in_work) {
        node_score = heap[1].v;
        node_index = heap[1].idx;
        cur_level  = g_data[node_index].level + 1;

#ifdef VERBOSE_DEBUG
        printf(PR"pop node: %d:%.3f\n", node_index, heap[1].v);
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

    __local int  l_radix[NUM_LEVEL];
    __local uint l_optimal;
    l_optimal = UINT_MAX;
    if (lid < NUM_LEVEL) {
        l_radix[lid] = 0;
        rot_per_level[lid] = g_rot_per_level[lid];
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    /* =========================== update answer ============================ */
    if (in_work) {
        if (cur_level == NUM_LEVEL) {
            /* add to answer array if we are leaves */
            int answer_index = atomic_inc(g_answer_size);
            g_answer[answer_index] = g_data[node_index];
            if (heap_size > 0)
                atomic_min(&l_optimal, flip_float(heap[1].v));

            in_work = false;
            node_index = -1;
        } else
            atomic_min(&l_optimal, flip_float(node_score));
    }
    barrier(CLK_LOCAL_MEM_FENCE);
    if (lid == 0) {
        atomic_min(g_optimal, l_optimal);
    }
    /* ======================== radix sort first part ======================== */
    g_father[id] = node_index;

    if (in_work)
        atomic_add(&l_radix[cur_level], rot_per_level[cur_level]);
    barrier(CLK_LOCAL_MEM_FENCE);
    if (lid < NUM_LEVEL) {
        g_node_cnt[id] = l_radix[lid];
    }

#pragma unroll
    for (int i = 1; i < NUM_LEVEL; i *= 2) {
        if (lid - i >= 0 && lid < NUM_LEVEL)
            l_radix[lid] += l_radix[lid - i];
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    if (lid < NUM_LEVEL) {
        atomic_add(&g_radix[lid], l_radix[lid]);
    }
}

__kernel __attribute__((reqd_work_group_size(GROUP_SIZE, 1, 1)))
void radix_sort(__global int  *g_father,
                __global int2 *g_input,
                __global int  *g_radix,
                __global int  *g_node_cnt,

                __global data_t *g_data,

                const __global uchar *g_rot_per_level)
{
    __local uchar rot_per_level[NUM_LEVEL];
    __local int   l_radix[NUM_LEVEL];

    int id  = get_global_id(0);
    int lid = get_local_id(0);

#ifdef VERBOSE_DEBUG
    printf(PR"\n\n================radix_sort[%d]=================\n", id);
#endif

    if (lid < NUM_LEVEL) {
        rot_per_level[lid] = g_rot_per_level[lid];

        int t = g_node_cnt[id];
        l_radix[lid] = atomic_sub(&g_radix[lid], t) - t;
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    int node_index = g_father[id];
#ifdef VERBOSE_DEBUG
    printf(PR"node_index: %d\n", node_index);
#endif
    if (node_index >= 0) {
        int cur_level = g_data[node_index].level + 1;
        int index = atomic_add(&l_radix[cur_level], rot_per_level[cur_level]);
#ifdef VERBOSE_DEBUG
        printf(PR"cur_level: %d\n", cur_level);
#endif
        for (int i = 0; i < rot_per_level[cur_level]; ++i) {
            g_input[index++] = (int2)(node_index, i);
#ifdef VERBOSE_DEBUG
            printf(PR"g_input[%d] = (index: %d, conf: %d)\n", index-1, node_index, i);
#endif
        }
    }
}



__kernel __attribute__((reqd_work_group_size(GROUP_SIZE, 1, 1)))
void compute_score(
        __global int2   *g_input,
        __global data_t *g_data,
        __global int    *g_data_size,
        __global int    *g_output_size,

        const __global ushort *g_node_offset,
        const __global uchar  *g_rot_per_level,
        const __global float  *g_self_energy,
        __read_only image2d_t  reduce_energy,
        __read_only image2d_t  pm_energy)
{
    __local uchar  rot_per_level[NUM_LEVEL];
    __local ushort node_offset[NUM_LEVEL];
    __local float  self_energy[NUM_ROTAMER];

    event_t event[3];
    event[0] = async_work_group_copy(rot_per_level, g_rot_per_level, NUM_LEVEL,   0);
    event[1] = async_work_group_copy(node_offset,   g_node_offset,   NUM_LEVEL,   0);
    event[2] = async_work_group_copy(self_energy,   g_self_energy,   NUM_ROTAMER, 0);

    int id = get_global_id(0);
    int global_size = get_global_size(0);

#ifdef VERBOSE_DEBUG
    printf(PR"\n\n================compute_score[%d]=================\n", id);
#endif

    int index  = *g_data_size + id;

    if (index >= *g_output_size)
        return;

    int2   input;
    int    level;
    float  oscore;
    float2 score;
    int    conf[NUM_LEVEL];
    
    input  = g_input[id];
    level  = g_data[input.x].level + 1;
    oscore = g_data[input.x].g_score;
#pragma unroll
    for (int i = 0; i < NUM_LEVEL; ++i)
        conf[i] = g_data[input.x].conf[i];
    conf[level] = input.y;

    wait_group_events(3, event);
    score = __compute_score(conf,
                            level,
                            oscore,
                            node_offset,
                            rot_per_level,
                            self_energy,
                            reduce_energy,
                            pm_energy);

    g_data[index].f_score = score.x;
    g_data[index].g_score = score.y;
    g_data[index].level   = level;
#pragma unroll
    for (int i = 0; i < NUM_LEVEL; ++i)
        g_data[index].conf[i] = conf[i];

#ifdef VERBOSE_DEBUG
    printf(PR"index[%d]: ", index);
    print_node(level, score.x, score.y, conf);
#endif
}


__kernel __attribute__((reqd_work_group_size(GROUP_SIZE, 1, 1)))
void push_back(
        __global heap_t *g_heap,
        __global int *g_heap_size,
        __global data_t *g_data,
        __global int *g_data_size,
        __global int *g_output_size,
        __global int *begin_index,
        __global int *begin_index2,
        __global uint *g_optimal)
{
    int global_size = get_global_size(0);
    int id  = get_global_id(0);
    int lid = get_local_id(0);

#ifdef VERBOSE_DEBUG
    printf(PR"\n\n================push_back[%d]=================\n", id);
    if (id == 0)
        printf(PR"g_data_size: %d g:output: %d begin: %d\n",
               *g_data_size, *g_output_size, *begin_index);
#endif

    int data_size = *g_data_size;
    int output_size = *g_output_size;

    int index = id - *begin_index;
    index = select(index, index + global_size, index < 0);
    index += data_size;

    __global heap_t *heap = g_heap + HEAP_SIZE*id;
    int heap_size = g_heap_size[id];

    uint optimal = -1;
    while (index < output_size) {
        heap_t val = { g_data[index].f_score, index };
        optimal = min(optimal, flip_float(val.v));
#ifdef VERBOSE_DEBUG
        printf(PR"assign node [%d] to this heap\n", index);
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
    if (index == output_size)
        *begin_index2 = id;

    __local uint l_optimal;
    l_optimal = -1;
    barrier(CLK_LOCAL_MEM_FENCE);
    atomic_min(&l_optimal, optimal);
    barrier(CLK_LOCAL_MEM_FENCE);
    if (lid == 0)
        atomic_min(g_optimal, l_optimal);
}
