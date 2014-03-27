#include "astar-cl.h"
#include "astar-opencl.h"
#include "clwrapper.h"
#include "pq.h"

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

static size_t num_local_item = 64;
static size_t num_global_item;

#define DATA_SIZE 10000000
#define HEAP_SIZE (DATA_SIZE / num_global_item / 5)

static cl_device_id     device;
static cl_context       context;
static cl_command_queue cmd_queue;
static cl_program       program;
static cl_kernel        kernel_initialize;
static cl_kernel        kernel_delete_min;
static cl_kernel        kernel_radix_sort;
static cl_kernel        kernel_compute_score;
static cl_kernel        kernel_push_back;

static size_t           data_t_size;

static cl_mem d_heap;          /* heap_t heap[]   */
static cl_mem d_heap_size;     /* int heap_size   */
static cl_mem d_data;          /* data_t data[]   */
static cl_mem d_data_size;     /* int data_size   */
static cl_mem d_output_size;   /* int output_size */
static cl_mem d_optimal;       /* float optimal   */
static cl_mem d_answer;        /* data_t answer[] */
static cl_mem d_answer_size;   /* int answer_size */
static cl_mem d_father;        /* int father[] */
static cl_mem d_input;         /* int input[] */
static cl_mem d_radix;         /* int radix[] */
static cl_mem d_node_cnt;      /* int node_cnt[] */
static cl_mem d_begin_index;
static cl_mem d_begin_index2;
static cl_mem d_node_offset;
static cl_mem d_rot_per_level;
static cl_mem d_self_energy;
static cl_mem d_reduce_energy; /* image type */
static cl_mem d_pm_energy;     /* image type */

#define heap_cmp(a, b) ((a)->f_score < (b)->f_score)
static priority_queue(data_t *) heap;

uint flip_float(float f)
{
    union {
        float f;
        int  u;
    } un;
    un.f = f;
    return un.u ^ ((un.u >> 31) | 0x80000000);
}

float reverse_flip_float(int u)
{
    union {
        float f;
        int u;
    } un;
    un.u = u ^ ((~u >> 31) | 0x80000000);
    return un.f;
}

void init_gpu(void)
{
    device    = cl_get_first_gpu_device();
    context   = cl_create_context(device);
    cmd_queue = cl_create_cmd_queue(
            context,
            device,
            CL_QUEUE_PROFILING_ENABLE | CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE
    );

    cl_uint num_compute_unit;

    cl_query_device(device, CL_DEVICE_MAX_COMPUTE_UNITS,
                    sizeof(num_compute_unit), &num_compute_unit, NULL);
    printf("Num of CU: %d\n", num_compute_unit);
    /* num_global_item = num_local_item; */
    num_global_item = num_local_item * num_compute_unit;
    printf("Num of Work-Item: %d\n", (int)num_global_item);

    data_t_size = 12 + sizeof(uchar)*(unsigned)tree_level;
    if (data_t_size % 4 != 0)
        data_t_size += 4 - data_t_size % 4;
    printf("data_t_size: %d\n", (int)data_t_size);

    int max_child = 0;
    func_reduce(max_child, rot_per_level, rot_per_level+tree_level, max);

    d_heap    = cl_create_buffer(context,
                                CL_MEM_READ_WRITE,
                                HEAP_SIZE * num_global_item * sizeof(heap_t),
                                NULL);
    d_data    = cl_create_buffer(context,
                                CL_MEM_READ_WRITE,
                                DATA_SIZE * data_t_size,
                                NULL);
    d_answer  = cl_create_buffer(context,
                                CL_MEM_READ_WRITE,
                                num_global_item * data_t_size,
                                NULL);
    d_father  = cl_create_buffer(context,
                                CL_MEM_READ_WRITE,
                                num_global_item * sizeof(int),
                                NULL);
    d_input   = cl_create_buffer(context,
                                CL_MEM_READ_WRITE,
                                num_global_item * max_child * sizeof(int) * 2,
                                NULL);
    d_radix   = cl_create_buffer(context,
                                CL_MEM_READ_WRITE,
                                tree_level * sizeof(int),
                                NULL);
    d_node_cnt= cl_create_buffer(context,
                                CL_MEM_READ_WRITE,
                                num_global_item * sizeof(int),
                                NULL);

    d_heap_size    = cl_create_buffer(context, CL_MEM_READ_WRITE,
                                     sizeof(int) * num_global_item, NULL);
    d_data_size    = cl_create_buffer(context, CL_MEM_READ_WRITE, sizeof(int), NULL);
    d_answer_size  = cl_create_buffer(context, CL_MEM_READ_WRITE, sizeof(int), NULL);
    d_output_size  = cl_create_buffer(context, CL_MEM_READ_WRITE, sizeof(int), NULL);
    d_begin_index  = cl_create_buffer(context, CL_MEM_READ_WRITE, sizeof(int), NULL);
    d_begin_index2 = cl_create_buffer(context, CL_MEM_READ_WRITE, sizeof(int), NULL);
    d_optimal      = cl_create_buffer(context, CL_MEM_READ_WRITE, sizeof(float), NULL);

    int buf = 0;
    cl_write_buffer(cmd_queue, d_answer_size, CL_TRUE, 0, sizeof(int), &buf, 0, NULL, NULL);
    cl_write_buffer(cmd_queue, d_begin_index, CL_TRUE, 0, sizeof(int), &buf, 0, NULL, NULL);

    buf = -1;
    cl_write_buffer(cmd_queue, d_optimal    , CL_TRUE, 0, sizeof(int), &buf, 0, NULL, NULL);

    buf = rot_per_level[0];
    cl_write_buffer(cmd_queue, d_data_size  , CL_TRUE, 0, sizeof(int), &buf, 0, NULL, NULL);
    cl_write_buffer(cmd_queue, d_output_size, CL_TRUE, 0, sizeof(int), &buf, 0, NULL, NULL);

    void *ptr = calloc(num_global_item, sizeof(int));
    cl_write_buffer(cmd_queue, d_heap_size, CL_TRUE, 0, sizeof(int)*num_global_item, ptr, 0, NULL, NULL);
    free(ptr);

    static ushort m_node_offset[MAX_LEVEL];
    static uchar  m_rot_per_level[MAX_LEVEL];
    static float  m_self_energy[MAX_ROTAMER];
    static float  m_reduce_energy[MAX_ROTAMER*MAX_ROTAMER];
    static float  m_pm_energy[MAX_ROTAMER*(MAX_LEVEL+1)];

    for (int i = 0; i < tree_level; ++i) {
        m_node_offset[i]   = (ushort)node_offset[i];
        m_rot_per_level[i] = (uchar)rot_per_level[i];
    }
    for (int i = 0; i < rot_cnt; ++i)
        m_self_energy[i] = reduce_energy[rot_cnt][i] + reduce_energy[i][rot_cnt];

    printf("rot_per_level:\n");
    for (int i = 0; i < tree_level; ++i)
        printf("%d ", rot_per_level[i]);
    printf("\nnode_offset:\n");
    for (int i = 0; i < tree_level; ++i)
        printf("%d ", node_offset[i]);
    printf("\nself_energy:\n");
    for (int i = 0; i < rot_cnt; ++i)
        printf("%.3f ", m_self_energy[i]);
    printf("\n");


    float *p;

    p = m_reduce_energy;
    for (int i = 0; i < rot_cnt; ++i)
        for (int j = 0; j < rot_cnt; ++j)
            *p++ = reduce_energy[i][j];

    p = m_pm_energy;
    for (int i = 0; i < rot_cnt; ++i)
        for (int j = 0; j <= tree_level; ++j)
            *p++ = pm_energy[i][j];

    d_node_offset   = cl_create_buffer(context,
                                       CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                       sizeof(ushort) * tree_level,
                                       m_node_offset);
    d_rot_per_level = cl_create_buffer(context,
                                       CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                       sizeof(uchar) * tree_level,
                                       m_rot_per_level);
    d_self_energy   = cl_create_buffer(context,
                                       CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                       sizeof(float) * rot_cnt,
                                       m_self_energy);
    printf("rotcnt: %d\n", rot_cnt);

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
    /* these functions are needed to run OpenCL on nvidia graphic card */
    cl_image_format image_format = { CL_R, CL_FLOAT };
    d_reduce_energy = cl_create_image_2D(context,
                                      CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                      &image_format,
                                      rot_cnt,
                                      rot_cnt,
                                      0,
                                      m_reduce_energy);

    d_pm_energy     = cl_create_image_2D(context,
                                      CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                      &image_format,
                                      tree_level+1,
                                      rot_cnt,
                                      0,
                                      m_pm_energy);
#pragma GCC diagnostic pop

    char compile_option[200];
    snprintf(compile_option, 200,
             "-D NUM_LEVEL=%d -D NUM_ROTAMER=%d -D HEAP_SIZE=%d -D NUM_MAX_CHILD=%d -D GROUP_SIZE=%d",
             tree_level, rot_cnt, (int)HEAP_SIZE, max_child, (int)num_local_item);
    printf("Compiler option: %s\n", compile_option);

    program = cl_create_program_from_src(context,
                                         device,
                                         compile_option,
                                         ASTAR_KERNEL_SOURCE);

    kernel_initialize    = cl_create_kernel(program, "initialize");
    kernel_delete_min    = cl_create_kernel(program, "delete_min");
    kernel_radix_sort    = cl_create_kernel(program, "radix_sort");
    kernel_compute_score = cl_create_kernel(program, "compute_score");
    kernel_push_back     = cl_create_kernel(program, "push_back");

    cl_set_args(kernel_initialize, 0 , sizeof(d_heap), &d_heap);
    cl_set_args(kernel_initialize, 1 , sizeof(d_heap_size), &d_heap_size);
    cl_set_args(kernel_initialize, 2 , sizeof(d_data), &d_data);
    cl_set_args(kernel_initialize, 3 , sizeof(d_data_size), &d_data_size);
    cl_set_args(kernel_initialize, 4 , sizeof(d_optimal), &d_optimal);
    cl_set_args(kernel_initialize, 5 , sizeof(d_node_offset), &d_node_offset);
    cl_set_args(kernel_initialize, 6 , sizeof(d_rot_per_level), &d_rot_per_level);
    cl_set_args(kernel_initialize, 7 , sizeof(d_self_energy), &d_self_energy);
    cl_set_args(kernel_initialize, 8 , sizeof(d_reduce_energy), &d_reduce_energy);
    cl_set_args(kernel_initialize, 9 , sizeof(d_pm_energy), &d_pm_energy);

    cl_set_args(kernel_delete_min, 0 , sizeof(d_heap), &d_heap);
    cl_set_args(kernel_delete_min, 1 , sizeof(d_heap_size), &d_heap_size);
    cl_set_args(kernel_delete_min, 2 , sizeof(d_data), &d_data);
    cl_set_args(kernel_delete_min, 3 , sizeof(d_data_size), &d_data_size);
    cl_set_args(kernel_delete_min, 4 , sizeof(d_optimal), &d_optimal);
    cl_set_args(kernel_delete_min, 5 , sizeof(d_father), &d_father);
    cl_set_args(kernel_delete_min, 6 , sizeof(d_radix), &d_radix);
    cl_set_args(kernel_delete_min, 7 , sizeof(d_node_cnt), &d_node_cnt);
    cl_set_args(kernel_delete_min, 8 , sizeof(d_answer), &d_answer);
    cl_set_args(kernel_delete_min, 9 , sizeof(d_answer_size), &d_answer_size);
    cl_set_args(kernel_delete_min, 10, sizeof(d_rot_per_level), &d_rot_per_level);

    cl_set_args(kernel_radix_sort, 0 , sizeof(d_father), &d_father);
    cl_set_args(kernel_radix_sort, 1 , sizeof(d_input), &d_input);
    cl_set_args(kernel_radix_sort, 2 , sizeof(d_radix), &d_radix);
    cl_set_args(kernel_radix_sort, 3 , sizeof(d_node_cnt), &d_node_cnt);
    cl_set_args(kernel_radix_sort, 4 , sizeof(d_data), &d_data);
    cl_set_args(kernel_radix_sort, 5 , sizeof(d_rot_per_level), &d_rot_per_level);

    cl_set_args(kernel_compute_score, 0 , sizeof(d_input), &d_input);
    cl_set_args(kernel_compute_score, 1 , sizeof(d_data), &d_data);
    cl_set_args(kernel_compute_score, 2 , sizeof(d_data_size), &d_data_size);
    cl_set_args(kernel_compute_score, 3 , sizeof(d_output_size), &d_output_size);
    cl_set_args(kernel_compute_score, 4 , sizeof(d_node_offset), &d_node_offset);
    cl_set_args(kernel_compute_score, 5 , sizeof(d_rot_per_level), &d_rot_per_level);
    cl_set_args(kernel_compute_score, 6 , sizeof(d_self_energy), &d_self_energy);
    cl_set_args(kernel_compute_score, 7 , sizeof(d_reduce_energy), &d_reduce_energy);
    cl_set_args(kernel_compute_score, 8 , sizeof(d_pm_energy), &d_pm_energy);

    cl_set_args(kernel_push_back, 0 , sizeof(d_heap), &d_heap);
    cl_set_args(kernel_push_back, 1 , sizeof(d_heap_size), &d_heap_size);
    cl_set_args(kernel_push_back, 2 , sizeof(d_data), &d_data);
    cl_set_args(kernel_push_back, 3 , sizeof(d_data_size), &d_data_size);
    cl_set_args(kernel_push_back, 4 , sizeof(d_output_size), &d_output_size);
    cl_set_args(kernel_push_back, 5 , sizeof(d_begin_index), &d_begin_index);
    cl_set_args(kernel_push_back, 6 , sizeof(d_begin_index2), &d_begin_index2);
    cl_set_args(kernel_push_back, 7 , sizeof(d_optimal), &d_optimal);

    pq_init(heap);

    cl_finish(cmd_queue);
    puts("GPU setup finishes");
}

int *astar_gpu(bool first_run)
{
    static int rounds = 0;
    static int curr_conf[MAX_LEVEL];
    static uint optimal = 0;

    puts("GPU native A* start");

    wall_time_begin();
    if (first_run) {
        size_t global_size = rot_per_level[0];
        cl_launch_kernel(cmd_queue, kernel_initialize, 1, NULL,
                         &global_size, NULL, 0, NULL, NULL);
        cl_finish(cmd_queue);
        /* printf("initialize: %d\n", wall_time_elapsed()); */
    }

#define check_return() \
    if (!pq_empty(heap) && flip_float(pq_top(heap)->f_score) <= optimal) { \
        printf("GPU best result: %.3f\n", pq_top(heap)->f_score); \
        for (int i = 0; i < tree_level; ++i) \
            curr_conf[i] = pq_top(heap)->conf[i]; \
        free(pq_pop(heap, heap_cmp)); \
        printf("GPU native A* finish in %d rounds\n", rounds); \
        fflush(stdout); \
        return curr_conf; \
    } else

    check_return();

    static cl_event evt_cleanup;
    static cl_event evt_delete_min;
    static cl_event evt_radix_sort;
    static cl_event evt_compute_score;
    static cl_event evt_push_back;;


    static cl_uint  num_event = 0;
    static cl_event evt_wait[10];
    for (;;) {
        static int buf[MAX_LEVEL];
        cl_int num_nodes;
        cl_int data_size;
        cl_int output_size;

        cl_event evt_read;
        cl_event evt_write;

        ++rounds;

        /* fill zero to radix */
        cl_write_buffer(cmd_queue, d_radix, CL_FALSE,
                        0, tree_level*sizeof(int), buf,
                        num_event, num_event == 0 ? NULL : evt_wait, &evt_cleanup);
        num_event = 0;
        /* delete min */
        cl_launch_kernel(cmd_queue, kernel_delete_min, 1, NULL,
                         &num_global_item, &num_local_item,
                         1, &evt_cleanup, &evt_delete_min);
        /* printf("del-min: %d\n", wall_time_elapsed()); */

        /* read num_nodes and data_size */
        cl_read_buffer(cmd_queue, d_data_size, CL_TRUE,
                       0, sizeof(int), &data_size,
                       1, &evt_cleanup, &evt_read);
        cl_read_buffer(cmd_queue, d_radix, CL_TRUE,
                       (tree_level-1)*sizeof(int), sizeof(int), &num_nodes,
                       1, &evt_delete_min, &evt_write);
        output_size = data_size + num_nodes;
        printf("GPU data size: %dGPU output size: %d\n", data_size, output_size);
        cl_write_buffer(cmd_queue, d_output_size, CL_TRUE,
                        0, sizeof(int), &output_size,
                        0, NULL, &evt_write);
        /* printf("rrw: %d\n", wall_time_elapsed()); */

        /* radix sort */
        cl_launch_kernel(cmd_queue, kernel_radix_sort, 1, NULL,
                         &num_global_item, &num_local_item,
                         1, &evt_delete_min, &evt_radix_sort);
        /* printf("radix: %d\n", wall_time_elapsed()); */
        /* compute score */
        size_t snum_nodes = num_nodes;
        if (snum_nodes % num_local_item != 0)
            snum_nodes += num_local_item - snum_nodes % num_local_item;
        cl_launch_kernel(cmd_queue, kernel_compute_score, 1, NULL,
                         &snum_nodes, &num_local_item,
                         1, &evt_radix_sort, &evt_compute_score);
        /* printf("computer score: %d\n", wall_time_elapsed()); */

        /* push back */
        cl_launch_kernel(cmd_queue, kernel_push_back, 1, NULL,
                         &num_global_item, &num_local_item,
                         1, &evt_compute_score, &evt_push_back);

        cl_copy_buffer(cmd_queue, d_output_size, d_data_size, 0, 0, sizeof(int),
                       1, &evt_push_back, &evt_wait[num_event++]);
        cl_copy_buffer(cmd_queue, d_begin_index2, d_begin_index, 0, 0, sizeof(int),
                       1, &evt_push_back, &evt_wait[num_event++]);

        /* ------------------read/write answer_size------------------- */
        void *answer;
        int   answer_size;
        int  *answer_size_p;
        answer_size_p = cl_map_buffer(cmd_queue, d_answer_size, CL_TRUE,
                                      CL_MAP_READ | CL_MAP_WRITE, 0, sizeof(int),
                                      1, &evt_delete_min, NULL);
        answer_size = *answer_size_p;
        *answer_size_p = 0;
        cl_unmap_buffer(cmd_queue, d_answer_size, answer_size_p, 0, NULL, &evt_wait[num_event++]);
        /* ------------------ read answer ------------------- */
        if (answer_size) {
            answer = cl_map_buffer(cmd_queue, d_answer, CL_TRUE,
                                   CL_MAP_READ, 0, data_t_size * answer_size,
                                   0, NULL, NULL);
            void *p = answer;
            for (int i = 0; i < answer_size; ++i) {
                data_t *node = malloc(sizeof(data_t));

                node->f_score = *(float *)p;
                node->g_score = *(float *)(p+4);
                node->level = *(uchar *)(p+8);
                assert(node->f_score < -1.f);

                memcpy(node->conf, p+12, sizeof(uchar)*tree_level);

                pq_push(heap, node, heap_cmp);

                p += data_t_size;
            }
            cl_unmap_buffer(cmd_queue, d_answer, answer, 0, NULL, &evt_wait[num_event++]);
        }
        /* ------------------ read/write optimal ------------------- */
        uint *optimal_p;
        optimal_p = cl_map_buffer(cmd_queue, d_optimal, CL_TRUE,
                                  CL_MAP_READ | CL_MAP_WRITE, 0, sizeof(uint),
                                  1, &evt_push_back, NULL);
        optimal = *optimal_p;
        printf("optimal: %.3f\n", reverse_flip_float(optimal));
        *optimal_p = -1;
        cl_unmap_buffer(cmd_queue, d_optimal, optimal_p, 0, NULL, &evt_wait[num_event++]);
        /* ------------------ read host heap ------------------- */
        check_return();
        /* ------------------ print debug info ----------------- */
        /*
        cl_wait_for_events(1, &evt_push_back);
        printf("delete_min: %.3fms\tradix_sort: %.3fms\n"
               "compute_score: %.3fms\tpush_back: %.3fms\n"
               "read: %.3fms\twrite: %.3fms\tcleanup: %.3fms\n",
               cl_get_kernel_time(evt_delete_min),
               cl_get_kernel_time(evt_radix_sort),
               cl_get_kernel_time(evt_compute_score),
               cl_get_kernel_time(evt_push_back),
               cl_get_kernel_time(evt_read),
               cl_get_kernel_time(evt_write),
               cl_get_kernel_time(evt_cleanup));
               */
        fflush(stdout);
    }
}
