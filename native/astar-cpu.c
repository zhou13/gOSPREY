// This file is not commented since it is similar to MSAStar.java.
// This unit is just used to compare the speed to astar-cuda.cu.  It is almost
// 100x faster than MSAStar.java.  But I haven't implement some feature such as
// first_run and single_seq because 1) It is slow. 2) I do not understand it.
#include "astar-cpu.h"
#include "qsort.h"

#define MAX_HEAP      20000000
#define HEAP_RESERVED 2000000
/*
#define MAX_HEAP      50
#define HEAP_RESERVED 20
*/

static int cur_conf[MAX_LEVEL];

static inline float get_pair_energy(int x, int y)
{
    return reduce_energy[x][y];
}
static inline float get_intra_energy(int idx)
{
    return reduce_energy[idx][rot_cnt];
}
static inline float get_shell_rot_energy(int idx)
{
    return reduce_energy[rot_cnt][idx];
}
static inline float get_total_energy(int idx)
{
    return get_intra_energy(idx) + get_shell_rot_energy(idx);
}

/* Heap <<< */
typedef struct node_t {
    int    level;     /* What level of tree we are in           */
    int    changes;   /* TODO How many rotamers we have changed */
    float  f_score;
    float  g_score;
    int   *conf;
} node_t;

static int heap_size;
static node_t *heap[MAX_HEAP+1];

static inline bool heap_cmp(node_t *a, node_t *b)
{
    return a->f_score < b->f_score;
}

QSORT_INIT(node_t *, heap_cmp, node);

static void heap_init(void)
{
    heap_size = 0;
}

static void heap_push(node_t *val)
{
    int now = ++heap_size;
    if (now >= MAX_HEAP)
        error("Cannot add more element.  Heap memory is too small!");
    while (now > 1) {
        int next = now / 2;
        if (heap_cmp(val, heap[next])) {
            heap[now] = heap[next];
            now = next;
        } else
            break;
    }
    heap[now] = val;
}

static node_t *heap_delete_min(void)
{
    if (heap_size == 0)
        return NULL;
    node_t *ret = heap[1];
    node_t *val = heap[heap_size--];
    int now = 1;

    int next;
    while ((next = now*2) <= heap_size) {
        if (next+1 <= heap_size && heap_cmp(heap[next+1], heap[next]))
            ++next;
        if (heap_cmp(heap[next], val)) {
            heap[now] = heap[next];
            now = next;
        } else
            break;
    }
    heap[now] = val;
    return ret;
}

/* >>> */

static bool is_pruned(int level, int idx)
{
    // TODO optimized by bitwise
    if (!do_pair_prune)
        return false;
    assert(false);
    for (int i = 0; i <= level; ++i)
        if (pair_pruned[idx][node_offset[i]+cur_conf[i]])
            return true;

    // TODO triple prune
    if (!do_triple_prune)
        return false;
    return false;
}

static float compute_h(int level) // predict energy
{
    float ret = 0;
    for (int i = level+1; i < tree_level; ++i) {
        float min_energy = FLT_MAX;
        for (int j = 0; j < rot_per_level[i]; ++j) {
            int idx = node_offset[i]+j;
            float cur_energy = get_total_energy(idx);

            for (int k = 0; k <= level; ++k)
                cur_energy += get_pair_energy(idx, node_offset[k]+cur_conf[k]);

            cur_energy += pm_energy[idx][i+1];

            min_energy = min(min_energy, cur_energy);
        }
        ret += min_energy;
    }
    return ret;
}

static void print_state(node_t *node)
{
    static int cnt = 0;
    ++cnt;
    if (cnt > 100)
        return;
    printf("level: %d\t f_score: %.9f g_score: %.9f h_score: %.9f\n",
           node->level, node->f_score, node->g_score, node->f_score - node->g_score);
    for (int i = 0; i <= node->level; ++i)
        printf("%d ", node->conf[i]);
    printf("\n");
}

static float compute_g(int level)  // current energy
{
    float ret = 0;

    for (int i = 0; i <= level; ++i) {
        int idx = node_offset[i] + cur_conf[i];
        for (int j = i+1; j <= level; ++j)
            ret += get_pair_energy(idx, node_offset[j]+cur_conf[j]);
        ret += get_total_energy(idx);
    }

    return ret;
}

static float compute_g_delta(int level)
{
    float ret = 0;
    int idx = node_offset[level] + cur_conf[level];

    for (int j = 0; j < level; ++j)
        ret += get_pair_energy(idx, node_offset[j]+cur_conf[j]);
    ret += get_total_energy(idx);

    return ret;
}

static node_t *new_node(int cur_level, float f_score, float g_score)
{
    node_t *ele = malloc(sizeof(node_t));

    ele->level = cur_level;
    ele->f_score = f_score;
    ele->g_score = g_score;

    size_t num_byte = sizeof(int)*(uint)(cur_level+1);
    ele->conf = malloc(num_byte);
    memcpy(ele->conf, cur_conf, num_byte);

    return ele;
}

static void free_node(node_t *node)
{
    free(node->conf);
    free(node);
}

static void shrink(void)
{
    if (heap_size + HEAP_RESERVED > MAX_HEAP) {
        puts("  >>> CPU run out of memory! Rescan to free memory");
        qsort_node(&heap[1], &heap[heap_size+1]);
        int old_heap_size = heap_size;
        int shrink_factor = max(5, rot_cnt / tree_level / 2);
        heap_size /= shrink_factor;
        printf("  >>> %d elements left\n", heap_size);

        for (int i = heap_size+1; i <= old_heap_size; ++i)
            free_node(heap[i]);
    }
}

void init_cpu(void)
{
    if (heap_size > 0)
        for (int i = 1; i <= heap_size; ++i)
            free_node(heap[i]);
    heap_init();
}

int *astar_cpu(bool first_run)
{
    puts("\n====== CPU Native A* begin ======");
    fflush(stdout);

    memset(cur_conf, -1, sizeof(int)*(uint)tree_level);

    if (first_run) {
        for (int i = 0; i < rot_per_level[0]; ++i) {
            cur_conf[0] = i;
            float g_score = compute_g(0);
            float f_score = g_score + compute_h(0);
            node_t *node = new_node(0, f_score, g_score);
            heap_push(node);
        }
        if (cur_conf[0] == -1)
            assert(false);
    }

    static int num_push = 0;
    static int max_level = 0;
    printf("current level:");
    for (;;) {
        node_t *cur_node = heap_delete_min();
        if (!cur_node) {
            memset(cur_conf, -1, sizeof(int)*(uint)tree_level);
            printf("\n");
            break;  // No solution
        }

        memcpy(cur_conf, cur_node->conf, sizeof(int)*(uint)tree_level);
        if (cur_node->level+1 == tree_level) {
            printf("\nbest result:%.9f\n", cur_node->f_score);
            print_state(cur_node);
            free(cur_node);
            break;  // Reach target
        }

        int cur_level = cur_node->level + 1;
        if (cur_level > max_level) {
            printf(" %d", cur_level);
	    fflush(stdout);
            max_level = cur_level;
        }
        for (int i = 0; i < rot_per_level[cur_level]; ++i) {
            if (is_pruned(cur_level-1, node_offset[cur_level]+i))
                continue;
            cur_conf[cur_level] = i;

            float g_score = cur_node->g_score + compute_g_delta(cur_level);
            float h_score = compute_h(cur_level);
            float f_score = g_score + h_score;

	    ++num_push;
            heap_push(new_node(cur_level, f_score, g_score));
        }

        free_node(cur_node);
        shrink();
    }
    printf(">>> CPU has elements: %d\n", num_push);
    puts("====== CPU native A* finished ======");

    return cur_conf;
}
