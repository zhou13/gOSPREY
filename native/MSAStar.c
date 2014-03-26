#include "astar.h"

#include "astar-cpu.h"

#include "astar-opencl.h"
#include "astar-cuda.h"

#include "MSAStar.h"

int tree_level;
int rot_cnt;

bool do_pert;

int node_offset[MAX_LEVEL];
int  rot_per_level[MAX_LEVEL];

bool do_pair_prune;
bool **pair_pruned;

bool do_triple_prune;
bool ***triple_pruned;

bool enable_cpu;
bool enable_gpu;
size_t cpu_memory;
size_t gpu_memory;
int num_gpu_group;
int num_gpu_item;
int num_gpu_item2;
double shrink_ratio;

/* Now we suppose reduce energy matrix is a N+1 x N+1 */
float reduce_energy[MAX_ROTAMER][MAX_ROTAMER];
float pm_energy[MAX_ROTAMER][MAX_LEVEL+1];

STATIC_ASSERT(sizeof(bool) == sizeof(jboolean));
STATIC_ASSERT(sizeof(int) == sizeof(jint));
STATIC_ASSERT(sizeof(float) == sizeof(jfloat));

JNIEXPORT void JNICALL Java_MSAStar_initNativeAStar(
        JNIEnv *env,
        jobject self,
        jint _tree_level,
        jintArray _rot_per_level,
        jobjectArray _arp_matrix,
        jobject _steric_check,     /* TODO */
        jobjectArray _sp_flags,    /* TODO */
        jobjectArray _trip_flags,  /* TODO */
        jboolean _do_pert,         /* TODO */
        jboolean _enable_cpu,
        jboolean _enable_gpu,
        jlong _cpu_memory,
        jlong _gpu_memory,
        int _num_gpu_group,
        int _num_gpu_item,
        int _num_gpu_item2,
        double _shrink_ratio)
{
    tree_level = _tree_level;

    enable_cpu = _enable_cpu;
    enable_gpu = _enable_gpu;
    cpu_memory = _cpu_memory;
    gpu_memory = _gpu_memory;
    printf("Allowed memory usage: CPU: %zuMB, GPU: %zuMB\n",
           cpu_memory/1024/1024, gpu_memory/1024/1024);
    num_gpu_group = _num_gpu_group;
    num_gpu_item = _num_gpu_item;
    num_gpu_item2 = _num_gpu_item2;
    printf("Number of work-group: %d\n", num_gpu_group);
    printf("Number of work-item: %d\n", num_gpu_item);
    shrink_ratio = _shrink_ratio;
    printf("Shrink ratio: %.2f\n", shrink_ratio);

    rot_cnt = (*env)->GetArrayLength(env, _arp_matrix)-1;
    printf("rot_cnt: %d\n", rot_cnt);fflush(stdout);
    for (int i = 0; i <= rot_cnt; ++i) {
        jfloatArray t = (jfloatArray)(*env)->
            GetObjectArrayElement(env, _arp_matrix, i);
        assert(t);
        jfloat *fptr = (*env)->GetFloatArrayElements(env, t, NULL);
        memcpy(reduce_energy+i, fptr, sizeof(float)*(uint)(rot_cnt+1));
        (*env)->ReleaseFloatArrayElements(env, t, fptr, 0);
    }

    jint *iptr = (*env)->GetIntArrayElements(env, _rot_per_level, NULL);
    memcpy(rot_per_level, iptr, sizeof(int)*(uint)(tree_level));
    (*env)->ReleaseIntArrayElements(env, _rot_per_level, iptr, 0);

    int64_t space_size = 1;
    for (int i = 0; i < tree_level; ++i)
        space_size *= rot_per_level[i];
    printf("size of conformation space: %" PRId64 "\n", space_size);
    fflush(stdout);

    node_offset[0] = 0;
    for (int i = 1; i < tree_level; ++i)
        node_offset[i] = node_offset[i-1] + rot_per_level[i-1];
    printf("node offset OK\n");
    fflush(stdout);

    for (int i = 0; i < tree_level; ++i)
        for (int j = 0; j < rot_per_level[i]; ++j) {
            int idx = node_offset[i]+j;
            for (int k = 0; k < tree_level; ++k) {
                pm_energy[idx][k] = FLT_MAX;
                for (int p = 0; p < rot_per_level[k]; ++p)
                    pm_energy[idx][k] = min(pm_energy[idx][k],
                                            reduce_energy[idx][node_offset[k]+p]);
            }
            pm_energy[idx][tree_level] = 0;
            for (int k = tree_level-2; k >= 0; --k)
                pm_energy[idx][k] += pm_energy[idx][k+1];
        }
    printf("pm_energy OK\n");
    fflush(stdout);

    /*
    printf("rot_cnt: %d\n", rot_cnt);
    puts("Node offset:");
    for (int i = 0; i < tree_level; ++i)
        printf("%d %d|", node_offset[i], rot_per_level[i]);
    puts("");
    */

    if (enable_cpu)
        init_cpu();
    if (enable_gpu)
        init_cuda();
    // opencl module is not maintained due to lack of devices
    // init_opencl();

    if (!_sp_flags)
        do_pair_prune = false;
    else {
        do_pair_prune = true;
        for (int i = 0; i < rot_cnt; ++i) {
            jbooleanArray t = (jbooleanArray)(*env)->
                GetObjectArrayElement(env, _sp_flags, i);
            jboolean *bptr = (*env)->GetBooleanArrayElements(env, t, NULL);
            memcpy(pair_pruned+i, bptr, sizeof(bool)*(uint)(rot_cnt));
            (*env)->ReleaseBooleanArrayElements(env, t, bptr, 0);
        }
    }
    do_pert = _do_pert;

    puts("Return from JNI!");
    fflush(stdout);
    fflush(stderr);
}


// Note: Almost all the parameter is not used
JNIEXPORT jintArray JNICALL Java_MSAStar_doNativeAStar(
        JNIEnv        *env,
        jobject        self,
        jboolean       first_run,
        jint           max_change,          /* TODO */
        jintArray      _nodes_def,          /* TODO */
        jbooleanArray  _pruned_nodes,       /* TODO */
        jobjectArray   _strand_rot,         /* TODO */
        jobjectArray   _strand_def,         /* TODO */
        jintArray      _res_cnt,            /* TODO */
        jobjectArray   _strand_mut,         /* TODO */
        jboolean       single_seq,          /* TODO */
        jintArray      _mut_strand,         /* TODO */
        jintArray      _mut_mut)            /* TODO */
{
    puts("Comming back to MSAStar.c!");
    fflush(stdout);

    static int cntEnter = 0;
    if (++cntEnter > MAX_ROUND) {
        printf("MSAStar.c Line %d:  exceed MAX_ROUND, exit!\n", __LINE__);
        exit(0);
    }

    jint *conf_cpu  = NULL;
    jint *conf_cuda = NULL;

    if (enable_gpu) {
        wall_time_begin();
        conf_cuda = (jint*)astar_cuda((bool)first_run);
        printf(">>> Native CUDA runs in %" PRId64 " ms \n\n",
               wall_time_elapsed());
        conf_cpu = conf_cuda;
    }

    if (enable_cpu) {
        wall_time_begin();
        conf_cpu = (jint*)astar_cpu((bool)first_run);
        printf(">>> Native CPU runs in %" PRId64 " ms \n\n",
               wall_time_elapsed());

        if (enable_gpu) {
            for (int i = 0; i < tree_level; ++i) {
                if (conf_cpu[i] != conf_cuda[i]) {
                    printf("WARNING: Two Native A* result is different.  It may be due to float precision.\n");
                    printf("conf[%d], cpu: %d; gpu: %d\n", i, conf_cpu[i], conf_cuda[i]);
                }
	    }
        }
    }


    // OpenCL modules are not maintained
    /*
    jint *conf_opencl = (jint*)astar_opencl((bool)first_run);
    for (int i = 0; i < tree_level; ++i)
        assert(conf_cpu[i] == conf_opencl[i]);
    */

    jintArray ret = (*env)->NewIntArray(env, tree_level);
    if (ret == NULL)
        error("Failed to alloc a new java int array");
    (*env)->SetIntArrayRegion(env, ret, 0, tree_level, (jint*)conf_cpu);

    fflush(stdout);
    fflush(stderr);
    return ret;
}
