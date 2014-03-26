#include "clwrapper.h"

const char *cl_get_error_msg(cl_int status) // <<<
{
    switch (status) {
    case 0:
        return "CL_SUCCESS";
    case -1:
        return "CL_DEVICE_NOT_FOUND";
    case -2:
        return "CL_DEVICE_NOT_AVAILABLE";
    case -3:
        return "CL_COMPILER_NOT_AVAILABLE";
    case -4:
        return "CL_MEM_OBJECT_ALLOCATION_FAILURE";
    case -5:
        return "CL_OUT_OF_RESOURCES";
    case -6:
        return "CL_OUT_OF_HOST_MEMORY";
    case -7:
        return "CL_PROFILING_INFO_NOT_AVAILABLE";
    case -8:
        return "CL_MEM_COPY_OVERLAP";
    case -9:
        return "CL_IMAGE_FORMAT_MISMATCH";
    case -10:
        return "CL_IMAGE_FORMAT_NOT_SUPPORTED";
    case -11:
        return "CL_BUILD_PROGRAM_FAILURE";
    case -12:
        return "CL_MAP_FAILURE";
    case -13:
        return "CL_MISALIGNED_SUB_BUFFER_OFFSET";
    case -14:
        return "CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST";
    case -15:
        return "CL_COMPILE_PROGRAM_FAILURE";
    case -16:
        return "CL_LINKER_NOT_AVAILABLE";
    case -17:
        return "CL_LINK_PROGRAM_FAILURE";
    case -18:
        return "CL_DEVICE_PARTITION_FAILED";
    case -19:
        return "CL_KERNEL_ARG_INFO_NOT_AVAILABLE";
    case -30:
        return "CL_INVALID_VALUE";
    case -31:
        return "CL_INVALID_DEVICE_TYPE";
    case -32:
        return "CL_INVALID_PLATFORM";
    case -33:
        return "CL_INVALID_DEVICE";
    case -34:
        return "CL_INVALID_CONTEXT";
    case -35:
        return "CL_INVALID_QUEUE_PROPERTIES";
    case -36:
        return "CL_INVALID_COMMAND_QUEUE";
    case -37:
        return "CL_INVALID_HOST_PTR";
    case -38:
        return "CL_INVALID_MEM_OBJECT";
    case -39:
        return "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR";
    case -40:
        return "CL_INVALID_IMAGE_SIZE";
    case -41:
        return "CL_INVALID_SAMPLER";
    case -42:
        return "CL_INVALID_BINARY";
    case -43:
        return "CL_INVALID_BUILD_OPTIONS";
    case -44:
        return "CL_INVALID_PROGRAM";
    case -45:
        return "CL_INVALID_PROGRAM_EXECUTABLE";
    case -46:
        return "CL_INVALID_KERNEL_NAME";
    case -47:
        return "CL_INVALID_KERNEL_DEFINITION";
    case -48:
        return "CL_INVALID_KERNEL";
    case -49:
        return "CL_INVALID_ARG_INDEX";
    case -50:
        return "CL_INVALID_ARG_VALUE";
    case -51:
        return "CL_INVALID_ARG_SIZE";
    case -52:
        return "CL_INVALID_KERNEL_ARGS";
    case -53:
        return "CL_INVALID_WORK_DIMENSION";
    case -54:
        return "CL_INVALID_WORK_GROUP_SIZE";
    case -55:
        return "CL_INVALID_WORK_ITEM_SIZE";
    case -56:
        return "CL_INVALID_GLOBAL_OFFSET";
    case -57:
        return "CL_INVALID_EVENT_WAIT_LIST";
    case -58:
        return "CL_INVALID_EVENT";
    case -59:
        return "CL_INVALID_OPERATION";
    case -60:
        return "CL_INVALID_GL_OBJECT";
    case -61:
        return "CL_INVALID_BUFFER_SIZE";
    case -62:
        return "CL_INVALID_MIP_LEVEL";
    case -63:
        return "CL_INVALID_GLOBAL_WORK_SIZE";
    case -64:
        return "CL_INVALID_PROPERTY";
    case -65:
        return "CL_INVALID_IMAGE_DESCRIPTOR";
    case -66:
        return "CL_INVALID_COMPILER_OPTIONS";
    case -67:
        return "CL_INVALID_LINKER_OPTIONS";
    case -68:
        return "CL_INVALID_DEVICE_PARTITION_COUNT";
    }
    return "WTF";
} // >>>

cl_device_id cl_get_first_gpu_device(void)
{
    cl_uint num_platforms;
    cl_platform_id platform;

    cl_uint num_devices;
    cl_device_id device;
    

    CHK_CL(clGetPlatformIDs(0, NULL, &num_platforms));
    CHK_CD(num_platforms == 0, "Error: Number of platform is zero!\n");

    CHK_CL(clGetPlatformIDs(1, &platform, NULL));
    CHK_CL(clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 0, NULL, &num_devices));
    CHK_CD(num_devices == 0, "Error: No GPU device for OpenCL!\n");

    CHK_CL(clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 1, &device, NULL));
    return device;
}

cl_program cl_create_program_from_src(cl_context context,
                                      cl_device_id device,
                                      const char *options,
                                      const char *source)
{
    cl_int status;
    cl_program program= clCreateProgramWithSource(context, 1, &source, 0, &status);
    CHK_CL(status);

    status = clBuildProgram(program, 1, &device, options, NULL, NULL);

    if (status != CL_BUILD_PROGRAM_FAILURE)
        CHK_CL(status);

    cl_build_status build_status;
    CHK_CL(clGetProgramBuildInfo(program, device,
                                 CL_PROGRAM_BUILD_STATUS,
                                 sizeof(cl_build_status),
                                 &build_status,
                                 NULL));
    switch (build_status) {
    case CL_BUILD_SUCCESS:
        return program;

    case CL_BUILD_NONE:
    case CL_BUILD_IN_PROGRESS:
        assert(0);

    case CL_BUILD_ERROR: ;
        size_t size;
        CHK_CL(clGetProgramBuildInfo(program, device,
                                     CL_PROGRAM_BUILD_LOG,
                                     0,
                                     NULL,
                                     &size));
        char *log_msg = malloc(size+1);
        CHK_CL(clGetProgramBuildInfo(program, device,
                                     CL_PROGRAM_BUILD_LOG,
                                     size,
                                     log_msg,
                                     NULL));
        log_msg[size] = 0;
        printf("Build log:\n");
        puts(log_msg);
        exit(-1);
    }
    exit(-1);
}

cl_double cl_get_kernel_time(cl_event event)
{
    cl_ulong time_start;
    cl_ulong time_end;

    CHK_CL(clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(time_start), &time_start, NULL));
    CHK_CL(clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(time_end), &time_end, NULL));

    return (cl_double)(time_end - time_start) / 1000000.;
}

cl_double cl_get_kernel_bandwidth(cl_event event, size_t data_size)
{
    return (cl_double)(data_size) / 1024. / 1024. / 1024. / (cl_get_kernel_time(event) / 1000.);
}
