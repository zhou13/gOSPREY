#pragma once

#include <CL/cl.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

// check condition
#define CHK_CD(expression, msg, ...) \
    if (expression) { \
        fprintf(stderr, msg, ##__VA_ARGS__); \
        exit(-1); \
    } else ;

// check opencl error
#define CHK_CL(__cl_wrapper_status) ({ \
    cl_int __my_cl___cl_wrapper_status__ = (__cl_wrapper_status); \
    if (__my_cl___cl_wrapper_status__ != CL_SUCCESS) { \
        fprintf(stderr, "[ERROR] (%s: %d): %s\n", __FILE__, __LINE__, \
                cl_get_error_msg(__my_cl___cl_wrapper_status__)); \
        exit(-1); \
    } \
})

const char *cl_get_error_msg(cl_int __cl_wrapper_status);
cl_device_id cl_get_first_gpu_device(void);
cl_program cl_create_program_from_src(cl_context context,
                                      cl_device_id device,
                                      const char *options,
                                      const char *source);

#define cl_query_device(device, param, buf_size, buf, real_size) \
    CHK_CL(clGetDeviceInfo(device, param, buf_size, buf, real_size));

#define cl_create_context(device) ({ \
    cl_int __cl_wrapper_status; \
    cl_context context = clCreateContext(NULL, 1, &(device), NULL, NULL, &__cl_wrapper_status); \
    CHK_CL(__cl_wrapper_status); \
    context; })

#define cl_create_cmd_queue(context, device, properties) ({ \
    cl_int __cl_wrapper_status; \
    cl_command_queue cmd_queue = clCreateCommandQueue(context, device, properties, &__cl_wrapper_status); \
    CHK_CL(__cl_wrapper_status); \
    cmd_queue; })

#define cl_create_kernel(program, kernel_name) ({ \
    cl_int __cl_wrapper_status; \
    cl_kernel kernel = clCreateKernel(program, kernel_name, &__cl_wrapper_status); \
    CHK_CL(__cl_wrapper_status); \
    kernel; })

#define cl_create_buffer(context, flags, size, host_ptr) ({ \
    cl_int __cl_wrapper_status; \
    cl_mem mem = clCreateBuffer(context, flags, size, host_ptr, &__cl_wrapper_status); \
    CHK_CL(__cl_wrapper_status); \
    mem; })

#define cl_fill_buffer(cmd_queue, buffer, pattern, pattern_size, \
                       offset, size, num_wait, wait, event) \
    CHK_CL(clEnqueueFillBuffer(cmd_queue, buffer, pattern, pattern_size, \
                               offset, size, num_wait, wait, event))

#define cl_set_args(kernel, arg_index, arg_size, arg_value) ({ \
    CHK_CL(clSetKernelArg(kernel, arg_index, arg_size, arg_value)); }) \

#define cl_launch_kernel(cmd_queue, kernel, dim, global_work_offset, \
                         global_work_size, local_work_size, \
                         num_wait, wait, event) \
    CHK_CL(clEnqueueNDRangeKernel(cmd_queue, kernel, dim, global_work_offset, \
                                  global_work_size, local_work_size, \
                                  num_wait, wait, event)); 

#define cl_finish(cmd_queue) CHK_CL(clFinish(cmd_queue))

#define cl_map_buffer(cmd_queue, buffer, blocking_map, map_flags, offset, cb, \
                      num_wait, wait, event) ({ \
    cl_int __cl_wrapper_status; \
    void *__cl_wrapper_ret = clEnqueueMapBuffer(cmd_queue, buffer, \
                                   blocking_map, map_flags, offset, cb, \
                                   num_wait, wait, event, &__cl_wrapper_status); \
    CHK_CL(__cl_wrapper_status); \
    __cl_wrapper_ret; })

#define cl_read_buffer(cmd_queue, buffer, blocking_read, offset, cb, ptr, \
                       num_wait, wait, event) ({ \
    CHK_CL(clEnqueueReadBuffer(cmd_queue, buffer, blocking_read, offset, cb, ptr, \
                                num_wait, wait, event)); })

#define cl_write_buffer(cmd_queue, buffer, blocking_write, offset, cb, ptr, \
                        num_wait, wait, event) ({ \
    CHK_CL(clEnqueueWriteBuffer(cmd_queue, buffer, blocking_write, offset, cb, ptr, \
                                num_wait, wait, event)); })

#define cl_copy_buffer(cmd_queue, src, dst, src_off, dst_off, cb, \
                       num_wait, wait, event) ({ \
    CHK_CL(clEnqueueCopyBuffer(cmd_queue, src, dst, src_off, dst_off, cb, \
                               num_wait, wait, event)); })

#define cl_unmap_buffer(cmd_queue, buffer, ptr, num_wait, wait, event) \
    CHK_CL(clEnqueueUnmapMemObject(cmd_queue, buffer, ptr, num_wait, wait, event))

#define cl_wait_for_events(num_events, event_list) \
    CHK_CL(clWaitForEvents(num_events, event_list));

#define cl_create_image(context, flags, format, desc, ptr) ({ \
    cl_int __cl_wrapper_status; \
    cl_mem __cl_wrapper_ret = clCreateImage(context, flags, format, desc, ptr, &__cl_wrapper_status); \
    CHK_CL(__cl_wrapper_status); \
    __cl_wrapper_ret; })

#define cl_create_image_2D(context, flags, format, \
                          width, height, row_pitch, ptr) ({ \
    cl_int __cl_wrapper_status; \
    cl_mem __cl_wrapper_ret = clCreateImage2D(context, flags, format, \
                                              width, height, row_pitch, \
                                              ptr, &__cl_wrapper_status); \
    CHK_CL(__cl_wrapper_status); \
    __cl_wrapper_ret; })

#define cl_create_image_3D(context, flags, format, \
                          width, height, depth, row_pitch, slice_pitch, ptr) ({ \
    cl_int __cl_wrapper_status; \
    cl_mem __cl_wrapper_ret = clCreateImage3D(context, flags, format, \
                                              width, height, depth, \
                                              row_pitch, slice_pitch, \
                                              ptr, &__cl_wrapper_status); \
    CHK_CL(__cl_wrapper_status); \
    __cl_wrapper_ret; })

/* in ms */
cl_double cl_get_kernel_time(cl_event event);
/* in GB/s */
cl_double cl_get_kernel_bandwidth(cl_event event, size_t data_size);
