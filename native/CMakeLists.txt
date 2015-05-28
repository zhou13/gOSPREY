project(OSPREY-Native)

find_package(JNI REQUIRED)
find_package(CUDA REQUIRED)

set(CUDA_NVCC_FLAGS "-arch=sm_20 -Xcompiler -fPIC")
set(JNI_HEADER ${CMAKE_CURRENT_BINARY_DIR}/MSAStar.h)
cuda_add_library(cuda_astar STATIC astar-cuda.cu)

add_definitions(-std=gnu99)
set(osprey_native_SRC
    astar-cpu.c
    MSAStar.c
)

include_directories(${JNI_INCLUDE_DIRS})

include_directories(${CMAKE_CURRENT_BINARY_DIR})
set_source_files_properties(${JNI_HEADER} PROPERTIES GENERATED ON)
add_library(MSAStar SHARED ${osprey_native_SRC})
target_link_libraries(MSAStar cuda_astar)

install(TARGETS MSAStar DESTINATION lib)
