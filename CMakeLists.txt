cmake_minimum_required(VERSION 2.8)

project(gOSPREY)

include(UseJava)
find_package(Java REQUIRED)

file(GLOB osprey_java_SRC
    "src/*.java"
    "src/mpi/*.java")

add_jar(osprey ${osprey_java_SRC})

# Try to check javah
FIND_PROGRAM( JAVA_HEADER_COMPILE javah
    PATHS
    "[HKEY_LOCAL_MACHINE\\SOFTWARE\\JavaSoft\\Java Development Kit\\1.4;JavaHome]/bin"
    "[HKEY_LOCAL_MACHINE\\SOFTWARE\\JavaSoft\\Java Development Kit\\1.3;JavaHome]/bin"
    /usr/bin
    /usr/lib/java/bin
    /usr/share/java/bin
    /usr/local/bin
    /usr/local/java/bin
)
MARK_AS_ADVANCED( JAVA_HEADER_COMPILE )

IF( NOT JAVA_HEADER_COMPILE )
    MESSAGE( FATAL_ERROR "javah cannot be found on you system." )
ENDIF( NOT JAVA_HEADER_COMPILE )
# End of the check


add_custom_command(
    OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/native/MSAStar.h
    COMMAND javah -classpath ${CMAKE_CURRENT_BINARY_DIR}/osprey.jar -o ${CMAKE_CURRENT_BINARY_DIR}/native/MSAStar.h MSAStar
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/MSAStar.java)
add_custom_target(
    native_header
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/native/MSAStar.h)
add_subdirectory(native)
add_dependencies(native_header osprey)
add_dependencies(MSAStar native_header)
