cmake_minimum_required(VERSION 3.5)

project(ispc_kernels)

find_program (ISPC_EXECUTABLE ispc)
if (NOT ISPC_EXECUTABLE)
    message(FATAL_ERROR "Failed to find ispc" )
endif()

set(ISPC_FLAGS "" CACHE STRING "ISPC Flags")
set(ISPC_IA_TARGETS "avx1-i32x8,avx2-i32x8,avx512knl-i32x16,avx512skx-i32x16" CACHE STRING "ISPC IA targets")
set(ISPC_KNOWN_TARGETS "sse2" "sse4" "avx1-" "avx1.1" "avx2" "avx512knl" "avx512skx")

if (UNIX)
    execute_process( COMMAND bash "-c" "uname -m | sed -e s/x86_64/x86/ -e s/i686/x86/ -e s/arm.*/arm/ -e s/sa110/arm/" OUTPUT_VARIABLE ARCH)
    string(STRIP ${ARCH} ARCH)
    execute_process( COMMAND getconf LONG_BIT OUTPUT_VARIABLE ARCH_BIT)
    string(STRIP ${ARCH_BIT} ARCH_BIT)
    if (${ARCH_BIT} EQUAL 32)
        set(ISPC_ARCH "x86")
    else()
        set(ISPC_ARCH "x86-64")
    endif()
else()
    set(ARCH "x86")
    if (CMAKE_SIZEOF_VOID_P EQUAL 8 )
        set(ISPC_ARCH "x86-64")
    else()
        set(ISPC_ARCH "x86")
    endif()
endif()

# collect ispc files
file(GLOB ISPC_FILES *.ispc)

# for each .ispc file
foreach(ISPC_FILE ${ISPC_FILES})
	get_filename_component(ISPC_SRC_NAME ${ISPC_FILE} NAME_WE)

	set(ISPC_HEADER_NAME "${CMAKE_CURRENT_BINARY_DIR}/${ISPC_SRC_NAME}_ispc.h")
	set(ISPC_OBJ_NAME "${CMAKE_CURRENT_BINARY_DIR}/${ISPC_SRC_NAME}_ispc${CMAKE_CXX_OUTPUT_EXTENSION}")
	
	#message("${ISPC_SRC_NAME}")
	#message("${ISPC_HEADER_NAME}")

	# Collect list of expected outputs
	list(APPEND ISPC_BUILD_OUTPUT ${ISPC_HEADER_NAME} ${ISPC_OBJ_NAME})
	list(APPEND ISPC_OBJ_FILES ${ISPC_OBJ_NAME})
    if ("${ARCH}" STREQUAL "x86")
        set(ISPC_TARGETS ${ISPC_IA_TARGETS})
        string(FIND ${ISPC_IA_TARGETS} "," MULTI_TARGET)
        if (${MULTI_TARGET} GREATER -1)
            foreach (ispc_target ${ISPC_KNOWN_TARGETS})
                string(FIND ${ISPC_IA_TARGETS} ${ispc_target} FOUND_TARGET)
                if (${FOUND_TARGET} GREATER -1)
                    set(OUTPUT_TARGET ${ispc_target})
                    if (${ispc_target} STREQUAL "avx1-")
                        set(OUTPUT_TARGET "avx")
                    elseif (${ispc_target} STREQUAL "avx1.1")
                        set(OUTPUT_TARGET "avx11")
                    endif()
                    list(APPEND ISPC_BUILD_OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${ISPC_SRC_NAME}_ispc_${OUTPUT_TARGET}.h"
                                "${CMAKE_CURRENT_BINARY_DIR}/${ISPC_SRC_NAME}_ispc_${OUTPUT_TARGET}${CMAKE_CXX_OUTPUT_EXTENSION}")
                    list(APPEND ISPC_OBJ_FILES "${CMAKE_CURRENT_BINARY_DIR}/${ISPC_SRC_NAME}_ispc_${OUTPUT_TARGET}${CMAKE_CXX_OUTPUT_EXTENSION}")
                endif()
            endforeach()
        endif()
    elseif ("${ARCH}" STREQUAL "arm")
        set(ISPC_TARGETS ${example_ISPC_ARM_TARGETS})
    else()
        message(FATAL_ERROR "Unknown architecture ${ARCH}")
    endif()
	add_custom_command(OUTPUT ${ISPC_OBJ_NAME}
	    COMMAND ${ISPC_EXECUTABLE} ${ISPC_FILE} ${ISPC_FLAGS} --target=${ISPC_TARGETS} --arch=${ISPC_ARCH}
	                                -h ${ISPC_HEADER_NAME} -o ${ISPC_OBJ_NAME}
	    VERBATIM
	    DEPENDS ${ISPC_FILE})
endforeach(ISPC_FILE)
add_custom_target(ispc_objects ALL DEPENDS ${ISPC_OBJ_FILES})

add_library(ispc_kernels STATIC EXCLUDE_FROM_ALL ${ISPC_FILES} ${ISPC_OBJ_FILES})
add_dependencies(ispc_kernels ispc_objects)

target_include_directories(ispc_kernels
	PUBLIC
		$<BUILD_INTERFACE: ${CMAKE_CURRENT_BINARY_DIR}>
)

set_target_properties(ispc_kernels PROPERTIES LINKER_LANGUAGE C)

# To show ispc source in VS solution:
if (WIN32)
    set_source_files_properties(${ISPC_FILES} PROPERTIES HEADER_FILE_ONLY TRUE)
endif()

if(MSVC)
    # Group ISPC files inside Visual Studio
    source_group("ispc" FILES ${ISPC_FILES})
endif()