# from here http://programming-and-stuff.blogspot.com/2013/10/intel-spmd-ispc-with-cmake.html

find_program (ISPC_EXECUTABLE ispc)
if (NOT ISPC_EXECUTABLE)
    message(FATAL_ERROR "Failed to find ispc" )
endif()

#set (ISPC_IA_TARGETS "sse2-i32x4,sse4-i32x4,avx1-i32x8,avx2-i32x8,avx512knl-i32x16,avx512skx-i32x16" CACHE STRING "ISPC IA targets")
set (ISPC_IA_TARGETS "avx1" CACHE STRING "ISPC IA targets")

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

# given ispc file(s), append it and the generated .h file to 'outList'
macro( add_ispc_src outList outListObj srcLocationsList )
  set (BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}/CMakeFiles/ispc)
  file(MAKE_DIRECTORY ${BUILD_DIR}/ispc)

  foreach(ispc_file ${srcLocationsList})
    set(srcLocation2 ${ispc_file})
    find_file(srcLocation2 NAMES ${ispc_file} PATHS ${CMAKE_CURRENT_SOURCE_DIR} )
    get_source_file_property(srcLocation2 ${srcLocation2} LOCATION)
    string(FIND "${ispc_file}" "/" lastSlash REVERSE) #just get the file name, strip the path
    if(NOT ${lastSlash} STREQUAL "-1" )
      math(EXPR lastSlash ${lastSlash}+1)
      string(SUBSTRING "${ispc_file}" ${lastSlash} -1 objLocation)
    else()
      set( objLocation ${ispc_file} )
    endif()
    string(REPLACE ".ispc" ".obj" objLocation ${objLocation})
    set(objLocation ${BUILD_DIR}/${objLocation}) # put the obj in the intermediate dir
    string(REPLACE ".obj" "_ispc.h" hdrLocation ${objLocation})

    message("ISPC source file " ${ispc_file})
    #message("ISPC header " ${hdrLocation})
    #message("ISPC object " ${objLocation})
    #message("${ISPC_EXECUTABLE} ${srcLocation2} -o ${objLocation} -h ${hdrLocation} --arch=${ISPC_ARCH} --target=${ISPC_IA_TARGETS}")
    ADD_CUSTOM_COMMAND(
       OUTPUT ${objLocation}
       COMMAND ${ISPC_EXECUTABLE} ${srcLocation2} -o ${objLocation} -h ${hdrLocation} --arch=${ISPC_ARCH} --target=${ISPC_IA_TARGETS}
       IMPLICIT_DEPENDS C ${ispc_file}
       DEPENDS ${ispc_file}
       MAIN_DEPENDENCY  ${ispc_file}
       )
    set_source_files_properties( ${objLocation} PROPERTIES GENERATED TRUE EXTERNAL_OBJECT TRUE )
    set_source_files_properties( ${hdrLocation} PROPERTIES GENERATED TRUE EXTERNAL_OBJECT TRUE )
    list( APPEND ${outList} ${srcLocation} ${hdrLocation} )
    list( APPEND ${outListObj} ${objLocation})
  endforeach(ispc_file)
endmacro( add_ispc_src )