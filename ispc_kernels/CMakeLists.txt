cmake_minimum_required(VERSION 3.5)

include(cmake/ispc.cmake)

project(libispc_objs LANGUAGES C)

# collect ispc files
file(GLOB_RECURSE ISPC_FILES  ${PROJECT_SOURCE_DIR}/*.ispc)

# call ispc preprocessor, collect *obj files in outListObj
add_ispc_src(outList outListObj "${ISPC_FILES}")

#add_custom_target(ISPC_OBJECTS ALL DEPENDS ${outListObj})

# make a library of all ispc *obj file
add_library(ispc_objs STATIC EXCLUDE_FROM_ALL ${outListObj})

set(ISPC_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}/CMakeFiles/ispc)
message(${ISPC_BUILD_DIR})

target_include_directories(ispc_objs
	PUBLIC
		$<BUILD_INTERFACE: ${ISPC_BUILD_DIR}>
)

#add_dependencies(ispc_objs ISPC_OBJECTS)
set_target_properties(ispc_objs PROPERTIES LINKER_LANGUAGE C)