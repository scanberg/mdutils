cmake_minimum_required(VERSION 3.5)

project(libmdutils VERSION 0.7.7 LANGUAGES CXX)

option(MDUTILS_USE_RELATIVE_PATHS "Use relative paths for resource directories" OFF)

file(GLOB CORE_FILES src/core/*)
file(GLOB GFX_FILES src/gfx/*)
file(GLOB MOL_FILES src/mol/*)

source_group("core" FILES ${CORE_FILES})
source_group("gfx" FILES ${GFX_FILES})
source_group("mol" FILES ${MOL_FILES})

add_subdirectory(ispc_kernels)
add_subdirectory(ext/gl3w)
add_subdirectory(ext/libxdrfile)
add_subdirectory(ext/stb)

message("${ISPC_BUILD_DIR}")

add_library(mdutils ${CORE_FILES} ${GFX_FILES} ${MOL_FILES})

target_include_directories(mdutils
	PUBLIC
		$<INSTALL_INTERFACE:src;ext/glm;ext/gl3w>
		$<BUILD_INTERFACE:
			${CMAKE_CURRENT_SOURCE_DIR}/src;
			${CMAKE_CURRENT_SOURCE_DIR}/ext/glm;
			${CMAKE_CURRENT_SOURCE_DIR}/ext/gl3w
		>
	PRIVATE
		${CMAKE_CURRENT_SOURCE_DIR}/src
		${CMAKE_CURRENT_SOURCE_DIR}/ext/glm
		${CMAKE_CURRENT_SOURCE_DIR}/ext/libxdrfile/src
		${CMAKE_CURRENT_SOURCE_DIR}/ext/meow_hash
		${ISPC_BUILD_DIR}
)

set(RESOURCE_DIR "${PROJECT_SOURCE_DIR}/")
if (MDUTILS_USE_RELATIVE_PATHS)
    set(RESOURCE_DIR "")
endif()

target_compile_definitions(mdutils PRIVATE
	MDUTILS_SHADER_DIR=\"${RESOURCE_DIR}shaders\"
	MDUTILS_IMAGE_DIR=\"${RESOURCE_DIR}images\")

target_compile_features(mdutils PRIVATE cxx_std_17)

target_link_libraries(mdutils
	PUBLIC
		gl3w
    PRIVATE
		libxdrfile
		stb
		ispc_objs
)