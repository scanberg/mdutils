cmake_minimum_required(VERSION 3.5)

project(libmdutils VERSION 0.1.1 LANGUAGES CXX)

file(GLOB CORE_FILES src/core/*)
file(GLOB MOL_FILES src/mol/*)

source_group("core" FILES ${CORE_FILES})
source_group("mol" FILES ${MOL_FILES})

add_subdirectory(ext/xdrfile)
add_subdirectory(ext/railgun)

#message("${ISPC_BUILD_DIR}")

add_library(mdutils ${CORE_FILES} ${MOL_FILES})

target_include_directories(mdutils
	PUBLIC
		$<INSTALL_INTERFACE:src;ext/glm;>
		$<BUILD_INTERFACE:
			${CMAKE_CURRENT_SOURCE_DIR}/src;
			${CMAKE_CURRENT_SOURCE_DIR}/ext/glm;
			${CMAKE_CURRENT_SOURCE_DIR}/ext/svd3
		>
	PRIVATE
		${CMAKE_CURRENT_SOURCE_DIR}/src
		${CMAKE_CURRENT_SOURCE_DIR}/ext/glm
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
    PRIVATE
		xdrfile
		railgun
)