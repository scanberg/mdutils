### LIBXDRFILE ###

project(libxdrfile LANGUAGES C)

set(HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/src/xdrfile.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/xdrfile_trr.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/xdrfile_xtc.h
)
                 
set(SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/xdrfile.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/xdrfile_trr.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/xdrfile_xtc.c
)

if(MSVC)
add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif()

add_library(libxdrfile ${HEADERS} ${SOURCES})

target_include_directories(libxdrfile
	PUBLIC
		$<INSTALL_INTERFACE:/src>
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
)