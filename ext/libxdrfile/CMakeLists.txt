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

add_library(libxdrfile ${HEADERS} ${SOURCES})

if(MSVC)
	target_compile_definitions(libxdrfile
		PRIVATE _CRT_SECURE_NO_WARNINGS)
	target_compile_options(libxdrfile
		PRIVATE /W3 /wd4244 /wd4267 /wd4018 /wd4700)
endif()

target_include_directories(libxdrfile
	PUBLIC
		$<INSTALL_INTERFACE:/src>
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
)