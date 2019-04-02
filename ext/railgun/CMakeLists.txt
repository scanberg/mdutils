### RAILGUN ###
project(railgun LANGUAGES C)

add_library(railgun STATIC ${CMAKE_CURRENT_SOURCE_DIR}/Railgun_Trolldom.c)
target_include_directories(railgun PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>)
target_compile_options(railgun PRIVATE /W3 /wd4101)