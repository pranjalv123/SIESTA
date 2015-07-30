cmake_minimum_required (VERSION 2.6)
project (wASTRAL)

link_directories(/usr/local/lib )
file(GLOB SOURCES *.cpp)
add_definitions(-std=c++14 -g -O3)

add_executable(wASTRAL ${SOURCES})
target_link_libraries(wASTRAL profiler pthread )





