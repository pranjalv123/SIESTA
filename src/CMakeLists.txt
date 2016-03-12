cmake_minimum_required (VERSION 2.6)
project (wASTRAL)

find_package(OpenMP)

link_directories(/usr/local/lib )
file(GLOB SOURCES *.cpp)
add_definitions(-std=c++14 -g -O3 -Wall ${OpenMP_CXX_FLAGS})

option(ENABLE_PROFILING
  "Enable google profiler"
  OFF)

add_executable(wASTRAL ${SOURCES})
target_link_libraries(wASTRAL pthread python2.7 ${OpenMP_CXX_FLAGS})

if(ENABLE_PROFILING)
  add_definitions(-DENABLE_PROFILING)
  target_link_libraries(wASTRAL profiler)
endif(ENABLE_PROFILING)

install(TARGETS wASTRAL
		RUNTIME DESTINATION bin)


