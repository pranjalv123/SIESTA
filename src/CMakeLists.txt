cmake_minimum_required (VERSION 2.6)
project (FastRFS)


link_directories(/usr/local/lib )
file(GLOB SOURCES *.cpp)
add_definitions(-std=c++14 -g -O3 -Wall)

option(ENABLE_PROFILING
  "Enable google profiler"
  OFF)

add_executable(FastRFS ${SOURCES})
target_link_libraries(FastRFS phylokit )

if(ENABLE_PROFILING)
  add_definitions(-DENABLE_PROFILING)
  target_link_libraries(FastRFS profiler)
endif(ENABLE_PROFILING)

if(UNIX)
  SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  SET(BUILD_SHARED_LIBRARIES OFF)
  SET(CMAKE_EXE_LINKER_FLAGS "-static -static-libgcc")
endif()

install(TARGETS FastRFS
		RUNTIME DESTINATION bin)
