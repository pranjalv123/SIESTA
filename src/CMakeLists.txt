cmake_minimum_required (VERSION 2.6)
project (wASTRAL)

find_package(OpenMP REQUIRED)
find_package(PythonLibs 2.7 REQUIRED)

link_directories(/usr/local/lib )
file(GLOB SOURCES *.cpp)
add_definitions(-std=c++11 -g -O3 -Wall ${OpenMP_CXX_FLAGS})

option(ENABLE_PROFILING
  "Enable google profiler"
  OFF)

add_executable(wASTRAL ${SOURCES})
target_link_libraries(wASTRAL pthread ${OpenMP_CXX_FLAGS} ${PYTHON_LIBRARIES})

if(ENABLE_PROFILING)
  add_definitions(-DENABLE_PROFILING)
  target_link_libraries(wASTRAL profiler)
endif(ENABLE_PROFILING)

install(TARGETS wASTRAL
		RUNTIME DESTINATION bin)


