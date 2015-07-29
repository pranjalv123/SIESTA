cmake_minimum_required (VERSION 2.6)
project (wASTRAL)


find_package(Boost 1.55.0 COMPONENTS program_options REQUIRED)
link_directories(/usr/local/lib )
link_directories(${Boost_LIBRARY_DIRS})
message("FOOBAR Boost found" ${Boost_FOUND})
message("FOOBAR Boost library dir" ${Boost_LIBRARY_DIRS})


file(GLOB SOURCES *.cpp)


add_definitions(-std=c++14 -g -O3 -DBOOST_CONSTEXPR=)



add_executable(wASTRAL ${SOURCES})
target_link_libraries(wASTRAL profiler ${Boost_LIBRARIES} pthread )





