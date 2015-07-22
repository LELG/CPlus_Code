cmake_minimum_required (VERSION 2.8)

project (TestParser)

if(DEFINED ENV{PBS_SERVER})
    message("PBS_SERVER:" $ENV{PBS_SERVER})
    if (ENV{PBS_SERVER} strequal "bioinf-head.petermac.org.au")
	    message(STATUS "Building on the PMCC Cluster")
	    message(STATUS "Setting CXX compiler ...")
	    SET(CMAKE_CXX_COMPILER "g++")
	    message(STATUS "CXX compiler is: ${CMAKE_CXX_COMPILER}")
	    message(STATUS "Setting Boost root directory ...")
	    set(BOOST_ROOT "/usr/local/cluster/all_arch/boost/boost_1_54_0/")
    endif()
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11")

find_package(Boost COMPONENTS program_options REQUIRED)
message(status "** Boost Include: ${Boost_INCLUDE_DIR}")
message(status "** Boost Libraries: ${Boost_LIBRARY_DIRS}")
message(status "** Boost Libraries: ${Boost_LIBRARIES}")

add_executable (testParser test_parser.cpp params.cpp)

include_directories("${PROJECT_SOURCE_DIR}")
include_directories("${Boost_INCLUDE_DIR}")
link_directories("${Boost_LIBRARY_DIRS}")
target_link_libraries(testParser ${Boost_LIBRARIES})

install (TARGETS testParser DESTINATION bin)
