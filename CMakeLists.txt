cmake_minimum_required (VERSION 2.8)

project (TestParser)

if(DEFINED ENV{PBS_SERVER})
    message("PBS_SERVER:" $ENV{PBS_SERVER})
    if (ENV{PBS_SERVER} strequal "bioinf-head.petermac.org.au")
	    message(status "Building on the PMCC Cluster")
	    message(status "If you have problems, try running cmake with:")
	    message(status "    CC=gcc CXX=g++ cmake /path/to/source")
	    message(status "Setting Boost root directory ...")
	    set(BOOST_ROOT "/usr/local/cluster/all_arch/boost/boost_1_54_0/")
    endif()
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11")

# Search for required libraries - Boost, OpenMP, GSL

find_package(Boost COMPONENTS program_options REQUIRED)
message(status "** Boost Include: ${Boost_INCLUDE_DIR}")
message(status "** Boost Libraries Dir: ${Boost_LIBRARY_DIRS}")
message(status "** Boost Libraries: ${Boost_LIBRARIES}")

find_package(OpenMP REQUIRED)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")
find_package(GSL REQUIRED)
message(status "** GSL Include: ${GSL_INCLUDE_DIRS}")
message(status "** GSL Libraries: ${GSL_LIBRARIES}")

add_executable (testParser test_parser.cpp params.cpp)
add_executable (tumourSim V1_Main_Tumour_Evoulution.cpp params.cpp)

include_directories("${Boost_INCLUDE_DIR}")
include_directories("${GSL_INCLUDE_DIRS}")
link_directories("${Boost_LIBRARY_DIRS}")

set_target_properties(tumourSim PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS})

target_link_libraries(testParser ${Boost_LIBRARIES})
target_link_libraries(tumourSim ${Boost_LIBRARIES})
target_link_libraries(tumourSim ${GSL_LIBRARIES})

install (TARGETS testParser DESTINATION bin)
