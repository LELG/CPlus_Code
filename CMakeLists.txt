cmake_minimum_required (VERSION 2.8)

project (TumourSim)

if(DEFINED ENV{PBS_SERVER})
    message("PBS_SERVER:" $ENV{PBS_SERVER})
    if ($ENV{PBS_SERVER} STREQUAL "bioinf-head.petermac.org.au")
	    message(status "Building on the PMCC Cluster")
	    message(status "Setting Boost root directory ...")
	    set(BOOST_ROOT "/usr/local/cluster/all_arch/boost/boost_1_54_0/")
    endif()
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11")

set(LOCAL_INCLUDES ${PROJECT_SOURCE_DIR}/include)

# Search for required libraries - Boost, OpenMP, GSL

find_package(Boost COMPONENTS program_options filesystem system REQUIRED)
message(status "** Boost Include: ${Boost_INCLUDE_DIR}")
message(status "** Boost Libraries Dir: ${Boost_LIBRARY_DIRS}")
message(status "** Boost Libraries: ${Boost_LIBRARIES}")

find_package(OpenMP REQUIRED)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")
find_package(GSL REQUIRED)
message(status "** GSL Include: ${GSL_INCLUDE_DIR}")
message(status "** GSL Libraries: ${GSL_LIBRARIES}")

add_executable (testParser test/test_parser.cpp test/test_main.cpp test/catch.hpp parser.cpp)
add_executable (tumourSim V1_Main_Tumour_Evoulution.cpp parser.cpp Random.cpp)

include_directories("${Boost_INCLUDE_DIR}")
include_directories("${GSL_INCLUDE_DIR}")
include_directories("${LOCAL_INCLUDES}")
include_directories("${PROJECT_SOURCE_DIR}/test")
link_directories("${Boost_LIBRARY_DIRS}")

set_target_properties(tumourSim PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS})

target_link_libraries(testParser ${Boost_LIBRARIES})
target_link_libraries(tumourSim ${Boost_LIBRARIES})
target_link_libraries(tumourSim ${GSL_LIBRARIES})
