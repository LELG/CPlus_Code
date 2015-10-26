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
# and MPI

find_package(Boost COMPONENTS program_options filesystem system REQUIRED)
message(status "** Boost Include: ${Boost_INCLUDE_DIR}")
message(status "** Boost Libraries Dir: ${Boost_LIBRARY_DIRS}")
message(status "** Boost Libraries: ${Boost_LIBRARIES}")

find_package(OpenMP REQUIRED)

find_package(MPI REQUIRED)
message(status "** MPI Include: ${MPI_INCLUDE_PATH}")
message(status "** MPI Libraries: ${MPI_LIBRARIES}")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")
find_package(GSL REQUIRED)
message(status "** GSL Include: ${GSL_INCLUDE_DIR}")
message(status "** GSL Libraries: ${GSL_LIBRARIES}")

add_executable (testParser test/test_parser.cpp test/test_main.cpp test/catch.hpp parser.cpp)
add_executable (tumourSim V1_Main_Tumour_Evoulution.cpp parser.cpp Random.cpp)
add_executable (tumourSimTreatment TEM_V1.cpp parser.cpp)

include_directories("${Boost_INCLUDE_DIR}")
include_directories("${GSL_INCLUDE_DIR}")
include_directories("${MPI_INCLUDE_PATH}")
include_directories("${LOCAL_INCLUDES}")
include_directories("${PROJECT_SOURCE_DIR}/test")
link_directories("${Boost_LIBRARY_DIRS}")

set_target_properties(tumourSim PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS})

target_link_libraries(testParser ${Boost_LIBRARIES})
target_link_libraries(tumourSim ${Boost_LIBRARIES})
target_link_libraries(tumourSim ${GSL_LIBRARIES})
target_link_libraries(tumourSimTreatment ${GSL_LIBRARIES} ${MPI_LIBRARIES} ${Boost_LIBRARIES})

if(MPI_COMPILE_FLAGS)
  set_target_properties(tumourSimTreatment PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(tumourSimTreatment PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()
