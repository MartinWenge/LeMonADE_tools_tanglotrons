cmake_minimum_required (VERSION 2.6)

SET (LEMONADE_INCLUDE_DIR "/scratch/localuser/lemonade/lemonadeInstall/include/")
SET (LEMONADE_LIBRARY_DIR "/scratch/localuser/lemonade/lemonadeInstall/lib")

include_directories("../updater" "../analyzer" "../feature")

if (NOT DEFINED LEMONADE_INCLUDE_DIR)
message("LEMONADE_INCLUDE_DIR is not provided. If build fails, use -DLEMONADE_INCLUDE_DIR=/path/to/LeMonADE/headers/ or install to default location")
endif()

if (NOT DEFINED LEMONADE_LIBRARY_DIR)
message("LEMONADE_LIBRARY_DIR is not provided. If build fails, use -DLEMONADE_LIBRARY_DIR=/path/to/LeMonADE/lib/ or install to default location")
endif()

find_package( Boost REQUIRED COMPONENTS program_options)
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

include_directories (${LEMONADE_INCLUDE_DIR})
link_directories (${LEMONADE_LIBRARY_DIR})

SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_DEBUG} -O2 -std=c++11")
SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS_DEBUG} -O2 ")

## ###############  test executable  ############# ##
add_executable(testTanglotron test_main.cpp test_createChainInSlit.cpp test_analyzerForce.cpp)
target_link_libraries(testTanglotron LeMonADE ${Boost_LIBRARIES})

