cmake_minimum_required (VERSION 2.6)

set (LEMONADE_DIR "/scratch/localuser/lemonade/lemonadeInstall/")
#set (LEMONADE_DIR "/scratch/localuser/lemonade/lemonadeToniMueller/install/")
#set (LEMONADE_DIR "/home/s7146169/softwareSCS5/LeMonADE/LeMonADE-install/")

SET (LEMONADE_INCLUDE_DIR ${LEMONADE_DIR}include/)
SET (LEMONADE_LIBRARY_DIR ${LEMONADE_DIR}lib/)

message( ${LEMONADE_INCLUDE_DIR} )

SET (CMAKE_RUNTIME_OUTPUT_DIRECTORY "./bin/")

include_directories("../updater" "../analyzer" "../features")

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

SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_DEBUG} -O3 -std=c++11")
SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS_DEBUG} -O3 ")

## ###############  System Creators ############# ##

add_executable(createFixedChainInSlit createChainInSlit.cpp)
target_link_libraries(createFixedChainInSlit LeMonADE ${Boost_LIBRARIES})

## ###############  Analyzers ############# ##

#add_executable(evaluateSelectCloseMonomers evaluateSelectCloseMonomers.cpp)
#target_link_libraries(evaluateSelectCloseMonomers LeMonADE ${Boost_LIBRARIES})

## ###############  Simulators ############# ##

add_executable(SimualtorChainInSlitForce simulatorSlitChain.cpp)
target_link_libraries(SimualtorChainInSlitForce LeMonADE ${Boost_LIBRARIES})

## ###############  Modifiers ############# ##

#add_executable(setUpCUDANNInteraction setUpCUDANNInteractions.cpp)
#target_link_libraries(setUpCUDANNInteraction LeMonADE ${Boost_LIBRARIES})
