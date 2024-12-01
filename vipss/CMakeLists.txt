
cmake_minimum_required(VERSION 3.10)
project(vipss)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++2a -O3")

find_package (Eigen3 REQUIRED)
find_package(NLopt REQUIRED)
find_package(PicoTree REQUIRED)
find_package(CLI11)
# set(NLopt_LIB_DIR "")
# set(NLopt_LIB ${NLOPT_LIBRARIES})
# add_definitions(-w)

find_package(OpenMP REQUIRED)
find_package(Armadillo REQUIRED)
# set(ARMADILLO_LIB_DIRS "")
set(ARMADILLO_LIB ${ARMADILLO_LIBRARIES})
include_directories( ${EIGEN3_INCLUDE_DIR} )

# set(PicoTree_INCLUDE_DIR "/usr/local/include/pico_tree")
# include(CMakePrintHelpers)
# cmake_print_variables(${EIGEN3_INCLUDE_DIR})
# cmake_print_variables(${QHULL_INCLUDE_DIR})

# include_directories(${PicoTree_INCLUDE_DIR})

include_directories(${NLOPT_INCLUDE_DIRS} ${ARMADILLO_INCLUDE_DIRS} ./src/surfacer)
aux_source_directory(. MAIN)
aux_source_directory(./src SRC_LIST)
aux_source_directory(./src/surfacer SURFACER_LIST)



LINK_DIRECTORIES(${ARMADILLO_LIB_DIRS} ${NLOPT_LIB_DIR} )
add_executable(${PROJECT_NAME} ${SRC_LIST} ${MAIN} ${SURFACER_LIST})

target_link_libraries(${PROJECT_NAME} ${ARMADILLO_LIB} ${NLOPT_LIBRARIES} 
            OpenMP::OpenMP_CXX)
