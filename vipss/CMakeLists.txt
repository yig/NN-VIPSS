cmake_minimum_required(VERSION 3.14...3.30)
project(nnvipss)
set (CMAKE_CXX_STANDARD 20)
find_package(OpenMP REQUIRED)

## Load packages
include(FetchContent)
set(FETCHCONTENT_QUIET FALSE)

FetchContent_Declare(
    eigen
    # URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.zip
    # URL_HASH SHA256=1ccaabbfe870f60af3d6a519c53e09f3dcf630207321dffa553564a8e75c4fc8
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen
    GIT_TAG 3.3.9
    ## We want a newer tag to avoid Eigen's testing polluting our testing.
    # GIT_TAG 99ffad1971fb03de65f36e5401add8f872550e5f
    ## We can't use `GIT_SHALLOW TRUE` with an unnamed tag.
    GIT_SHALLOW FALSE
    GIT_PROGRESS TRUE
    PATCH_COMMAND  ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/patch/Eigen_99ffad197_CMakeLists.txt CMakeLists.txt
    )
FetchContent_MakeAvailable( eigen )
## We need C++20 features, but std::result was removed from the standard.
target_compile_definitions(eigen INTERFACE EIGEN_HAS_STD_RESULT_OF=0)

FetchContent_Declare(
    json
    GIT_REPOSITORY https://github.com/nlohmann/json.git
    GIT_TAG v3.12.0
    GIT_SHALLOW TRUE
    GIT_PROGRESS TRUE
    )
FetchContent_MakeAvailable( json )

FetchContent_Declare(
    NLopt
    GIT_REPOSITORY https://github.com/stevengj/nlopt
    GIT_TAG v2.10.0
    GIT_SHALLOW TRUE
    GIT_PROGRESS TRUE
    )
set(NLOPT_PYTHON OFF CACHE BOOL "" FORCE)
set(NLOPT_OCTAVE OFF CACHE BOOL "" FORCE)
set(NLOPT_MATLAB OFF CACHE BOOL "" FORCE)
set(NLOPT_GUILE OFF CACHE BOOL "" FORCE)
set(NLOPT_JAVA OFF CACHE BOOL "" FORCE)
set(NLOPT_SWIG OFF CACHE BOOL "" FORCE)
set(NLOPT_LUKSAN OFF CACHE BOOL "" FORCE)
set(NLOPT_TESTS OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable( NLopt )

FetchContent_Declare(
    PicoTree
    GIT_REPOSITORY https://github.com/Jaybro/pico_tree
    GIT_TAG v0.8.3
    GIT_SHALLOW TRUE
    GIT_PROGRESS TRUE
    ## A patch to skip pycotree
    PATCH_COMMAND  ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/patch/PicoTree_CMakeLists.txt src/CMakeLists.txt
    )
set(SKBUILD ON CACHE BOOL "" FORCE)
# FetchContent_Populate( PicoTree )
# file(WRITE "${PicoTree_SOURCE_DIR}/src/CMakeLists.txt" "add_subdirectory(pico_tree)")
FetchContent_MakeAvailable( PicoTree )

FetchContent_Declare(
    CLI11
    GIT_REPOSITORY https://github.com/CLIUtils/CLI11
    GIT_TAG v2.5.0
    GIT_SHALLOW TRUE
    GIT_PROGRESS TRUE
    )
FetchContent_MakeAvailable( CLI11 )

FetchContent_Declare(
    Armadillo
    GIT_REPOSITORY https://gitlab.com/conradsnicta/armadillo-code
    GIT_TAG 14.6.1
    GIT_SHALLOW TRUE
    GIT_PROGRESS TRUE
    )
FetchContent_MakeAvailable( Armadillo )

FetchContent_Declare(
    OpenBLAS
    GIT_REPOSITORY https://github.com/OpenMathLib/OpenBLAS
    GIT_TAG v0.3.30
    GIT_SHALLOW TRUE
    GIT_PROGRESS OpenBLAS
    )
FetchContent_MakeAvailable( OpenBLAS )

# find_library(ATOMIC_OPS_LIB atomic_ops REQUIRED)
# aux_source_directory(./adgrid ADGRID_LIST)
# aux_source_directory(./adgrid/nanothread NANOTHREAD_LIST)
# aux_source_directory(./adgrid/mshio MSHIO_LIST)
# aux_source_directory(./adgrid/predicates PREDICATES_LIST)
file(GLOB ADGRID_LIST ./adgrid/*.c*)
file(GLOB NANOTHREAD_LIST ./adgrid/nanothread/*.c*)
file(GLOB MSHIO_LIST ./adgrid/mshio/*.c*)
file(GLOB PREDICATES_LIST ./adgrid/predicates/predicates.c)
# set(PREDICATES_LIST ./adgrid/predicates/predicates.c)

# find_library(ATOMIC_LIB atomic_ops REQUIRED)
# message(STATUS "MY_VAR is: ${ATOMIC_LIB}")
# message(STATUS  ${ATOMIC_LIB})

## The adgrid library
add_library(adgrid ${ADGRID_LIST}  ${NANOTHREAD_LIST} ${MSHIO_LIST} ${PREDICATES_LIST} )
target_include_directories(adgrid PUBLIC ./adgrid/include )
target_include_directories(adgrid PUBLIC ./adgrid/include/mshio )
# target_link_libraries(adgrid  nlohmann_json::nlohmann_json atomic)
target_link_libraries(adgrid PUBLIC pico_tree nlopt nlohmann_json::nlohmann_json Eigen3::Eigen CLI11 )

## The local_vipss binary
aux_source_directory(./src SRC_LIST)
aux_source_directory(./src/surfacer SURFACER_LIST)
add_library(local_vipss ${SRC_LIST}  ${SURFACER_LIST})
target_include_directories(local_vipss PUBLIC ./src )
# include_directories(${NLOPT_INCLUDE_DIRS} ${ARMADILLO_INCLUDE_DIRS} ./src/surfacer)
target_include_directories(local_vipss PUBLIC ./src/surfacer )
target_link_libraries(local_vipss PUBLIC armadillo adgrid OpenMP::OpenMP_CXX )

# LINK_DIRECTORIES(${ARMADILLO_LIB_DIRS} ${NLOPT_LIB_DIR} )
# target_link_libraries(local_vipss ${ARMADILLO_LIB} ${NLOPT_LIBRARIES} OpenMP::OpenMP_CXX adgrid )

add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} PUBLIC local_vipss )
