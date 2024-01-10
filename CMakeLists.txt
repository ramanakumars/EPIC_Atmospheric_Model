cmake_minimum_required(VERSION 3.16)
cmake_policy(SET CMP0060 NEW)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

project(epic VERSION 5.23)

set(ENV{EPIC_PRECISION} 8)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DEPIC_PRECISION=$ENV{EPIC_PRECISION} -DEPIC_PATH=\\\"$ENV{EPIC_PATH}\\\"")

file(GLOB_RECURSE SRC_FILES 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/core/*.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/clouds/*.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/chemistry/*.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/subgrid/*.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/rt/rt_heating.c 
)

list (REMOVE_ITEM SRC_FILES "${CMAKE_SOURCE_DIR}/src/core/epic_initial.c")
list (REMOVE_ITEM SRC_FILES "${CMAKE_SOURCE_DIR}/src/core/epic_change.c")
list (REMOVE_ITEM SRC_FILES "${CMAKE_SOURCE_DIR}/src/core/epic_main.c")

file(GLOB_RECURSE SINGLE_FILES 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/single/*.c 
)

file(GLOB_RECURSE MPI_FILES 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/mpi/m_epic_functions.c
)

add_custom_target(
    buildcdisort ALL
    COMMAND make -f Makefile_cdisort install
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/src/rt/cdisort
)
add_library(cdisort::cdisort STATIC IMPORTED)
set_target_properties(cdisort::cdisort
    PROPERTIES
  IMPORTED_LOCATION ${CMAKE_CURRENT_SOURCE_DIR}/src/rt/cdisort/cdisort.o 
  INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR}/src/rt/cdisort
)
add_library(cdisort::locate STATIC IMPORTED)
set_target_properties(cdisort::locate
    PROPERTIES
  IMPORTED_LOCATION ${CMAKE_CURRENT_SOURCE_DIR}/src/rt/cdisort/locate.o 
  INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR}/src/rt/cdisort
)

add_library(cdisort INTERFACE IMPORTED)
set_property(
    TARGET cdisort PROPERTY 
    INTERFACE_LINK_LIBRARIES cdisort::cdisort cdisort::locate
)
add_dependencies(cdisort buildcdisort)

set(LIBS -lm cdisort)
set(INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_executable(epic ${SRC_FILES} ${SINGLE_FILES} ${CMAKE_CURRENT_SOURCE_DIR}/src/core/epic_main.c)
add_executable(change ${SRC_FILES} ${SINGLE_FILES} ${CMAKE_CURRENT_SOURCE_DIR}/src/core/epic_change.c)
add_executable(initial ${SRC_FILES} ${SINGLE_FILES} ${CMAKE_CURRENT_SOURCE_DIR}/src/core/epic_initial.c)
add_executable(mpi_epic ${SRC_FILES} ${MPI_FILES} ${CMAKE_CURRENT_SOURCE_DIR}/src/core/epic_main.c)

find_library(NETCDF_LIB NAMES netcdf)
if(NOT NETCDF_LIB)
    include(ExternalProject)
    ExternalProject_Add(netCDF 
        URL https://github.com/Unidata/netcdf-c/archive/refs/tags/v4.9.2.zip
        PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/netcdf/
        CMAKE_ARGS 
            -DENABLE_DAP=0
            -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_SOURCE_DIR}/netcdf/
    )
    
    add_dependencies(epic netCDF)
    add_dependencies(mpi_epic netCDF)
    add_dependencies(change netCDF)
    add_dependencies(initial netCDF)

    SET(INCLUDES ${INCLUDES} ${CMAKE_CURRENT_SOURCE_DIR}/netcdf/include)
    target_link_libraries(epic ${CMAKE_CURRENT_SOURCE_DIR}/netcdf/lib/libnetcdf.so)
    target_link_libraries(mpi_epic ${CMAKE_CURRENT_SOURCE_DIR}/netcdf/lib/libnetcdf.so)
    target_link_libraries(initial ${CMAKE_CURRENT_SOURCE_DIR}/netcdf/lib/libnetcdf.so)
    target_link_libraries(change ${CMAKE_CURRENT_SOURCE_DIR}/netcdf/lib/libnetcdf.so)
else()
    SET(LIBS ${LIBS} netcdf)
endif()

message("LIST OF INCLUDES: ")
message(${INCLUDES})

message("LIST OF LIBS: ")
message(${LIBS})

target_link_libraries (epic ${LIBS})
target_include_directories(epic PUBLIC ${INCLUDES})

target_link_libraries (change ${LIBS})
target_include_directories(change PUBLIC ${INCLUDES})

target_link_libraries (initial ${LIBS})
target_include_directories(initial PUBLIC ${INCLUDES})


install(TARGETS epic DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/bin/)
install(TARGETS initial DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/bin/)
install(TARGETS change DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/bin/)

SET(MPG_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/src/mpi/mpg/mpg.h;${INCLUDES}")
find_package(MPI REQUIRED)
add_custom_target(
    buildmpg ALL
    COMMAND make -f Makefile.openmpi install
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/src/mpi/mpg/
)
add_library(libmpg STATIC IMPORTED)
set_target_properties(libmpg
    PROPERTIES
  IMPORTED_LOCATION ${CMAKE_CURRENT_SOURCE_DIR}/src/mpi/mpg/libmpg.a
  INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_CURRENT_SOURCE_DIR}/src/mpi/mpg/;${INCLUDES}"
)
add_dependencies(libmpg buildmpg)

# file(GLOB_RECURSE MPG_FILES 
#     ${CMAKE_CURRENT_SOURCE_DIR}/src/mpi/mpg/*.c
# )
# add_executable(libmpg.o ${MPG_FILES})
# target_include_directories(libmpg.o PUBLIC ${INCLUDES} ${CMAKE_CURRENT_SOURCE_DIR}/src/mpi/mpg/mpg.h)
# target_link_libraries(libmpg.o MPI::MPI_C -lm)
# install(TARGETS libmpg.o DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/lib)

# add_library(libmpg SHARED IMPORTED)
# set_target_properties(libmpg
#     PROPERTIES
#   IMPORTED_LOCATION ${CMAKE_CURRENT_SOURCE_DIR}/lib/libmpg.o
#   INTERFACE_INCLUDE_DIRECTORIES "${INCLUDES};${CMAKE_CURRENT_SOURCE_DIR}/src/mpi/mpg/mpg.h"
# )
# add_dependencies(libmpg libmpg.o)

target_compile_definitions(mpi_epic PRIVATE EPIC_MPI=1)
target_link_libraries(mpi_epic MPI::MPI_C libmpg ${LIBS})
target_include_directories(mpi_epic PUBLIC ${INCLUDES} ${MPG_INCLUDES})
install(TARGETS mpi_epic DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/bin/)