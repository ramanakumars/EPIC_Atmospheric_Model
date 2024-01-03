cmake_minimum_required(VERSION 3.16)

project(epic VERSION 5.23)

set(ENV{EPIC_PATH} "${CMAKE_SOURCE_DIR}")
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
    ${CMAKE_CURRENT_SOURCE_DIR}/src/mpi/*.c 
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

add_executable(epic ${SRC_FILES} ${SINGLE_FILES} ${CMAKE_CURRENT_SOURCE_DIR}/src/core/epic_main.c)

set(LIBS netcdf -lm cdisort)
set(INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries (epic ${LIBS})
target_include_directories(epic PUBLIC ${INCLUDES})

add_executable(change ${SRC_FILES} ${SINGLE_FILES} ${CMAKE_CURRENT_SOURCE_DIR}/src/core/epic_change.c)
target_link_libraries (change ${LIBS})
target_include_directories(change PUBLIC ${INCLUDES})

add_executable(initial ${SRC_FILES} ${SINGLE_FILES} ${CMAKE_CURRENT_SOURCE_DIR}/src/core/epic_initial.c)
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

add_executable(mpi_epic ${SRC_FILES} ${MPI_FILES} ${CMAKE_CURRENT_SOURCE_DIR}/src/core/epic_main.c)
target_compile_definitions(mpi_epic PRIVATE EPIC_MPI=1)
target_link_libraries(mpi_epic MPI::MPI_C libmpg ${LIBS})
target_include_directories(mpi_epic PUBLIC ${INCLUDES} ${MPG_INCLUDES})
install(TARGETS mpi_epic DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/bin/)