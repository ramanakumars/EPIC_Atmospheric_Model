cmake_minimum_required(VERSION 3.10)
cmake_policy(SET CMP0060 NEW)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

project(epic VERSION 5.30)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DEPIC_PATH=\\\"$ENV{EPIC_PATH}\\\" -W")
Message("CMAKE_C_FLAGS: ")
Message(${CMAKE_C_FLAGS})

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

add_subdirectory(src/mpi/mpg)
add_subdirectory(src/rt/cdisort)

set(LIBS -lm cdisort)

add_executable(epic ${SRC_FILES} ${SINGLE_FILES} ${CMAKE_CURRENT_SOURCE_DIR}/src/core/epic_main.c)
add_executable(change ${SRC_FILES} ${SINGLE_FILES} ${CMAKE_CURRENT_SOURCE_DIR}/src/core/epic_change.c)
add_executable(initial ${SRC_FILES} ${SINGLE_FILES} ${CMAKE_CURRENT_SOURCE_DIR}/src/core/epic_initial.c)
add_executable(mpi_epic ${SRC_FILES} ${MPI_FILES} ${CMAKE_CURRENT_SOURCE_DIR}/src/core/epic_main.c)

find_library(NETCDF_LIB NAMES netcdf)
find_path(NETCDF_INCLUDES netcdf.h)
SET(LIBS ${LIBS} ${NETCDF_LIB})
set(INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/include ${NETCDF_INCLUDES})

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


install(TARGETS epic 
    LIBRARY DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/lib/
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/bin/)
install(TARGETS initial 
    LIBRARY DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/lib/
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/bin/)
install(TARGETS change
    LIBRARY DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/lib/
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/bin/)
find_package(MPI REQUIRED)

target_compile_definitions(mpi_epic PRIVATE EPIC_MPI=1)
target_link_libraries(mpi_epic MPI::MPI_C mpg ${LIBS})
target_include_directories(mpi_epic PUBLIC ${INCLUDES} ${MPG_INCLUDES})
install(TARGETS mpi_epic 
    LIBRARY DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/lib/
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/bin/)
