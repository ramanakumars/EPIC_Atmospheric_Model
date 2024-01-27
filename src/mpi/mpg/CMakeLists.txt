file(GLOB SOURCES "*.c")

find_package(MPI REQUIRED)

add_library(mpg STATIC ${SOURCES})

target_include_directories(mpg PUBLIC include/)
target_link_libraries(mpg MPI::MPI_C)