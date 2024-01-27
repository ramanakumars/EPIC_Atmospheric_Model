file(GLOB SOURCES "*.c")

add_library(cdisort STATIC ${SOURCES})

target_include_directories(cdisort PUBLIC .)
target_link_libraries(cdisort -lm)
