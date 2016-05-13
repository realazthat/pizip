cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)

project(pizip)

#relative link directories are to be turned into full paths
cmake_policy(SET CMP0015 NEW)

include_directories(./libs/CTPL/CTPL ./libs/maginatics-threadpool/threadpool/src)

################################################################################
add_executable(gmp-chudnovsky
    gmp-chudnovsky.c
    )

target_link_libraries(gmp-chudnovsky
    gmp)

set_property(TARGET gmp-chudnovsky PROPERTY C_STANDARD 99)
set_property(TARGET gmp-chudnovsky PROPERTY C_STANDARD_REQUIRED ON)
################################################################################


add_executable(pizip
    main.cpp
)


target_link_libraries(pizip )


set_property(TARGET pizip PROPERTY CXX_STANDARD 11)
set_property(TARGET pizip PROPERTY CXX_STANDARD_REQUIRED ON)


################################################################################
