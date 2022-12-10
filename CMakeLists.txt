cmake_minimum_required(VERSION 3.0.0)
project(magma_fields VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_EXPORT_COMPILE_COMMANDS true)
set(CMAKE_CXX_FLAGS "-Ofast -flto -ffast-math -Wall -Wextra -Wno-unused-parameter -pedantic")
set(CMAKE_EXE_LINKER_FLAGS "-Ofast -flto -ffast-math")

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR})
find_package(GMP REQUIRED)

include_directories(magma_fields PUBLIC ${CMAKE_SOURCE_DIR} ${GMP_INCLUDE_DIRS})

add_executable(magma_fields main.cpp)
target_link_libraries(magma_fields gmp gmpxx)
# set_target_properties(magma_fields PROPERTIES COMPILE_FLAGS "-DNDEBUG")
