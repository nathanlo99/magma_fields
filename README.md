# magma_fields

An implementation of Finite Fields in C++ as they (hopefully) appear in the Magma CAS

## Prerequisites

- CMake
- GMP (https://gmplib.org/)

## Build instructions

```
mkdir -p build
cd build
cmake ..
cd ../
cmake --build build
```

## Troubleshooting 

If you are running from Linux or compiling with certain compilers, you may get a message asking you to add the flag `-fconcepts`. If this happens, replace the lines in the `CMakeLists.txt` file with 
```
set(CMAKE_CXX_FLAGS "-fconcepts -Ofast -flto -ffast-math -Wall -Wextra -Wno-unused-parameter -pedantic")
set(CMAKE_EXE_LINKER_FLAGS "-fconcepts -Ofast -flto -ffast-math")
```

## To run the executable

```
build/magma_fields
```

## Resources:

**BCS97**  
Wieb Bosma, John Cannon, and Allan Steel.  
Lattices of Compatibly Embedded Finite Fields.  
J. Symbolic Comp., 24(3):351--369, 1997.   
