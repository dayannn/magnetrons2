cmake_minimum_required(VERSION 3.5)

project(magnetrons LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fopenmp")


include_directories(
    C:/Qt/Tools/mingw810_64/lib/gcc/x86_64-w64-mingw32/8.1.0/include
)

add_executable(magnetrons main.cpp)
