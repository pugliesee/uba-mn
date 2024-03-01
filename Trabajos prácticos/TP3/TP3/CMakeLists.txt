cmake_minimum_required(VERSION 3.5)
project(TP3)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -std=gnu++17 -O3")

set(SOURCE_FILES main.cpp JACOBI.cpp DLU.cpp GAUSS-SEIDEL.cpp)

set(HEADER_FILES JACOBI.h GAUSS-SEIDEL.h DLU.h)

# include_directories(${CMAKE_SOURCE_DIR}/eigen)


message("Compilando:" ${SOURCE_FILES})

# Ejecutable para correr main.cpp

add_executable(run ${SOURCE_FILES})

# Target para correr valgrind

add_custom_target(CorrerValgrind
        valgrind --leak-check=full -v ./run
        DEPENDS run)