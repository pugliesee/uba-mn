cmake_minimum_required(VERSION 3.5)
project(2TP_Procesamiento_de_imágenes)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -std=gnu++17 -O3")

set(SOURCE_FILES main.cpp)

# include_directories(${CMAKE_SOURCE_DIR}/eigen)


message("Compilando:" ${SOURCE_FILES})

# Ejecutable para correr main.cpp

add_executable(run ${SOURCE_FILES})

# Target para correr valgrind

add_custom_target(CorrerValgrind
        valgrind --leak-check=full -v ./run
        DEPENDS run)