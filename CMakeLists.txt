cmake_minimum_required(VERSION 3.30)
project(Sem6Laba3_FDTD)

set(CMAKE_CXX_STANDARD 20)

add_executable(Sem6Laba3_FDTD
        FDTD1D.cpp
        PMLCoefficients.cpp
        Monitor.cpp
        PMLAnalysis.cpp
        Material.cpp
        main.cpp)
