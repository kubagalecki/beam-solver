cmake_minimum_required(VERSION 3.1)

project(konkurs)

set(CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_BUILD_TYPE RELEASE)

if (UNIX)
    set(THREADS_PREFER_PTHREAD_FLAG ON)
    find_package(Threads REQUIRED)
    find_package(TBB REQUIRED)
endif()

include_directories(${CMAKE_SOURCE_DIR}/include/)

add_executable(benchmarks
        ${CMAKE_SOURCE_DIR}/src/solve.cpp
        ${CMAKE_SOURCE_DIR}/src/benchmark_beam_solution.cpp
        ${CMAKE_SOURCE_DIR}/src/MesLib.cpp
        ${CMAKE_SOURCE_DIR}/src/ParLib.cpp
        ${CMAKE_SOURCE_DIR}/src/ThreadPool.cpp
	)

if (UNIX)
    target_link_libraries(konkurs Threads::Threads)
    target_link_libraries(konkurs ${TBB_IMPORTED_TARGETS})
endif()
