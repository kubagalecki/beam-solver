cmake_minimum_required(VERSION 3.1)

project(konkurs)

SET(CXX_STANDARD_REQUIRED ON)
SET(CMAKE_CXX_STANDARD 20)
SET(CMAKE_BUILD_TYPE RELEASE)

if (UNIX)
    SET(THREADS_PREFER_PTHREAD_FLAG ON)
    find_package(Threads REQUIRED)
    find_package(TBB REQUIRED)
endif()

add_executable(konkurs
        solve.cpp
        benchmark_beam_solution.cpp
        MesLib.cpp
        ParLib.cpp ThreadPool.cpp)

if (UNIX)
    target_link_libraries(konkurs Threads::Threads)
    target_link_libraries(konkurs ${TBB_IMPORTED_TARGETS})
endif()