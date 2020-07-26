# beam-solver
Toy beam solver project for the numerical methods class I'm teaching. Includes a very basic benchmark suite.

### Prerequisites
- C++ compiler supporting select features of C++20 (tested using gcc 10.1 and MSVC 19.26)
- On Linux, TBB is required for parallel std::algorithms (solve.cpp only)

### Compilation
To compile this project on Linux follow the usual steps:
	
	cd beam-solver
	mkdir build
	cd build
	cmake ..
	make
	
On Windows you can simply point Visual Studio to the folder containing CMakeLists.txt and work with the project as you would with any other.

### Explanation
The general idea is for students to try their hand at writing a simple parallel finite element solver. This repository contains basic utilities for both FEM (`MesLib.h`) and parallel programming (`ParLib.hpp`). Note that these utilities are only meant to illustrate some basic concepts, they are far from anything one would find in real world production code. `ParLib` is really just a C-style interface for the most elementary STL multithreading functionality.

To benchmark your own solver you need to substitute the `solve` function with your implementation. The function signature is `void solve(double*, const size_t)`. The first argument specifies the address where the solution vector is to be written, and the second one defines the number of threads which are to be used for the computations. All other parameters (geometry, load, BCs, etc.) are set using global variables.

The solver computes the solution of a linear elasticity problem in 2D. The domain is a rectangular one, discretized using a cartesian grid of unit square 1st order quadrilateral finite elements (basically the simplest case imaginable). For additional details, as well as some documentation of `MesLib` and `ParLib` see [this page](http://ccfd.github.io/courses/metnum_lab2.html) and [this other page](http://ccfd.github.io/courses/metnum_lab4.html) (unfortunately only in Polish at the time of writing).
