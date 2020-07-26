#ifndef KONKURS_SOLVE_H
#define KONKURS_SOLVE_H

#include <cmath>
#include <cstddef>
#include <cstdint>

#include <algorithm>
#include <array>
#include <execution>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <numeric>
#include <utility>
#include <vector>

#include "MesLib.h"
#include "ParLib.hpp"
#include "ThreadPool.hpp"

void solve(double*, const std::size_t);

#endif // KONKURS_SOLVE_H
