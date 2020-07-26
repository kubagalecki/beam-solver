#include "BenchSuite.hpp"
#include "solve.h"

std::uint_fast16_t mx;
std::uint_fast16_t my;
bool*              fix;
double*            F;
double*            d_global;

bool validate(const double* d)
{
    std::vector< double > res((mx + 1) * (my + 1) * 2);

    for (size_t elx = 0; elx < mx; ++elx)
    {
        for (size_t ely = 0; ely < my; ++ely)
        {
            for (size_t i = 0; i < 8; ++i)
            {
                const size_t dofi = DOF(elx, ely, i);

                if (fix[dofi])
                    continue;

                for (size_t j = 0; j < 8; ++j)
                {
                    const size_t dofj = DOF(elx, ely, j);
                    res[dofi] += K[i][j] * d[dofj];
                }
            }
        }
    }

    for (size_t i = 0; i < res.size(); ++i)
        res[i] -= F[i];

    constexpr double tol = 1.01e-6;
    return sqrt(std::inner_product(res.cbegin(), res.cend(), res.cbegin(), 0.)) <= tol;
}

void allocate_globals()
{
    const size_t problem_size = (mx + 1) * (my + 1) * 2;
    F                         = new double[problem_size];
    d_global                  = new double[problem_size];
    fix                       = new bool[problem_size];

    for (size_t i = 0; i < problem_size; ++i)
    {
        F[i]   = 0.;
        fix[i] = false;
    }
}

void deallocate_globals()
{
    delete[] F;
    delete[] d_global;
    delete[] fix;
}

int main()
{
    // Define beam sizes
    const std::vector< std::pair< size_t, size_t > > size_span = {{100, 100}, {20, 500}, {500, 20}};

    // Define load configurations
    const std::vector< std::function< void() > > case_span = {[]() {
                                                                  for (size_t i = 0; i <= my; ++i)
                                                                  {
                                                                      fix[P(0, i, 0)] = true;
                                                                      fix[P(0, i, 1)] = true;
                                                                  }

                                                                  F[P(mx, 0, 1)] = -0.01;
                                                              },

                                                              []() {
                                                                  for (size_t i = 0; i <= mx; ++i)
                                                                  {
                                                                      fix[P(i, my, 0)] = true;
                                                                      fix[P(i, my, 1)] = true;
                                                                  }

                                                                  for (size_t i = 0; i <= mx; ++i)
                                                                      F[P(i, 0, 1)] = -0.001;
                                                              },

                                                              []() {
                                                                  for (size_t i = 0; i <= my; ++i)
                                                                  {
                                                                      fix[P(0, i, 0)]  = true;
                                                                      fix[P(0, i, 1)]  = true;
                                                                      fix[P(mx, i, 0)] = true;
                                                                      fix[P(mx, i, 1)] = true;
                                                                  }

                                                                  for (size_t i = 1; i < mx; ++i)
                                                                      F[P(i, 0, 1)] = -0.001;
                                                              }};

    // Define thread allotments
    const std::vector< size_t > thread_span = {1, 2, 4};

    BenchSuite benchmarks(Verbosity::medium);

    // Push all size-load-thread_no combinations
    for (const auto& size : size_span)
        for (const auto& load : case_span)
            for (const auto& number_of_threads : thread_span)
                benchmarks.push_benchmark(
                    [mx_ = size.first, my_ = size.second, l = load]() {
                        mx = mx_;
                        my = my_;
                        allocate_globals();
                        l();
                    },
                    deallocate_globals,
                    [d = d_global](double* a1, const size_t a2) {
                        solve(a1, a2);
                        return a1;
                    },
                    validate,
                    std::reference_wrapper(d_global),
                    number_of_threads);

    const double total_time = benchmarks.run();
    printf("\n\nCumulative median execution time (in seconds):\t%f\n\n", total_time);
}
