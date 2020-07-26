#ifndef KONKURS_BENCHSUITE_HPP
#define KONKURS_BENCHSUITE_HPP

#include <algorithm>
#include <chrono>
#include <cmath>
#include <functional>
#include <iostream>
#include <memory>
#include <numeric>
#include <optional>
#include <queue>
#include <utility>
#include <vector>

enum class Verbosity
{
    none,
    low,
    medium,
    high,
};

class BenchSuite
{
public:
    BenchSuite(Verbosity _v      = Verbosity::low,
               size_t    minRuns = 10,
               size_t    maxRuns = 50,
               double    varCoef = 0.05)
        : verbosity(_v), min_runs(minRuns), max_runs(maxRuns), var_coef(varCoef)
    {}

    BenchSuite(BenchSuite&&) = default;
    BenchSuite& operator=(BenchSuite&&) = default;
    BenchSuite(const BenchSuite&)       = delete;
    BenchSuite& operator=(const BenchSuite&) = delete;
    ~BenchSuite()                            = default;

    template < typename S, typename T, typename F, typename V, typename... Args >
    void push_benchmark(S&& setup, T&& teardown, F&& fun, V&& val, Args&&... args)
    {
        bench_queue.emplace([this,
                             s     = std::forward< S >(setup),
                             t     = std::forward< T >(teardown),
                             f     = std::forward< F >(fun),
                             v     = std::forward< V >(val),
                             ... a = std::forward< Args >(args)]() {
            s();
            bench(std::move(f), std::move(v), std::move(a)...);
            t();
        });
    }

    double run()
    {
        while (!bench_queue.empty())
        {
            bench_queue.front()();
            bench_queue.pop();
        }

        double ret;

        if (std::any_of(bench_times.cbegin(),
                        bench_times.cend(),
                        [](const std::optional< double >& t) { return !t; }))
            ret = std::numeric_limits< double >::infinity();
        else
            ret = std::accumulate(bench_times.cbegin(),
                                  bench_times.cend(),
                                  0.,
                                  [](const auto& a, const auto& b) { return a + (*b); });

        if (verbosity != Verbosity::none)
        {
            puts("\n=====================================");

            puts("SUITE SUMMARY:\n");
            std::cout << "Benchmarks performed:     " << bench_n << '\n';
            std::cout << "Benchmarks passed:        "
                      << std::count_if(bench_times.cbegin(),
                                       bench_times.cend(),
                                       [](const auto& b) { return b.has_value(); })
                      << '\n';
            if (ret != std::numeric_limits< double >::infinity())
                std::cout << "Median cumulative time:   " << ret << '\n';

            puts("=====================================");
        }

        return ret;
    }

private:
    template < typename Fun, typename Valid, typename... Args >
    void bench(const Fun& fun, const Valid& val, const Args&... arg)
    {
        ++bench_n;

        bool   keep_running = true;
        size_t validated    = 0;
        size_t n_runs       = 0;

        std::vector< double > exec_time;
        exec_time.reserve(max_runs);

        const auto mean_comp = [&exec_time]() {
            return std::accumulate(exec_time.cbegin(), exec_time.cend(), 0.) /
                   static_cast< double >(exec_time.size());
        };

        const auto std_comp = [&]() {
            const double mean = mean_comp();
            return sqrt(
                std::transform_reduce(exec_time.cbegin(),
                                      exec_time.cend(),
                                      0.,
                                      std::plus<>{},
                                      [&mean](const double& a) { return pow(a - mean, 2); }) /
                static_cast< double >(exec_time.size()));
        };

        while (keep_running && n_runs < max_runs)
        {
            ++n_runs;

            const auto t_start = std::chrono::steady_clock::now();
            const auto res     = fun(arg...);
            const auto t_end   = std::chrono::steady_clock::now();
            exec_time.push_back(std::chrono::duration< double >(t_end - t_start).count());

            if (verbosity == Verbosity::high)
            {
                std::cout << "Run " << n_runs << ":\t\t";
            }

            if (!val(res))
            {
                if (verbosity == Verbosity::high)
                    std::cout << "failed,\texecution time: " << exec_time.back() << '\n';

                exec_time.pop_back();
                continue;
            }
            ++validated;
            if (verbosity == Verbosity::high)
                std::cout << "success,\texecution time: " << exec_time.back() << '\n';

            if (n_runs >= min_runs)
            {
                const double mean = mean_comp();

                const double std_dev = std_comp();

                if (std_dev / mean <= var_coef)
                    keep_running = false;

                if (n_runs > max_runs)
                    keep_running = false;
            }
        }

        std::optional< double > ret;

        if (validated == 0 && (verbosity == Verbosity::medium || verbosity == Verbosity::high))
        {
            puts("--------------------------");
            std::cout << "Benchmark " << bench_n << " failed on all runs\n";
            puts("--------------------------");
        }
        else
        {
            // Return median
            std::sort(exec_time.begin(), exec_time.end());
            ret = exec_time[exec_time.size() / 2];

            if (verbosity == Verbosity::medium || verbosity == Verbosity::high)
            {
                puts("--------------------------");
                std::cout << "Benchmark " << bench_n << " summary:\n\n";
                std::cout << "Number of runs:               " << n_runs << '\n';
                std::cout << "Number of successful runs:    " << validated << '\n';
                std::cout << "Mean run time:                " << mean_comp() << '\n';
                std::cout << "Median run time:              " << ret.value() << '\n';
                std::cout << "Run time std. dev.:           " << std_comp() << '\n';
                puts("--------------------------");
            }

            if (validated != n_runs)
                ret.reset();
        }

        bench_times.push_back(ret);
    }

    const Verbosity verbosity = Verbosity::low;
    const size_t    min_runs;
    const size_t    max_runs;
    const double    var_coef;

    size_t bench_n = 0;

    std::queue< std::function< void() > >  bench_queue = {};
    std::vector< std::optional< double > > bench_times;
};

#endif // KONKURS_BENCHSUITE_HPP
