#ifndef KONKURS_PARLIB_HPP
#define KONKURS_PARLIB_HPP

#include <thread>
#include <mutex>
#include <condition_variable>
#include <chrono>
#include <memory>
#include <vector>
#include <map>
#include <algorithm>
#include <iostream>

void sync();
unsigned int self_id();
unsigned int no_threads();
void init_mutex(const size_t&);
void mutex_lock(const size_t&);
void mutex_unlock(const size_t&);
void tic(const size_t&);
double toc(const size_t&);

namespace ParLib
{
    // Thread and mutex pools
    extern std::vector<std::thread>	                                                thread_pool;
    extern std::unique_ptr<std::mutex[]>	                                        mutex_pool;
    extern size_t                                                                   number_of_threads;
    extern std::vector<std::thread::id>                                             thread_ids;

    // Syncing stuff
    extern std::mutex                                                               sync_mutex;
    extern size_t                                                                   sync_counter;
    extern std::condition_variable                                                  sync_var;
    extern bool                                                                     sync_flag;

    // Timing stuff
    extern std::map< size_t, std::chrono::time_point<std::chrono::steady_clock> >   tics;
    extern std::mutex                                                               tic_mutex;
}

template <typename Fun, typename ... Args>
void execute_in_parallel(const size_t& n_threads, const Fun& f, const Args& ... args)
{
    ParLib::number_of_threads = n_threads;
    ParLib::thread_ids.resize(n_threads);
    ParLib::thread_pool.resize(n_threads - 1);

    for (size_t i = 0; i < n_threads - 1; i++)
    {
        ParLib::thread_pool[i] = std::thread(f, args ...);
        ParLib::thread_ids[i] = ParLib::thread_pool[i].get_id();
    }

    ParLib::thread_ids.back() = std::this_thread::get_id();
    std::invoke(f, args ...);

    for (auto& t : ParLib::thread_pool)
        t.join();

    ParLib::thread_pool.clear();
    ParLib::thread_ids.clear();
    ParLib::number_of_threads = 0;
}

#endif //KONKURS_PARLIB_HPP
