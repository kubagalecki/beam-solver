#ifndef KONKURS_THREADPOOL_HPP
#define KONKURS_THREADPOOL_HPP

#include <cstddef>

#include <atomic>
#include <condition_variable>
#include <mutex>
#include <thread>

#include <queue>
#include <vector>

#include <algorithm>

#include <functional>

class ThreadPool
{
public:
    ThreadPool(const size_t nt = 1);
    ThreadPool(const ThreadPool&) = delete;
    ThreadPool& operator=(const ThreadPool&) = delete;
    ThreadPool(ThreadPool&&)                 = default;
    ThreadPool& operator=(ThreadPool&&) = default;
    ~ThreadPool();

    template < typename F, typename... Args >
    void push_job(F&& fun, Args&&... args);

    void complete();

    size_t size() const;

private:
    std::condition_variable   wait_var;
    std::mutex                pool_mutex;
    std::atomic_bool          terminate_flag;
    std::atomic_uint_fast16_t n_busy;
    std::condition_variable   complete_var;

    std::queue< std::function< void() > > work_queue;
    std::vector< std::thread >            thread_pool;

    auto get_exec_fun();
};

template < typename F, typename... Args >
void ThreadPool::push_job(F&& fun, Args&&... args)
{
    {
        std::lock_guard< std::mutex > lock(pool_mutex);
        work_queue.emplace([f = std::forward< F >(fun), ... a = std::forward< Args >(args)]() {
            f(std::move(a)...);
        });
    }
    wait_var.notify_one();
}

#endif // KONKURS_THREADPOOL_HPP
