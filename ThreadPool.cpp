#include "ThreadPool.hpp"

auto ThreadPool::get_exec_fun()
{
    return [this](){
        while (true)
        {
            std::function<void()> work;

            {
                std::unique_lock<std::mutex> lock(pool_mutex);
                wait_var.wait(lock, [&]() {return terminate_flag || !work_queue.empty();});

                if (!work_queue.empty())
                {
                    work = std::move(work_queue.front());
                    work_queue.pop();
                    ++n_busy;
                }

                if (terminate_flag)
                    break;
            }

            if (work)
            {
                work();
                --n_busy;

                if (n_busy == 0)
                {
                    std::unique_lock<std::mutex> lock(pool_mutex);
                    complete_var.notify_all();
                }
            }
        }

        wait_var.notify_one();
    };
}

ThreadPool::ThreadPool(size_t nt) : thread_pool(nt), terminate_flag(false), n_busy(0)
{
    std::generate(thread_pool.begin(), thread_pool.end(),
            [exec_fun = get_exec_fun()](){return std::thread{exec_fun};});
}

ThreadPool::~ThreadPool()
{
    terminate_flag = true;

    wait_var.notify_one();

    for (auto& t : thread_pool)
        t.join();
}

void ThreadPool::complete()
{
    if (thread_pool.empty())
        return;

    std::unique_lock<std::mutex> lock(pool_mutex);
    complete_var.wait(lock, [this](){return n_busy == 0 && work_queue.empty();});
}

size_t ThreadPool::size() const
{
    return thread_pool.size();
}
