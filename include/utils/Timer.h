#include <chrono>

class Timer {
public:
    /**
    * @brief Starts the timer.
    */
    void start() { begin = std::chrono::steady_clock::now(); }

    /**
    * @brief Returns the elapsed time since the start of the timer.
    */
    double count() const {
        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> et{end - begin};
        return et.count();
    }

private:
    std::chrono::steady_clock::time_point begin;
};