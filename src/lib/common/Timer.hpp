#pragma once

#include <boost/chrono/chrono_io.hpp>

#include <cstddef>
#include <ostream>
#include <string>
#include <sstream>



template<typename Clock>
class Timer {
public:
    Timer() : _start(Clock::now()) {}

    template<typename T>
    T elapsed() const {
        return boost::chrono::duration_cast<T>(Clock::now() - _start);
    }

private:
    typename Clock::time_point _start;
};
