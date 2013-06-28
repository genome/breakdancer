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

template<typename Clock, typename Units>
class ScopedTimer {
public:
    ScopedTimer(
            std::ostream& stream,
            std::string const& message,
            char indent_char = ' ',
            unsigned short indent_width = 4u,
            unsigned short indent_level = 0u
            )
        : _stream(stream)
        , _message(message)
        , _indent_char(indent_char)
        , _indent_width(indent_width)
        , _indent_level(indent_level)
        , _indent(indent_width * indent_level, indent_char)
    {
        _stream << _indent << "BEGIN: " << message << "\n";
    }

    ~ScopedTimer() {
        printExtra();
        _stream << _indent<< "END: " << _message << ", duration = "
            << _timer.template elapsed<Units>()
            << ".\n";
    }

    ScopedTimer subtimer(std::string const& message) {
        return ScopedTimer(_stream, message, _indent_char, _indent_width,
            _indent_level + 1);
    }

    template<typename T>
    ScopedTimer& operator<<(T const& x) {
        _extra << x;
        return *this;
    }

private:
    void printExtra() {
        std::string line;
        while (getline(_extra, line))
            _stream << _indent << " >" << line << "\n";
    }

private:
    Timer<Clock> _timer;
    std::ostream& _stream;
    std::string _message;
    std::stringstream _extra;
    char _indent_char;
    unsigned short _indent_width;
    unsigned short _indent_level;
    std::string _indent;
};
