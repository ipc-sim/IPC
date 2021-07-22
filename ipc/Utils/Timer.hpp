/* Timer.h - Utility class for timing different parts of a program
 Copyright (C) 2012-2014 The University of Reading
 Copying and distribution of this file, with or without modification,
 are permitted in any medium without royalty provided the copyright
 notice and this notice are preserved.  This file is offered as-is,
 without any warranty.
 */

#ifndef Timer_H
#define Timer_H 1

#ifdef _WIN32
#include <windows.h>
#include <time.h>
#else
#include <sys/time.h>
#endif

#include <map>
#include <string>
#include <sstream>
#include <vector>
#include <iostream>

// The Timer class: all functions are inline
class Timer {
public:
    typedef int TimerInt;

    // Constructor can specify a number of unnamed activities
    Timer(TimerInt n_activities = 0)
        : current_activity_(-1), timer_on_(false), print_on_exit_(false)
    {
#ifdef _WIN32
        win_last_time_.QuadPart = 0;
#else
        last_time_.tv_sec = 0;
        last_time_.tv_usec = 0;
#endif
        timings_.reserve(100);
        names_.reserve(100);
        for (TimerInt i = 0; i < n_activities; i++) {
            std::stringstream s;
            s << "Activity " << i;
            timings_.push_back(0.0);
            names_.push_back(s.str());
        }
    }

    // When the timer is destructed (typically at program exit), print
    // out the times spent in each activity
    ~Timer()
    {
        if (print_on_exit_) {
            print();
        }
    }

    // Print out the times spent in each activity
    void print(std::ostream& os = std::cout)
    {
        double sum = 0.0;
        os << timings_.size() << " activities:\n";
        for (TimerInt i = 0; i < timings_.size(); i++) {
            os.width(10);
            os << std::right << timings_[i] << " s: " << names_[i] << "\n";
            sum += timings_[i];
        }
        os.width(10);
        os << std::right << sum << " s: Total\n";
    }

    // Register a new activity with the specified name, returning the
    // tag to be used to specify it in future, as a TimerInt
    TimerInt new_activity(const std::string& name)
    {
        TimerInt tag = timings_.size();
        names_.push_back(name);
        timings_.push_back(0.0);
        return tag;
    }

    // Stop timing current activity
    void stop()
    {
        if (timer_on_) {
            timings_[current_activity_] += split_();
        }
        timer_on_ = false;
    };

    // Start timing specified activity
    void start(TimerInt activity)
    {
        if (timer_on_) {
            timings_[current_activity_] += split_();
        }
        else {
            split_();
        }

        if (activity >= 0 && activity < timings_.size()) {
            current_activity_ = activity;
            timer_on_ = true;
        }
        else {
            // Activity out of range - to keep this inline function fast we
            // don't throw an exception but just don't record the time for
            // this event
            timer_on_ = false;
        }
    };

    // Set the timing for a specific activity back to zero
    void reset(TimerInt activity)
    {
        if (activity >= 0 && activity < timings_.size()) {
            timings_[activity] = 0.0;
        }
    }

    // Return the list of timings in seconds as a constant reference to
    // a vector of doubles
    const std::vector<double>& timings() { return timings_; }

    // Return a single timing
    double timing(TimerInt activity)
    {
        if (activity >= 0 && activity < timings_.size()) {
            return timings_[activity];
        }
        else {
            return 0.0;
        }
    }

    double timing_total(void)
    {
        double total = 0.0;
        for (int i = 0; i < timings_.size(); i++) {
            total += timings_[i];
        }
        return total;
    }

    // Decide whether the contents of the timer class will be printed
    // when it is destructed
    void print_on_exit(bool b = true)
    {
        print_on_exit_ = b;
    }

    Timer& operator+=(const Timer& other)
    {
        assert(timings_.size() == other.timings_.size());
        for (int tI = 0; tI < timings_.size(); ++tI) {
            timings_[tI] += other.timings_[tI];
        }
        return *this;
    }

private:
    // Use Unix system call to get the time accurately
    double split_()
    {
#ifdef _WIN32
        using namespace std;
        QueryPerformanceFrequency(&frequency);
        QueryPerformanceCounter(&win_time_);
        double dsec = (double)(win_time_.QuadPart - win_last_time_.QuadPart)
            / (double)frequency.QuadPart;
        win_last_time_ = win_time_;
        return dsec;
#else
        struct timeval time;
        gettimeofday(&time, NULL);
        double dsec = time.tv_sec - last_time_.tv_sec
            + 0.000001 * (time.tv_usec - last_time_.tv_usec);
        last_time_ = time;
        return dsec;
#endif
    }
    // Data
    std::vector<double> timings_;
    std::vector<std::string> names_;
    TimerInt current_activity_;
#ifdef _WIN32
    LARGE_INTEGER frequency; // ticks per second
    LARGE_INTEGER win_time_, win_last_time_; // ticks
#else
    timeval last_time_;
#endif
    bool timer_on_;
    bool print_on_exit_;
};

#endif
