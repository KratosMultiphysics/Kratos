#ifndef AMGCL_PROFILER_H
#define AMGCL_PROFILER_H

/*
The MIT License

Copyright (c) 2012-2017 Denis Demidov <dennis.demidov@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

/**
 * \file   amgcl/profiler.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Profiler class.
 */

#include <iostream>
#include <iomanip>
#include <map>
#include <string>
#include <vector>

#include <boost/type_traits.hpp>
#include <boost/io/ios_state.hpp>
#include <amgcl/perf_counter/clock.hpp>


namespace amgcl {

/// Profiler class.
/**
 * \param Counter     Performance counter to use for profiling.
 * \param SHIFT_WIDTH Indentation for output of profiling results.
 *
 * Provides simple to use, hierarchical profile with nicely formatted output.
 */
template <class Counter = amgcl::perf_counter::clock, unsigned SHIFT_WIDTH = 2>
class profiler {
    public:
        typedef typename Counter::value_type value_type;

        /// Initialization.
        /**
         * \param name Profile title to use with output.
         */
        profiler(const std::string &name = "Profile") : name(name) {
            stack.reserve(128);
            stack.push_back(&root);
            root.begin = counter.current();
        }

        /// Starts measurement.
        /**
         * \param name interval name.
         */
        void tic(const std::string &name) {
            stack.back()->children[name].begin = counter.current();
            stack.push_back(&stack.back()->children[name]);
        }

        /// Stops measurement.
        /**
         * Returns delta in the measured value since the corresponding tic().
         */
        value_type toc(const std::string& /*name*/) {
            profile_unit *top = stack.back();
            stack.pop_back();

            value_type current = counter.current();
            value_type delta   = current - top->begin;

            top->length += delta;
            root.length = current - root.begin;

            return delta;
        }

        void reset() {
            stack.clear();
            root.length = 0;
            root.children.clear();

            stack.push_back(&root);
            root.begin = counter.current();
        }

    private:
        struct profile_unit {
            profile_unit() : length(0) {}

            value_type children_time() const {
                value_type s = value_type();
                for(typename std::map<std::string, profile_unit>::const_iterator c = children.begin(); c != children.end(); c++)
                    s += c->second.length;
                return s;
            }

            size_t total_width(const std::string &name, int level) const {
                size_t w = name.size() + level;
                for(typename std::map<std::string, profile_unit>::const_iterator c = children.begin(); c != children.end(); c++)
                    w = std::max(w, c->second.total_width(c->first, level + SHIFT_WIDTH));
                return w;
            }

            void print(std::ostream &out, const std::string &name,
                    int level, value_type total, size_t width) const
            {
                using namespace std;

                out << "[" << setw(level) << "";
                print_line(out, name, length, 100 * length / total, width - level);

                if (children.size()) {
                    value_type val = length - children_time();
                    double perc = 100.0 * val / total;

                    if (perc > 1e-1) {
                        out << "[" << setw(level + 1) << "";
                        print_line(out, "self", val, perc, width - level - 1);
                    }
                }

                for(typename std::map<std::string, profile_unit>::const_iterator c = children.begin(); c != children.end(); c++)
                    c->second.print(out, c->first, level + SHIFT_WIDTH, total, width);
            }

            void print_line(std::ostream &out, const std::string &name,
                    value_type time, double perc, size_t width) const
            {
                using namespace std;

                out << name << ":"
                    << setw(width - name.size()) << ""
                    << setw(10)
                    << fixed << setprecision(3) << time << " " << Counter::units()
                    << "] (" << fixed << setprecision(2) << setw(6) << perc << "%)"
                    << endl;
            }

            value_type begin;
            value_type length;

            std::map<std::string, profile_unit> children;
        };

        Counter counter;
        std::string name;
        profile_unit root;
        std::vector<profile_unit*> stack;

        void print(std::ostream &out) {
            if (stack.back() != &root)
                out << "Warning! Profile is incomplete." << std::endl;

            boost::io::ios_all_saver stream_state(out);
            root.print(out, name, 0, root.length, root.total_width(name, 0));
        }

        /// Sends formatted profiling data to an output stream.
        /**
         * \param out  Output stream.
         * \param prof Profiler.
         */
        friend std::ostream& operator<<(std::ostream &out, profiler &prof) {
            out << std::endl;
            prof.print(out);
            return out << std::endl;
        }
};

/// Scoped ticker.
/** Calls prof.tic(name) on construction, and prof.toc(name) on destruction. */
template <class Profiler>
struct scoped_tic {
    Profiler &prof;
    std::string name;

    scoped_tic(Profiler &prof, const std::string &name)
        : prof(prof), name(name)
    {
        prof.tic(name);
    }

    ~scoped_tic() {
        prof.toc(name);
    }
};

} // namespace amgcl

#endif
