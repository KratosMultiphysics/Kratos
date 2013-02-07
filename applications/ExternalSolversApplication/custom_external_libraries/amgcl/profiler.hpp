#ifndef AMGCL_PROFILER_H
#define AMGCL_PROFILER_H

/*
The MIT License

Copyright (c) 2012-2013 Denis Demidov <ddemidov@ksu.ru>

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
 * \file   profiler.hpp
 * \author Denis Demidov <ddemidov@ksu.ru>
 * \brief  Profiler class.
 */

#include <iostream>
#include <iomanip>
#include <map>
#include <string>
#include <stack>

#include <boost/chrono.hpp>
#include <boost/typeof/typeof.hpp>


namespace amgcl {

/// Simple profiler class.
/**
 * \param clock       Clock to use for profiling.
 * \param SHIFT_WIDTH Indentation for output of profiling results.
 *
 * Provides simple to use, possibly nested timers with nicely formatted output.
 * The implementation was inspired by a HOOMD-Blue code
 * (http://codeblue.umich.edu/hoomd-blue)
 */
template <
    class clock = boost::chrono::high_resolution_clock,
    unsigned SHIFT_WIDTH = 2
    >
class profiler {
    public:
        /// Initialization.
        /**
         * \param name Profile title to use with output.
         */
        profiler(const std::string &name = "Profile") : name(name) {
            stack.push(&root);
            root.start_time = clock::now();
        }

        /// Starts named timer.
        /**
         * \param name Timer name.
         */
        void tic(const std::string &name) {
            stack.top()->children[name].start_time = clock::now();
            stack.push(&stack.top()->children[name]);
        }

        /// Stops named timer.
        double toc(const std::string& /*name*/) {
            BOOST_AUTO(top, stack.top());
            stack.pop();

            double delta = boost::chrono::duration<double>(
                    clock::now() - top->start_time).count();

            top->length += delta;

            return delta;
        }

    private:
        struct profile_unit {
            profile_unit() : length(0) {}

            double children_time() const {
                double s = 0;
                for(BOOST_AUTO(c, children.begin()); c != children.end(); c++)
                    s += c->second.length;
                return s;
            }

            size_t total_width(const std::string &name, int level) const {
                size_t w = name.size() + level;
                for(BOOST_AUTO(c, children.begin()); c != children.end(); c++)
                    w = std::max(w, c->second.total_width(c->first, level + SHIFT_WIDTH));
                return w;
            }

            void print(std::ostream &out, const std::string &name,
                    int level, double total, size_t width) const
            {
                using namespace std;

                out << "[" << setw(level) << "";
                print_line(out, name, length, 100 * length / total, width - level);

                if (children.size()) {
                    double sec = length - children_time();
                    double perc = 100 * sec / total;

                    if (perc > 1e-1) {
                        out << "[" << setw(level + 1) << "";
                        print_line(out, "self", sec, perc, width - level - 1);
                    }
                }

                for(BOOST_AUTO(c, children.begin()); c != children.end(); c++)
                    c->second.print(out, c->first, level + SHIFT_WIDTH, total, width);
            }

            void print_line(std::ostream &out, const std::string &name,
                    double time, double perc, size_t width) const
            {
                using namespace std;

                out << name << ":"
                    << setw(width - name.size()) << ""
                    << setw(10)
                    << fixed << setprecision(3) << time << " sec."
                    << "] (" << fixed << setprecision(2) << setw(6) << perc << "%)"
                    << endl;
            }

            boost::chrono::time_point<clock> start_time;

            double length;

            std::map<std::string, profile_unit> children;
        };

        std::string name;
        profile_unit root;
        std::stack<profile_unit*> stack;

        void print(std::ostream &out) {
            if (stack.top() != &root)
                out << "Warning! Profile is incomplete." << std::endl;

            root.length += boost::chrono::duration<double>(
                    clock::now() - root.start_time).count();

            BOOST_AUTO(ff, out.flags());
            BOOST_AUTO(pp, out.precision());

            root.print(out, name, 0, root.length, root.total_width(name, 0));

            out.flags(ff);
            out.precision(pp);
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

} // namespace amgcl

#endif
