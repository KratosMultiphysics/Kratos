#ifndef AMGCL_PROFILER_H
#define AMGCL_PROFILER_H

/*
The MIT License

Copyright (c) 2012-2015 Denis Demidov <dennis.demidov@gmail.com>

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
#include <amgcl/clock.hpp>


namespace amgcl {

/// Returns difference (in seconds) between two time points.
/** The time points should come either from boost::chrono or std::chrono */
template <class TP>
inline
typename boost::enable_if<boost::is_class<typename TP::duration>, double>::type
seconds(TP tic, TP toc) {
    return static_cast<double>(TP::duration::period::num)
        * typename TP::duration(toc - tic).count()
        / TP::duration::period::den;
}

/// Profiler class.
/**
 * \param clock       Clock to use for profiling.
 * \param SHIFT_WIDTH Indentation for output of profiling results.
 *
 * Provides simple to use, hierarchical timers with nicely formatted output.
 */
template <class clock = amgcl::clock, unsigned SHIFT_WIDTH = 2>
class profiler {
    public:
        /// Initialization.
        /**
         * \param name Profile title to use with output.
         */
        profiler(const std::string &name = "Profile") : name(name) {
            stack.reserve(128);
            stack.push_back(&root);
            root.start_time = clock::now();
        }

        /// Starts named timer.
        /**
         * \param name Timer name.
         */
        void tic(const std::string &name) {
            stack.back()->children[name].start_time = clock::now();
            stack.push_back(&stack.back()->children[name]);
        }

        /// Stops named timer.
        double toc(const std::string& /*name*/) {
            profile_unit *top = stack.back();
            stack.pop_back();

            double delta = seconds(top->start_time, clock::now());
            top->length += delta;
            return delta;
        }

        void reset() {
            stack.clear();
            root.length = 0;
            root.children.clear();

            stack.push_back(&root);
            root.start_time = clock::now();
        }

    private:
        struct profile_unit {
            profile_unit() : length(0) {}

            double children_time() const {
                double s = 0;
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

                for(typename std::map<std::string, profile_unit>::const_iterator c = children.begin(); c != children.end(); c++)
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

            typename clock::time_point start_time;

            double length;

            std::map<std::string, profile_unit> children;
        };

        std::string name;
        profile_unit root;
        std::vector<profile_unit*> stack;

        void print(std::ostream &out) {
            if (stack.back() != &root)
                out << "Warning! Profile is incomplete." << std::endl;

            root.length += seconds(root.start_time, clock::now());

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

} // namespace amgcl

#endif
