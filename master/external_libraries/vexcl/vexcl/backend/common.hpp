#ifndef VEXCL_BACKEND_COMMON_HPP
#define VEXCL_BACKEND_COMMON_HPP

/*
The MIT License

Copyright (c) 2012-2018 Denis Demidov <dennis.demidov@gmail.com>

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
 * \file   vexcl/backend/common.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Common backend utilities.
 */

#include <vector>
#include <map>

#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdlib>

#include <boost/version.hpp>
#if BOOST_VERSION >= 106600
#  include <boost/uuid/detail/sha1.hpp>
#else
#  include <boost/uuid/sha1.hpp>
#endif
#include <boost/optional.hpp>
#include <boost/filesystem.hpp>
#include <boost/thread.hpp>

namespace vex {

enum device_options_kind {
    compile_options,    ///< Options sent to the compute kernel compiler.
    program_header      ///< Program header prepended to all compute kernel source.
};

/// Global program options holder
template <device_options_kind kind>
struct device_options {
    static std::string get(const backend::command_queue &q) {
        auto dev = backend::get_device_id(q);

        std::ostringstream s;

        boost::lock_guard<boost::mutex> lock(options_mx);
        for(auto o = options[dev].begin(); o != options[dev].end(); ++o)
            s << (kind == compile_options ? " " : "\n") << (*o);
        return s.str();
    }

    static void push(const backend::command_queue &q, const std::string &str) {
        auto dev = backend::get_device_id(q);

        boost::lock_guard<boost::mutex> lock(options_mx);
        options[dev].push_back(str);
    }

    static void pop(const backend::command_queue &q) {
        auto dev = backend::get_device_id(q);

        boost::lock_guard<boost::mutex> lock(options_mx);
        if (!options[dev].empty()) options[dev].pop_back();
    }

    private:
        static std::map<backend::device_id, std::vector<std::string> > options;
        static boost::mutex options_mx;
};

template <device_options_kind kind>
std::map<backend::device_id, std::vector<std::string> > device_options<kind>::options;

template <device_options_kind kind>
boost::mutex device_options<kind>::options_mx;

inline std::string get_compile_options(const backend::command_queue &q) {
    return device_options<compile_options>::get(q);
}

inline std::string get_program_header(const backend::command_queue &q) {
    return device_options<program_header>::get(q);
}

/// Set global compute kernel compilation options for a given device.
/**
 * This replaces any previously set options. To roll back, call
 * pop_compile_options().
 */
inline void push_compile_options(const backend::command_queue &q, const std::string &str) {
    device_options<compile_options>::push(q, str);
}

/// Rolls back changes to compile options.
inline void pop_compile_options(const backend::command_queue &q) {
    device_options<compile_options>::pop(q);
}

/// Set global compute kernel header for a given device.
/**
 * This replaces any previously set header. To roll back, call
 * pop_program_header().
 */
inline void push_program_header(const backend::command_queue &q, const std::string &str) {
    device_options<program_header>::push(q, str);
}

/// Rolls back changes to compile options.
inline void pop_program_header(const backend::command_queue &q) {
    device_options<program_header>::pop(q);
}

/// Set global compute kernel compilation options for each device in queue list.
inline void push_compile_options(const std::vector<backend::command_queue> &queue, const std::string &str) {
    for(auto q = queue.begin(); q != queue.end(); ++q)
        device_options<compile_options>::push(*q, str);
}

/// Rolls back changes to compile options for each device in queue list.
inline void pop_compile_options(const std::vector<backend::command_queue> &queue) {
    for(auto q = queue.begin(); q != queue.end(); ++q)
        device_options<compile_options>::pop(*q);
}

/// Pushes compile options on construction, pops them on destruction
struct scoped_compile_options {
    std::vector<backend::command_queue> q;

    scoped_compile_options(
            const std::vector<backend::command_queue> &q,
            const std::string &str
            ) : q(q)
    {
        push_compile_options(q, str);
    }

    scoped_compile_options(
            const backend::command_queue &q,
            const std::string &str
            ) : q(1, q)
    {
        push_compile_options(this->q, str);
    }

    ~scoped_compile_options() {
        pop_compile_options(q);
    }
};

/// Set global OpenCL program header for each device in queue list.
inline void push_program_header(const std::vector<backend::command_queue> &queue, const std::string &str) {
    for(auto q = queue.begin(); q != queue.end(); ++q)
        device_options<program_header>::push(*q, str);
}

/// Rolls back changes to compile options for each device in queue list.
inline void pop_program_header(const std::vector<backend::command_queue> &queue) {
    for(auto q = queue.begin(); q != queue.end(); ++q)
        device_options<program_header>::pop(*q);
}

/// Pushes program header on construction, pops it on destruction
struct scoped_program_header {
    std::vector<backend::command_queue> q;

    scoped_program_header(
            const std::vector<backend::command_queue> &q,
            const std::string &str
            ) : q(q)
    {
        push_program_header(q, str);
    }

    scoped_program_header(
            const backend::command_queue &q,
            const std::string &str
            ) : q(1, q)
    {
        push_program_header(this->q, str);
    }

    ~scoped_program_header() {
        pop_program_header(q);
    }
};

/// Path delimiter symbol.
inline const std::string& path_delim() {
    static const std::string delim = boost::filesystem::path("/").make_preferred().string();
    return delim;
}

/// Path to appdata folder.
inline const std::string& appdata_path() {
#ifdef _WIN32
    static const std::string appdata = getenv("APPDATA") + path_delim() + "vexcl";
#else
    static const std::string appdata = getenv("HOME") + path_delim() + ".vexcl";
#endif
    return appdata;
}

/// Path to cached binaries.
inline std::string program_binaries_path(const std::string &hash, bool create = false)
{
    std::string dir = appdata_path()    + path_delim()
                    + hash.substr(0, 2) + path_delim()
                    + hash.substr(2);
    if (create) boost::filesystem::create_directories(dir);
    return dir + path_delim();
}

/// SHA1 hasher.
class sha1_hasher {
    public:
        sha1_hasher(const std::string &s = "") {
            if (!s.empty()) this->process(s);
        }

        sha1_hasher& process(const std::string &s) {
            h.process_bytes(s.c_str(), s.size());
            return *this;
        }

        operator std::string() {
            boost::uuids::detail::sha1::digest_type digest;
            h.get_digest(digest);
            return serialize_digest(digest);
        }
    private:
        // Before boost 1.86.0, boost::uuids::detail::sha1 used to be unsigned int[5].
        std::string serialize_digest(unsigned int digest[5]) const {
            std::ostringstream buf;
            for(int i = 0; i < 5; ++i)
                buf << std::hex << std::setfill('0') << std::setw(8) << digest[i];
            return buf.str();
        }

        // Since boost 1.86.0, boost::uuids::detail::sha1 has been unsigned char[20].
        std::string serialize_digest(unsigned char digest[20]) const {
            // Convert digest to a string representation of 5 integers in hexadecimal format
            // (digest_type used to be unsigned int[5]). The reason for boost's change is to
            // fix a bug related to endianness, but to conform to old VexCL behavior, bytes
            // are reordered if necessary (whether conformance is required is up for debate).
            std::ostringstream buf;
            #if (defined(__BYTE_ORDER__) && (__BYTE_ORDER__ == __ORDER_BIG_ENDIAN__)) || defined(__BIG_ENDIAN__) || defined(__BIG_ENDIAN) || defined(_BIG_ENDIAN)
                const int int_size = sizeof(int);
                for (int i=int_size; i<=20; i+=int_size) {
                    for (int j=i-1; j>=i-int_size; --j) {
                        buf << std::hex << std::setfill('0') << std::setw(2) << (int)digest[j];
                    }
                }
            #else
                for (unsigned i=0; i<20; ++i) {
                    buf << std::hex << std::setfill('0') << std::setw(2) << (int)digest[i];
                }
            #endif

            return buf.str();
        }

        boost::uuids::detail::sha1 h;
};

} // namespace vex


#endif
