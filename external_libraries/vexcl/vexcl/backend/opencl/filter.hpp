#ifndef VEXCL_BACKEND_OPENCL_FILTER_HPP
#define VEXCL_BACKEND_OPENCL_FILTER_HPP

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
 * \file   vexcl/backend/opencl/filter.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Device filters for OpenCL backend.
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <functional>
#include <cstdlib>

#include <boost/interprocess/sync/file_lock.hpp>
#include <boost/date_time/posix_time/posix_time_types.hpp>
#include <boost/filesystem.hpp>
#include <boost/config.hpp>
#include <boost/version.hpp>

#include <vexcl/backend/opencl/defines.hpp>
#ifdef VEXCL_HAVE_OPENCL_HPP
#  include <CL/opencl.hpp>
#else
#  include <CL/cl2.hpp>
#endif


namespace vex {

/// Device filters.
namespace Filter {

    /// Selects devices whose vendor name match given value.
    struct Vendor {
        explicit Vendor(std::string name) : vendor(std::move(name)) {}

        bool operator()(const cl::Device &d) const {
            return d.getInfo<CL_DEVICE_VENDOR>().find(vendor) != std::string::npos;
        }

        private:
            std::string vendor;
    };

    /// Selects devices whose platform name match given value.
    struct Platform {
        explicit Platform(std::string name) : platform(std::move(name)) {}

        bool operator()(const cl::Device &d) const {
            return cl::Platform(d.getInfo<CL_DEVICE_PLATFORM>()).getInfo<CL_PLATFORM_NAME>().find(platform) != std::string::npos;
        }

        private:
            std::string platform;
    };

    /// Selects devices whose names match given value.
    struct Name {
        explicit Name(std::string name) : devname(std::move(name)) {}

        bool operator()(const cl::Device &d) const {
            return d.getInfo<CL_DEVICE_NAME>().find(devname) != std::string::npos;
        }

        private:
            std::string devname;
    };

    /// Selects devices by type.
    struct Type {
        explicit Type(cl_device_type t)    : type(t)              {}
        explicit Type(const std::string t) : type(device_type(t)) {}

        bool operator()(const cl::Device &d) const {
            return 0 != (d.getInfo<CL_DEVICE_TYPE>() & type);
        }

        private:
            cl_device_type type;

            static cl_device_type device_type(const std::string &t) {
                if (t.find("CPU") != std::string::npos)
                    return CL_DEVICE_TYPE_CPU;

                if (t.find("GPU") != std::string::npos)
                    return CL_DEVICE_TYPE_GPU;

                if (t.find("ACCELERATOR") != std::string::npos)
                    return CL_DEVICE_TYPE_ACCELERATOR;

                return CL_DEVICE_TYPE_ALL;
            }
    };

    /// Selects CPU devices
    const Type CPU        (CL_DEVICE_TYPE_CPU);
    /// Selects GPU devices
    const Type GPU        (CL_DEVICE_TYPE_GPU);
    /// Selects accelerator devices
    const Type Accelerator(CL_DEVICE_TYPE_ACCELERATOR);

    struct DoublePrecisionFilter {
        bool operator()(const cl::Device &d) const {
            std::string ext = d.getInfo<CL_DEVICE_EXTENSIONS>();
            return (
                    ext.find("cl_khr_fp64") != std::string::npos ||
                    ext.find("cl_amd_fp64") != std::string::npos
                   );
        }
    };

    /// Selects devices supporting double precision.
    const DoublePrecisionFilter DoublePrecision = {};

    /// Selects devices providing given OpenCL extensions.
    struct Extension {
        std::string extension;

        Extension(std::string ext) : extension(std::move(ext)) {}

        bool operator()(const cl::Device& d) const {
            std::string ext = d.getInfo<CL_DEVICE_EXTENSIONS>();
            return ext.find(extension) != std::string::npos;
        }
    };

    /// Selects devices supporting OpenGL Sharing Extension
    const Extension GLSharing(
    #if defined(__APPLE__)
        "cl_apple_gl_sharing"
    #else
        "cl_khr_gl_sharing"
    #endif
    );

    /// Selects device by OpenCL version.
    struct CLVersion {
        std::tuple<int, int> min_version;

        CLVersion(int major, int minor) : min_version(major, minor) {}

        bool operator()(const cl::Device &d) const {
            // version contains OpenCL n.m
            //                  0123456789
            std::string version = d.getInfo<CL_DEVICE_VERSION>();

            return min_version <= std::make_tuple(version[7]-'0', version[9]-'0');
        }
    };

    /// List of device filters based on environment variables.
    inline std::vector< std::function<bool(const cl::Device&)> >
    backend_env_filters()
    {
        std::vector< std::function<bool(const cl::Device&)> > filter;

        if (const char *platform = getenv("OCL_PLATFORM"))
            filter.push_back(Platform(platform));

        if (const char *vendor = getenv("OCL_VENDOR"))
            filter.push_back(Vendor(vendor));

        if (const char *name = getenv("OCL_DEVICE"))
            filter.push_back(Name(name));

        if (const char *devtype = getenv("OCL_TYPE"))
            filter.push_back(Type(devtype));

        if (const char *extension = getenv("OCL_EXTENSION"))
            filter.push_back(Extension(extension));

        return filter;
    }

    /// \internal Exclusive access to selected devices.
    class ExclusiveFilter {
        private:
            std::function<bool(const cl::Device&)> filter;

            static std::map<cl_device_id, std::string> get_uids() {
                std::map<cl_device_id, std::string> uids;

                std::vector<cl::Platform> platform;
                cl::Platform::get(&platform);

                const char *lock_dir = getenv("VEXCL_LOCK_DIR");

                for(size_t p_id = 0; p_id < platform.size(); p_id++) {
                    std::vector<cl::Device> device;

                    platform[p_id].getDevices(CL_DEVICE_TYPE_ALL, &device);

                    for(size_t d_id = 0; d_id < device.size(); d_id++) {
                        std::ostringstream id;
#ifdef _WIN32
                        id << (lock_dir ? lock_dir : getenv("TEMP")) << "\\";
#else
                        id << (lock_dir ? lock_dir : "/tmp") << "/";
#endif
                        id << "vexcl_device_" << p_id << "_" << d_id << ".lock";

                        uids[device[d_id]()] = id.str();
                    }
                }

                return uids;
            }

            struct locker {
                locker(std::string fname) : file(fname)
                {
                    if (!file.is_open() || file.fail()) {
                        std::cerr
                            << "WARNING: failed to open file \"" << fname << "\"\n"
                            << "  Check that target directory is exists and is writable.\n"
                            << "  Exclusive mode is off.\n"
                            << std::endl;
                    } else {
                        flock.reset(new boost::interprocess::file_lock(fname.c_str()));
#if BOOST_VERSION >= 105000
                        // In case we created the lock file,
                        // lets make it writable by others:
                        try {
                            boost::filesystem::permissions(fname, boost::filesystem::all_all);
                        } catch (const boost::filesystem::filesystem_error&) {
                            (void)0;
                        }
#endif
                    }
                }

                bool try_lock() {
                    if (flock) {
                        // Try and lock the file related to compute device.
                        // If the file is locked already, it could mean two
                        // things:
                        // 1. Somebody locked the file, and uses the device.
                        // 2. Somebody locked the file, and is in process of
                        //    checking the device. If device is not good (for
                        //    them) they will release the lock in a few
                        //    moments.
                        // To process case 2 correctly, we use timed_lock().

                        return flock->timed_lock(
                                boost::posix_time::microsec_clock::universal_time() +
                                boost::posix_time::milliseconds(100)
                                );
                    } else return true;
                }

                std::ofstream file;
                std::unique_ptr<boost::interprocess::file_lock> flock;
            };
        public:
            template <class Filter>
            ExclusiveFilter(Filter&& filter)
                : filter(std::forward<Filter>(filter)) {}

            bool operator()(const cl::Device &d) const {
                static std::map<cl_device_id, std::string> dev_uids = get_uids();
                static /*thread_local*/ std::vector<std::unique_ptr<locker>> locks;

                std::unique_ptr<locker> lck(new locker(dev_uids[d()]));

                if (lck->try_lock() && filter(d)) {
                    locks.push_back(std::move(lck));
                    return true;
                }

                return false;
            }

    };

    /// Allows exclusive access to compute devices across several processes.
    /**
     * Returns devices that pass through provided device filter and are not
     * locked.
     *
     * \param filter Compute device filter
     *
     * \note Depends on boost::interprocess library.
     *
     * lock files are created in directory specified in VEXCL_LOCK_DIR
     * environment variable. If the variable does not exist, /tmp is
     * used on Linux and %TMPDIR% on Windows. The lock directory should exist
     * and be writable by the running user.
     */
    template <class Filter>
    ExclusiveFilter Exclusive(Filter&& filter) {
        return ExclusiveFilter(std::forward<Filter>(filter));
    }

} // namespace Filter
} // namespace vex

#endif
