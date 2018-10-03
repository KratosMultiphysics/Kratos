#ifndef AMGCL_MPI_REPARTITION_RUNTIME_HPP
#define AMGCL_MPI_REPARTITION_RUNTIME_HPP

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
 * \file   amgcl/mpi/partition/runtime.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Runtime wrapper for distributed partitioners.
 */

#include <memory>

#ifdef AMGCL_NO_BOOST
#  error Runtime interface relies on Boost.PropertyTree!
#endif

#include <boost/property_tree/ptree.hpp>

#include <amgcl/util.hpp>
#include <amgcl/mpi/partition/merge.hpp>
#ifdef AMGCL_HAVE_SCOTCH
#  include <amgcl/mpi/partition/ptscotch.hpp>
#endif
#ifdef AMGCL_HAVE_PARMETIS
#  include <amgcl/mpi/partition/parmetis.hpp>
#endif

namespace amgcl {
namespace runtime {
namespace mpi {
namespace partition {

enum type {
    merge
#ifdef AMGCL_HAVE_SCOTCH
  , ptscotch
#endif
#ifdef AMGCL_HAVE_PARMETIS
  , parmetis
#endif
};

std::ostream& operator<<(std::ostream &os, type s)
{
    switch (s) {
        case merge:
            return os << "merge";
#ifdef AMGCL_HAVE_SCOTCH
        case ptscotch:
            return os << "ptscotch";
#endif
#ifdef AMGCL_HAVE_PARMETIS
        case parmetis:
            return os << "parmetis";
#endif
        default:
            return os << "???";
    }
}

std::istream& operator>>(std::istream &in, type &s)
{
    std::string val;
    in >> val;

    if (val == "merge")
        s = merge;
#ifdef AMGCL_HAVE_SCOTCH
    else if (val == "ptscotch")
        s = ptscotch;
#endif
#ifdef AMGCL_HAVE_PARMETIS
    else if (val == "parmetis")
        s = parmetis;
#endif
    else
        throw std::invalid_argument("Invalid partitioner value. Valid choices are: "
                "merge"
#ifdef AMGCL_HAVE_SCOTCH
                ", ptscotch"
#endif
#ifdef AMGCL_HAVE_PARMETIS
                ", parmetis"
#endif
                ".");

    return in;
}

template <class Backend>
struct wrapper {
    typedef amgcl::mpi::distributed_matrix<Backend> matrix;
    typedef boost::property_tree::ptree params;

    type t;
    void *handle;

    wrapper(params prm = params()) : t(prm.get("type",
#if defined(AMGCL_HAVE_SCOTCH)
                ptscotch
#elif defined(AMGCL_HAVE_PARMETIS)
                parmetis
#else
                merge
#endif
                )), handle(0)
    {
        if (!prm.erase("type")) AMGCL_PARAM_MISSING("type");

        switch (t) {
            case merge:
                {
                    typedef amgcl::mpi::partition::merge<Backend> R;
                    handle = static_cast<void*>(new R(prm));
                }
                break;
#ifdef AMGCL_HAVE_SCOTCH
            case ptscotch:
                {
                    typedef amgcl::mpi::partition::ptscotch<Backend> R;
                    handle = static_cast<void*>(new R(prm));
                }
                break;
#endif
#ifdef AMGCL_HAVE_PARMETIS
            case parmetis:
                {
                    typedef amgcl::mpi::partition::parmetis<Backend> R;
                    handle = static_cast<void*>(new R(prm));
                }
                break;
#endif
            default:
                throw std::invalid_argument("Unsupported partition type");
        }
    }

    ~wrapper() {
        switch(t) {
            case merge:
                {
                    typedef amgcl::mpi::partition::merge<Backend> R;
                    delete static_cast<R*>(handle);
                }
                break;
#ifdef AMGCL_HAVE_SCOTCH
            case ptscotch:
                {
                    typedef amgcl::mpi::partition::ptscotch<Backend> R;
                    delete static_cast<R*>(handle);
                }
                break;
#endif
#ifdef AMGCL_HAVE_PARMETIS
            case parmetis:
                {
                    typedef amgcl::mpi::partition::parmetis<Backend> R;
                    delete static_cast<R*>(handle);
                }
                break;
#endif
            default:
                break;
        }
    }

    bool is_needed(const matrix &A) const {
        switch (t) {
            case merge:
                {
                    typedef amgcl::mpi::partition::merge<Backend> R;
                    return static_cast<const R*>(handle)->is_needed(A);
                }
#ifdef AMGCL_HAVE_SCOTCH
            case ptscotch:
                {
                    typedef amgcl::mpi::partition::ptscotch<Backend> R;
                    return static_cast<const R*>(handle)->is_needed(A);
                }
#endif
#ifdef AMGCL_HAVE_PARMETIS
            case parmetis:
                {
                    typedef amgcl::mpi::partition::parmetis<Backend> R;
                    return static_cast<const R*>(handle)->is_needed(A);
                }
#endif
            default:
                throw std::invalid_argument("Unsupported partition type");
        }
    }

    std::shared_ptr<matrix> operator()(const matrix &A, unsigned block_size = 1) const {
        switch (t) {
            case merge:
                {
                    typedef amgcl::mpi::partition::merge<Backend> R;
                    return static_cast<const R*>(handle)->operator()(A, block_size);
                }
#ifdef AMGCL_HAVE_SCOTCH
            case ptscotch:
                {
                    typedef amgcl::mpi::partition::ptscotch<Backend> R;
                    return static_cast<const R*>(handle)->operator()(A, block_size);
                }
#endif
#ifdef AMGCL_HAVE_PARMETIS
            case parmetis:
                {
                    typedef amgcl::mpi::partition::parmetis<Backend> R;
                    return static_cast<const R*>(handle)->operator()(A, block_size);
                }
#endif
            default:
                throw std::invalid_argument("Unsupported partition type");
        }
    }
};

} // namespace partition
} // namespace mpi
} // namespace runtime
} // namespace amgcl

#endif
