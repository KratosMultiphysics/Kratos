#ifndef AMGCL_MPI_COARSENING_RUNTIME_HPP
#define AMGCL_MPI_COARSENING_RUNTIME_HPP

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
 * \file   amgcl/mpi/coarsening/aggregation.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Runtime wrapper for distributed coarsening schemes.
 */

#ifdef AMGCL_NO_BOOST
#  error Runtime interface relies on Boost.PropertyTree!
#endif

#include <boost/property_tree/ptree.hpp>

#include <amgcl/util.hpp>
#include <amgcl/mpi/util.hpp>
#include <amgcl/mpi/distributed_matrix.hpp>
#include <amgcl/mpi/coarsening/aggregation.hpp>
#include <amgcl/mpi/coarsening/smoothed_aggregation.hpp>

namespace amgcl {
namespace runtime {
namespace mpi {
namespace coarsening {

enum type {
    aggregation,
    smoothed_aggregation
};

std::ostream& operator<<(std::ostream &os, type s)
{
    switch (s) {
        case aggregation:
            return os << "aggregation";
        case smoothed_aggregation:
            return os << "smoothed_aggregation";
        default:
            return os << "???";
    }
}

std::istream& operator>>(std::istream &in, type &s)
{
    std::string val;
    in >> val;

    if (val == "aggregation")
        s = aggregation;
    else if (val == "smoothed_aggregation")
        s = aggregation;
    else
        throw std::invalid_argument("Invalid coarsening value. Valid choices are: "
                "aggregation, smoothed_aggregation.");

    return in;
}

template <class Backend>
struct wrapper {
    typedef amgcl::mpi::distributed_matrix<Backend> matrix;
    typedef boost::property_tree::ptree params;

    type c;
    void *handle;

    wrapper(params prm = params())
        : c(prm.get("type", smoothed_aggregation)), handle(0)
    {
        if (!prm.erase("type")) AMGCL_PARAM_MISSING("type");

        switch (c) {
            case aggregation:
                {
                    typedef amgcl::mpi::coarsening::aggregation<Backend> C;
                    handle = static_cast<void*>(new C(prm));
                }
                break;
            case smoothed_aggregation:
                {
                    typedef amgcl::mpi::coarsening::smoothed_aggregation<Backend> C;
                    handle = static_cast<void*>(new C(prm));
                }
                break;
            default:
                throw std::invalid_argument("Unsupported coarsening type");
        }
    }

    ~wrapper() {
        switch(c) {
            case aggregation:
                {
                    typedef amgcl::mpi::coarsening::aggregation<Backend> C;
                    delete static_cast<C*>(handle);
                }
                break;
            case smoothed_aggregation:
                {
                    typedef amgcl::mpi::coarsening::smoothed_aggregation<Backend> C;
                    delete static_cast<C*>(handle);
                }
                break;
            default:
                break;
        }
    }

    std::tuple< std::shared_ptr<matrix>, std::shared_ptr<matrix> >
    transfer_operators(const matrix &A) {
        switch (c) {
            case aggregation:
                {
                    typedef amgcl::mpi::coarsening::aggregation<Backend> C;
                    return static_cast<C*>(handle)->transfer_operators(A);
                }
            case smoothed_aggregation:
                {
                    typedef amgcl::mpi::coarsening::smoothed_aggregation<Backend> C;
                    return static_cast<C*>(handle)->transfer_operators(A);
                }
            default:
                throw std::invalid_argument("Unsupported partition type");
        }
    }

    std::shared_ptr<matrix>
    coarse_operator(const matrix &A, const matrix &P, const matrix &R) const {
        switch (c) {
            case aggregation:
                {
                    typedef amgcl::mpi::coarsening::aggregation<Backend> C;
                    return static_cast<C*>(handle)->coarse_operator(A, P, R);
                }
            case smoothed_aggregation:
                {
                    typedef amgcl::mpi::coarsening::smoothed_aggregation<Backend> C;
                    return static_cast<C*>(handle)->coarse_operator(A, P, R);
                }
            default:
                throw std::invalid_argument("Unsupported partition type");
        }
    }
};

template <class Backend>
unsigned block_size(const wrapper<Backend> &w) {
    switch (w.c) {
        case aggregation:
            {
                typedef amgcl::mpi::coarsening::aggregation<Backend> C;
                return block_size(*static_cast<const C*>(w.handle));
            }
        case smoothed_aggregation:
            {
                typedef amgcl::mpi::coarsening::smoothed_aggregation<Backend> C;
                return block_size(*static_cast<const C*>(w.handle));
            }
        default:
            throw std::invalid_argument("Unsupported coarsening type");
    }
}

} // namespace coarsening
} // namespace mpi
} // namespace runtime
} // namespace amgcl

#endif
