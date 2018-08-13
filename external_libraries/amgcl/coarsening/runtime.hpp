#ifndef AMGCL_COARSENING_RUNTIME_HPP
#define AMGCL_COARSENING_RUNTIME_HPP

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
 * \file   amgcl/coarsening/runtime.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Runtime-configurable coarsening.
 */

#include <iostream>
#include <stdexcept>
#include <type_traits>

#ifdef AMGCL_NO_BOOST
#  error Runtime interface relies on Boost.PropertyTree!
#endif

#include <boost/property_tree/ptree.hpp>

#include <amgcl/util.hpp>
#include <amgcl/backend/interface.hpp>
#include <amgcl/coarsening/ruge_stuben.hpp>
#include <amgcl/coarsening/aggregation.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/coarsening/smoothed_aggr_emin.hpp>

namespace amgcl {
namespace runtime {

/// Coarsening kinds.
namespace coarsening {

enum type {
    ruge_stuben,            ///< Ruge-Stueben coarsening
    aggregation,            ///< Aggregation
    smoothed_aggregation,   ///< Smoothed aggregation
    smoothed_aggr_emin      ///< Smoothed aggregation with energy minimization
};

inline std::ostream& operator<<(std::ostream &os, type c) {
    switch (c) {
        case ruge_stuben:
            return os << "ruge_stuben";
        case aggregation:
            return os << "aggregation";
        case smoothed_aggregation:
            return os << "smoothed_aggregation";
        case smoothed_aggr_emin:
            return os << "smoothed_aggr_emin";
        default:
            return os << "???";
    }
}

inline std::istream& operator>>(std::istream &in, type &c)
{
    std::string val;
    in >> val;

    if (val == "ruge_stuben")
        c = ruge_stuben;
    else if (val == "aggregation")
        c = aggregation;
    else if (val == "smoothed_aggregation")
        c = smoothed_aggregation;
    else if (val == "smoothed_aggr_emin")
        c = smoothed_aggr_emin;
    else
        throw std::invalid_argument("Invalid coarsening value. Valid choices are: "
                "ruge_stuben, aggregation, smoothed_aggregation, smoothed_aggr_emin.");

    return in;
}

template <class Backend>
struct wrapper {
    typedef boost::property_tree::ptree params;
    type c;
    void *handle;

    wrapper(params prm = params())
        : c(prm.get("type", runtime::coarsening::smoothed_aggregation)),
          handle(0)
    {
        if (!prm.erase("type")) AMGCL_PARAM_MISSING("type");

        switch(c) {

#define AMGCL_RUNTIME_COARSENING(type) \
            case type: \
                handle = call_constructor<amgcl::coarsening::type>(prm); \
                break

            AMGCL_RUNTIME_COARSENING(ruge_stuben);
            AMGCL_RUNTIME_COARSENING(aggregation);
            AMGCL_RUNTIME_COARSENING(smoothed_aggregation);
            AMGCL_RUNTIME_COARSENING(smoothed_aggr_emin);

#undef AMGCL_RUNTIME_COARSENING

            default:
                throw std::invalid_argument("Unsupported coarsening type");
        }
    }

    ~wrapper() {
        switch(c) {

#define AMGCL_RUNTIME_COARSENING(type) \
            case type: \
                call_destructor<amgcl::coarsening::type>(); \
                break

            AMGCL_RUNTIME_COARSENING(ruge_stuben);
            AMGCL_RUNTIME_COARSENING(aggregation);
            AMGCL_RUNTIME_COARSENING(smoothed_aggregation);
            AMGCL_RUNTIME_COARSENING(smoothed_aggr_emin);

#undef AMGCL_RUNTIME_COARSENING
        }
    }

    template <class Matrix>
    std::tuple<
        std::shared_ptr<Matrix>,
        std::shared_ptr<Matrix>
        >
    transfer_operators(const Matrix &A) {
        switch(c) {

#define AMGCL_RUNTIME_COARSENING(type) \
            case type: \
                return make_operators<amgcl::coarsening::type>(A)

            AMGCL_RUNTIME_COARSENING(ruge_stuben);
            AMGCL_RUNTIME_COARSENING(aggregation);
            AMGCL_RUNTIME_COARSENING(smoothed_aggregation);
            AMGCL_RUNTIME_COARSENING(smoothed_aggr_emin);

#undef AMGCL_RUNTIME_COARSENING

            default:
                throw std::invalid_argument("Unsupported coarsening type");
        }
    }

    template <class Matrix>
    std::shared_ptr<Matrix>
    coarse_operator(const Matrix &A, const Matrix &P, const Matrix &R) const {
        switch(c) {

#define AMGCL_RUNTIME_COARSENING(type) \
            case type: \
                return make_coarse<amgcl::coarsening::type>(A, P, R)

            AMGCL_RUNTIME_COARSENING(ruge_stuben);
            AMGCL_RUNTIME_COARSENING(aggregation);
            AMGCL_RUNTIME_COARSENING(smoothed_aggregation);
            AMGCL_RUNTIME_COARSENING(smoothed_aggr_emin);

#undef AMGCL_RUNTIME_COARSENING

            default:
                throw std::invalid_argument("Unsupported coarsening type");
        }
    }

    template <template <class> class Coarsening>
    typename std::enable_if<
        backend::coarsening_is_supported<Backend, Coarsening>::value,
        void*
    >::type
    call_constructor(const params &prm) {
        return static_cast<void*>(new Coarsening<Backend>(prm));
    }

    template <template <class> class Coarsening>
    typename std::enable_if<
        !backend::coarsening_is_supported<Backend, Coarsening>::value,
        void*
    >::type
    call_constructor(const params&) {
        throw std::logic_error("The coarsening is not supported by the backend");
    }

    template <template <class> class Coarsening>
    typename std::enable_if<
        backend::coarsening_is_supported<Backend, Coarsening>::value,
        void
    >::type
    call_destructor() {
        delete static_cast<Coarsening<Backend>*>(handle);
    }

    template <template <class> class Coarsening>
    typename std::enable_if<
        !backend::coarsening_is_supported<Backend, Coarsening>::value,
        void
    >::type
    call_destructor() {
    }

    template <template <class> class Coarsening, class Matrix>
    typename std::enable_if<
        backend::coarsening_is_supported<Backend, Coarsening>::value,
        std::tuple<
            std::shared_ptr<Matrix>,
            std::shared_ptr<Matrix>
            >
    >::type
    make_operators(const Matrix &A) const {
        return static_cast<Coarsening<Backend>*>(handle)->transfer_operators(A);
    }

    template <template <class> class Coarsening, class Matrix>
    typename std::enable_if<
        !backend::coarsening_is_supported<Backend, Coarsening>::value,
        std::tuple<
            std::shared_ptr<Matrix>,
            std::shared_ptr<Matrix>
            >
    >::type
    make_operators(const Matrix&) {
        throw std::logic_error("The coarsening is not supported by the backend");
    }

    template <template <class> class Coarsening, class Matrix>
    typename std::enable_if<
        backend::coarsening_is_supported<Backend, Coarsening>::value,
        std::shared_ptr<Matrix>
    >::type
    make_coarse(const Matrix &A, const Matrix &P, const Matrix &R) const {
        return static_cast<Coarsening<Backend>*>(handle)->coarse_operator(A, P, R);
    }

    template <template <class> class Coarsening, class Matrix>
    typename std::enable_if<
        !backend::coarsening_is_supported<Backend, Coarsening>::value,
        std::shared_ptr<Matrix>
    >::type
    make_coarse(const Matrix&, const Matrix&, const Matrix&) const {
        throw std::logic_error("The coarsening is not supported by the backend");
    }
};

} // namespace coarsening
} // namespace runtime
} // namespace amgcl

#endif
