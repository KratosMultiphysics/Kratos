#ifndef AMGCL_PRECONDITIONER_RUNTIME_HPP
#define AMGCL_PRECONDITIONER_RUNTIME_HPP

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
 * \file   amgcl/preconditioner/runtime.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Runtime-configurable wrappers around amgcl classes.
 */

#ifdef AMGCL_NO_BOOST
#  error Runtime interface relies on Boost.PropertyTree!
#endif

#include <boost/property_tree/ptree.hpp>

#include <amgcl/util.hpp>
#include <amgcl/solver/runtime.hpp>
#include <amgcl/coarsening/runtime.hpp>
#include <amgcl/relaxation/runtime.hpp>
#include <amgcl/relaxation/as_preconditioner.hpp>
#include <amgcl/preconditioner/dummy.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/amg.hpp>

namespace amgcl {
namespace runtime {

/// Preconditioner kinds.
namespace precond_class {
enum type {
    amg,            ///< AMG
    relaxation,     ///< Single-level relaxation
    dummy,          ///< Identity matrix as preconditioner.
    nested          ///< Nested solver as preconditioner.
};

inline std::ostream& operator<<(std::ostream &os, type p) {
    switch (p) {
        case amg:
            return os << "amg";
        case relaxation:
            return os << "relaxation";
        case dummy:
            return os << "dummy";
        case nested:
            return os << "nested";
        default:
            return os << "???";
    }
}

inline std::istream& operator>>(std::istream &in, type &p)
{
    std::string val;
    in >> val;

    if (val == "amg")
        p = amg;
    else if (val == "relaxation")
        p = relaxation;
    else if (val == "dummy")
        p = dummy;
    else if (val == "nested")
        p = nested;
    else
        throw std::invalid_argument("Invalid preconditioner class. Valid choices are: "
                "amg, relaxation, dummy, nested");

    return in;
}
} // namespace precond_class

template <class Backend>
class preconditioner {
    public:
        typedef Backend backend_type;

        typedef typename Backend::value_type value_type;
        typedef typename Backend::matrix     matrix;
        typedef typename Backend::vector     vector;
        typedef typename Backend::params     backend_params;

        typedef boost::property_tree::ptree params;

        template <class Matrix>
        preconditioner(
                const Matrix &A,
                params prm = params(),
                const backend_params &bprm = backend_params())
            : _class(prm.get("class", runtime::precond_class::amg)),
              handle(0)
        {
            if (!prm.erase("class")) AMGCL_PARAM_MISSING("class");

            switch(_class) {
                case precond_class::amg:
                    {
                        typedef
                            amgcl::amg<Backend, runtime::coarsening::wrapper, runtime::relaxation::wrapper>
                            Precond;

                        handle = static_cast<void*>(new Precond(A, prm, bprm));
                    }
                    break;
                case precond_class::relaxation:
                    {
                        typedef
                            amgcl::relaxation::as_preconditioner<Backend, runtime::relaxation::wrapper>
                            Precond;

                        handle = static_cast<void*>(new Precond(A, prm, bprm));
                    }
                    break;
                case precond_class::dummy:
                    {
                        typedef
                            amgcl::preconditioner::dummy<Backend>
                            Precond;

                        handle = static_cast<void*>(new Precond(A, prm, bprm));
                    }
                    break;
                case precond_class::nested:
                    {
                        typedef
                            make_solver<
                                preconditioner,
                                runtime::solver::wrapper<Backend>
                                >
                            Precond;

                        handle = static_cast<void*>(new Precond(A, prm, bprm));
                    }
                    break;
                default:
                    throw std::invalid_argument("Unsupported preconditioner class");
            }
        }

        ~preconditioner() {
            switch (_class) {
                case precond_class::amg:
                    {
                        typedef
                            amgcl::amg<Backend, runtime::coarsening::wrapper, runtime::relaxation::wrapper>
                            Precond;

                        delete static_cast<Precond*>(handle);
                    }
                    break;
                case precond_class::relaxation:
                    {
                        typedef
                            amgcl::relaxation::as_preconditioner<Backend, runtime::relaxation::wrapper>
                            Precond;

                        delete static_cast<Precond*>(handle);
                    }
                    break;
                case precond_class::dummy:
                    {
                        typedef
                            amgcl::preconditioner::dummy<Backend>
                            Precond;

                        delete static_cast<Precond*>(handle);
                    }
                    break;
                case precond_class::nested:
                    {
                        typedef
                            make_solver<
                                preconditioner,
                                runtime::solver::wrapper<Backend>
                                >
                            Precond;

                        delete static_cast<Precond*>(handle);
                    }
                    break;
                default:
                    break;
            }
        }

        template <class Vec1, class Vec2>
        void apply(const Vec1 &rhs, Vec2 &x) const {
            switch(_class) {
                case precond_class::amg:
                    {
                        typedef
                            amgcl::amg<Backend, runtime::coarsening::wrapper, runtime::relaxation::wrapper>
                            Precond;

                        static_cast<Precond*>(handle)->apply(rhs, x);
                    }
                    break;
                case precond_class::relaxation:
                    {
                        typedef
                            amgcl::relaxation::as_preconditioner<Backend, runtime::relaxation::wrapper>
                            Precond;

                        static_cast<Precond*>(handle)->apply(rhs, x);
                    }
                    break;
                case precond_class::dummy:
                    {
                        typedef
                            amgcl::preconditioner::dummy<Backend>
                            Precond;

                        static_cast<Precond*>(handle)->apply(rhs, x);
                    }
                    break;
                case precond_class::nested:
                    {
                        typedef
                            make_solver<
                                preconditioner,
                                runtime::solver::wrapper<Backend>
                                >
                            Precond;

                        static_cast<Precond*>(handle)->apply(rhs, x);
                    }
                    break;
                default:
                    throw std::invalid_argument("Unsupported preconditioner class");
            }
        }

        std::shared_ptr<matrix> system_matrix_ptr() const {
            switch(_class) {
                case precond_class::amg:
                    {
                        typedef
                            amgcl::amg<Backend, runtime::coarsening::wrapper, runtime::relaxation::wrapper>
                            Precond;

                        return static_cast<Precond*>(handle)->system_matrix_ptr();
                    }
                case precond_class::relaxation:
                    {
                        typedef
                            amgcl::relaxation::as_preconditioner<Backend, runtime::relaxation::wrapper>
                            Precond;

                        return static_cast<Precond*>(handle)->system_matrix_ptr();
                    }
                case precond_class::dummy:
                    {
                        typedef
                            amgcl::preconditioner::dummy<Backend>
                            Precond;

                        return static_cast<Precond*>(handle)->system_matrix_ptr();
                    }
                case precond_class::nested:
                    {
                        typedef
                            make_solver<
                                preconditioner,
                                runtime::solver::wrapper<Backend>
                                >
                            Precond;

                        return static_cast<Precond*>(handle)->system_matrix_ptr();
                    }
                default:
                    throw std::invalid_argument("Unsupported preconditioner class");
            }
        }

        const matrix& system_matrix() const {
            return *system_matrix_ptr();
        }

        size_t size() const {
            return backend::rows( system_matrix() );
        }

        friend std::ostream& operator<<(std::ostream &os, const preconditioner &p)
        {
            switch(p._class) {
                case precond_class::amg:
                    {
                        typedef
                            amgcl::amg<Backend, runtime::coarsening::wrapper, runtime::relaxation::wrapper>
                            Precond;

                        return os << *static_cast<Precond*>(p.handle);
                    }
                case precond_class::relaxation:
                    {
                        typedef
                            amgcl::relaxation::as_preconditioner<Backend, runtime::relaxation::wrapper>
                            Precond;

                        return os << *static_cast<Precond*>(p.handle);
                    }
                case precond_class::dummy:
                    {
                        typedef
                            amgcl::preconditioner::dummy<Backend>
                            Precond;

                        return os << *static_cast<Precond*>(p.handle);
                    }
                case precond_class::nested:
                    {
                        typedef
                            make_solver<
                                preconditioner,
                                runtime::solver::wrapper<Backend>
                                >
                            Precond;

                        return os << *static_cast<Precond*>(p.handle);
                    }
                default:
                    throw std::invalid_argument("Unsupported preconditioner class");
            }
        }
    private:
        const runtime::precond_class::type _class;

        void *handle;
};

} // namespace runtime
} // namespace amgcl


#endif
