#ifndef AMGCL_MPI_PRECONDITIONER_HPP
#define AMGCL_MPI_PRECONDITIONER_HPP

/*
The MIT License

Copyright (c) 2012-2019 Denis Demidov <dennis.demidov@gmail.com>

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
 * \file   amgcl/mpi/preconditioner.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Runtime wrapper around mpi preconditioners.
 */

#ifdef AMGCL_NO_BOOST
#  error Runtime interface relies on Boost.PropertyTree!
#endif

#include <iostream>

#include <boost/property_tree/ptree.hpp>

#include <amgcl/mpi/amg.hpp>
#include <amgcl/mpi/coarsening/runtime.hpp>
#include <amgcl/mpi/relaxation/runtime.hpp>
#include <amgcl/mpi/direct_solver/runtime.hpp>
#include <amgcl/mpi/partition/runtime.hpp>
#include <amgcl/mpi/relaxation/as_preconditioner.hpp>
#include <amgcl/mpi/distributed_matrix.hpp>
#include <amgcl/mpi/util.hpp>

namespace amgcl {
namespace runtime {
namespace mpi {

/// Preconditioner kinds.
namespace precond_class {
enum type {
    amg,            ///< AMG
    relaxation      ///< Single-level relaxation
};

inline std::ostream& operator<<(std::ostream &os, type p) {
    switch (p) {
        case amg:
            return os << "amg";
        case relaxation:
            return os << "relaxation";
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
    else
        throw std::invalid_argument("Invalid preconditioner class. "
                "Valid choices are: amg, relaxation");

    return in;
}
} // namespace precond_class

template <class Backend>
class preconditioner {
    public:
        typedef Backend backend_type;
        typedef typename backend_type::params backend_params;
        typedef boost::property_tree::ptree params;
        typedef typename backend_type::value_type value_type;
        typedef amgcl::mpi::distributed_matrix<backend_type> matrix;

        template <class Matrix>
        preconditioner(
                amgcl::mpi::communicator comm,
                const Matrix &Astrip,
                params prm = params(),
                const backend_params &bprm = backend_params()
                ) : _class(prm.get("class", precond_class::amg)), handle(0)
        {
            init(std::make_shared<matrix>(comm, Astrip, backend::rows(Astrip)), prm, bprm);
        }

        preconditioner(
                amgcl::mpi::communicator,
                std::shared_ptr<matrix> A,
                params prm = params(),
                const backend_params &bprm = backend_params()
                ) : _class(prm.get("class", precond_class::amg)), handle(0)
        {
            init(A, prm, bprm);
        }

        ~preconditioner() {
            switch (_class) {
                case precond_class::amg:
                    {
                        typedef
                            amgcl::mpi::amg<
                                Backend,
                                amgcl::runtime::mpi::coarsening::wrapper<Backend>,
                                amgcl::runtime::mpi::relaxation::wrapper<Backend>,
                                amgcl::runtime::mpi::direct::solver<value_type>,
                                amgcl::runtime::mpi::partition::wrapper<Backend>
                                >
                            Precond;

                        delete static_cast<Precond*>(handle);
                    }
                    break;
                case precond_class::relaxation:
                    {
                        typedef
                            amgcl::mpi::relaxation::as_preconditioner<
                                amgcl::runtime::mpi::relaxation::wrapper<Backend>
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
        void apply(const Vec1 &rhs, Vec2 &&x) const {
            switch(_class) {
                case precond_class::amg:
                    {
                        typedef
                            amgcl::mpi::amg<
                                Backend,
                                amgcl::runtime::mpi::coarsening::wrapper<Backend>,
                                amgcl::runtime::mpi::relaxation::wrapper<Backend>,
                                amgcl::runtime::mpi::direct::solver<value_type>,
                                amgcl::runtime::mpi::partition::wrapper<Backend>
                                >
                            Precond;

                        static_cast<Precond*>(handle)->apply(rhs, x);
                    }
                    break;
                case precond_class::relaxation:
                    {
                        typedef
                            amgcl::mpi::relaxation::as_preconditioner<
                                amgcl::runtime::mpi::relaxation::wrapper<Backend>
                                >
                            Precond;

                        static_cast<Precond*>(handle)->apply(rhs, x);
                    }
                    break;
                default:
                    throw std::invalid_argument("Unsupported preconditioner class");
            }
        }

        /// Returns the system matrix from the finest level.
        std::shared_ptr<matrix> system_matrix_ptr() const {
            switch(_class) {
                case precond_class::amg:
                    {
                        typedef
                            amgcl::mpi::amg<
                                Backend,
                                amgcl::runtime::mpi::coarsening::wrapper<Backend>,
                                amgcl::runtime::mpi::relaxation::wrapper<Backend>,
                                amgcl::runtime::mpi::direct::solver<value_type>,
                                amgcl::runtime::mpi::partition::wrapper<Backend>
                                >
                            Precond;

                        return static_cast<Precond*>(handle)->system_matrix_ptr();
                    }
                case precond_class::relaxation:
                    {
                        typedef
                            amgcl::mpi::relaxation::as_preconditioner<
                                amgcl::runtime::mpi::relaxation::wrapper<Backend>
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

        friend std::ostream& operator<<(std::ostream &os, const preconditioner &p) {
            switch(p._class) {
                case precond_class::amg:
                    {
                        typedef
                            amgcl::mpi::amg<
                                Backend,
                                amgcl::runtime::mpi::coarsening::wrapper<Backend>,
                                amgcl::runtime::mpi::relaxation::wrapper<Backend>,
                                amgcl::runtime::mpi::direct::solver<value_type>,
                                amgcl::runtime::mpi::partition::wrapper<Backend>
                                >
                            Precond;

                        return os << *static_cast<Precond*>(p.handle);
                    }
                case precond_class::relaxation:
                    {
                        typedef
                            amgcl::mpi::relaxation::as_preconditioner<
                                amgcl::runtime::mpi::relaxation::wrapper<Backend>
                                >
                            Precond;

                        return os << *static_cast<Precond*>(p.handle);
                    }
                default:
                    throw std::invalid_argument("Unsupported preconditioner class");
            }
        }

    private:
        precond_class::type _class;
        void *handle;

        void init(std::shared_ptr<matrix> A, params &prm, const backend_params &bprm) {
            if (!prm.erase("class")) AMGCL_PARAM_MISSING("class");

            switch(_class) {
                case precond_class::amg:
                    {
                        typedef
                            amgcl::mpi::amg<
                                Backend,
                                amgcl::runtime::mpi::coarsening::wrapper<Backend>,
                                amgcl::runtime::mpi::relaxation::wrapper<Backend>,
                                amgcl::runtime::mpi::direct::solver<value_type>,
                                amgcl::runtime::mpi::partition::wrapper<Backend>
                                >
                            Precond;

                        handle = static_cast<void*>(new Precond(A->comm(), A, prm, bprm));
                    }
                    break;
                case precond_class::relaxation:
                    {
                        typedef
                            amgcl::mpi::relaxation::as_preconditioner<
                                amgcl::runtime::mpi::relaxation::wrapper<Backend>
                                >
                            Precond;

                        handle = static_cast<void*>(new Precond(A->comm(), A, prm, bprm));
                    }
                    break;
                default:
                    throw std::invalid_argument("Unsupported preconditioner class");
            }
        }
};

} // namespace mpi
} // namespace runtime
} // namespace amgcl

#endif
