#ifndef AMGCL_MPI_UTIL_HPP
#define AMGCL_MPI_UTIL_HPP

/*
The MIT License

Copyright (c) 2012-2018 Denis Demidov <dennis.demidov@gmail.com>
Copyright (c) 2014, Riccardo Rossi, CIMNE (International Center for Numerical Methods in Engineering)

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
 * \file   amgcl/mpi/util.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  MPI utilities.
 */

#include <vector>
#include <numeric>

#include <type_traits>
#include <amgcl/value_type/interface.hpp>

#include <mpi.h>

namespace amgcl {
namespace mpi {

/// Converts C type to MPI datatype.
template <class T, class Enable = void>
struct datatype_impl {
    static MPI_Datatype get() {
        static const MPI_Datatype t = create();
        return t;
    }

    static MPI_Datatype create() {
        typedef typename math::scalar_of<T>::type S;
        MPI_Datatype t;
        int n = sizeof(T) / sizeof(S);
        MPI_Type_contiguous(n, datatype_impl<S>::get(), &t);
        MPI_Type_commit(&t);
        return t;
    }
};

template <>
struct datatype_impl<float> {
    static MPI_Datatype get() { return MPI_FLOAT; }
};

template <>
struct datatype_impl<double> {
    static MPI_Datatype get() { return MPI_DOUBLE; }
};

template <>
struct datatype_impl<long double> {
    static MPI_Datatype get() { return MPI_LONG_DOUBLE; }
};

template <>
struct datatype_impl<int> {
    static MPI_Datatype get() { return MPI_INT; }
};

template <>
struct datatype_impl<unsigned> {
    static MPI_Datatype get() { return MPI_UNSIGNED; }
};

template <>
struct datatype_impl<long long> {
    static MPI_Datatype get() { return MPI_LONG_LONG_INT; }
};

template <>
struct datatype_impl<unsigned long long> {
    static MPI_Datatype get() { return MPI_UNSIGNED_LONG_LONG; }
};

template <>
struct datatype_impl<ptrdiff_t>
    : std::conditional<
        sizeof(ptrdiff_t) == sizeof(int), datatype_impl<int>, datatype_impl<long long>
        >::type
{};

template <>
struct datatype_impl<size_t>
    : std::conditional<
        sizeof(size_t) == sizeof(unsigned), datatype_impl<unsigned>, datatype_impl<unsigned long long>
        >::type
{};

template <>
struct datatype_impl<char> {
    static MPI_Datatype get() { return MPI_CHAR; }
};

template <typename T>
MPI_Datatype datatype() {
    return datatype_impl<T>::get();
}

/// Convenience wrapper around MPI_Comm.
struct communicator {
    MPI_Comm comm;
    int      rank;
    int      size;

    communicator() {}

    communicator(MPI_Comm comm) : comm(comm) {
        MPI_Comm_rank(comm, &rank);
        MPI_Comm_size(comm, &size);
    };

    operator MPI_Comm() const {
        return comm;
    }

    /// Exclusive sum over mpi communicator
    template <typename T>
    std::vector<T> exclusive_sum(T n) const {
        std::vector<T> v(size + 1); v[0] = 0;
        MPI_Allgather(&n, 1, datatype<T>(), &v[1], 1, datatype<T>(), comm);
        std::partial_sum(v.begin(), v.end(), v.begin());
        return v;
    }

    template <typename T>
    T reduce(MPI_Op op, const T &lval) const {
        typedef typename math::scalar_of<T>::type S;

        const int elems = sizeof(T) / sizeof(S);
        T gval;

        MPI_Allreduce((void*)&lval, &gval, elems, datatype<T>(), op, comm);
        return gval;
    }

    /// Communicator-wise condition checking.
    /**
     * Checks conditions at each process in the communicator;
     *
     * If the condition is false on any of the participating processes, outputs the
     * provided message together with the ranks of the offending process.
     * After that each process in the communicator throws.
     */
    template <class Condition, class Message>
    void check(const Condition &cond, const Message &message) {
        int lc = static_cast<int>(cond);
        int gc = reduce(MPI_PROD, lc);

        if (!gc) {
            std::vector<int> c(size);
            MPI_Gather(&lc, 1, MPI_INT, &c[0], size, MPI_INT, 0, comm);
            if (rank == 0) {
                std::cerr << "Failed assumption: " << message << std::endl;
                std::cerr << "Offending processes:";
                for (int i = 0; i < size; ++i)
                    if (!c[i]) std::cerr << " " << i;
                std::cerr << std::endl;
            }
            MPI_Barrier(comm);
            throw std::runtime_error(message);
        }
    }

};

} // namespace mpi
} // namespace amgcl

#endif
