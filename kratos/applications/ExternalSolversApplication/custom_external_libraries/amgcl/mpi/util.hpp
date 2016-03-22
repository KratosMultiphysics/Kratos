#ifndef AMGCL_MPI_UTIL_HPP
#define AMGCL_MPI_UTIL_HPP

/*
The MIT License

Copyright (c) 2012-2016 Denis Demidov <dennis.demidov@gmail.com>
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

#include <mpi.h>
#include <boost/type_traits.hpp>

namespace amgcl {
namespace mpi {

/// Converts C type to MPI datatype.
template <class T, class Enable = void>
struct datatype;

template <>
struct datatype<float> {
    static MPI_Datatype get() { return MPI_FLOAT; }
};

template <>
struct datatype<double> {
    static MPI_Datatype get() { return MPI_DOUBLE; }
};

template <>
struct datatype<long double> {
    static MPI_Datatype get() { return MPI_LONG_DOUBLE; }
};

template <>
struct datatype<int> {
    static MPI_Datatype get() { return MPI_INT; }
};

template <>
struct datatype<long long> {
    static MPI_Datatype get() { return MPI_LONG_LONG_INT; }
};

template <>
struct datatype<ptrdiff_t>
    : boost::conditional<
        sizeof(ptrdiff_t) == sizeof(int), datatype<int>, datatype<long long>
        >::type
{};

/// Convenience wrapper around MPI_Comm.
struct communicator {
    MPI_Comm comm;
    int      rank;
    int      size;

    communicator(MPI_Comm comm) : comm(comm) {
        MPI_Comm_rank(comm, &rank);
        MPI_Comm_size(comm, &size);
    };

    operator MPI_Comm() const {
        return comm;
    }
};

/// Communicator-wise condition checking.
/**
 * Checks conditions at each process in the communicator;
 *
 * If the condition is false on any of the participating processes, outputs the
 * provided message together with the ranks of the offending process.
 * After that each process in the communicator throws.
 */
template <class Condition, class Message>
void precondition(communicator comm, const Condition &cond, const Message &message)
{
    int gc, lc = static_cast<int>(cond);
    MPI_Allreduce(&lc, &gc, 1, MPI_INT, MPI_PROD, comm);

    if (!gc) {
        std::vector<int> c(comm.size);
        MPI_Gather(&lc, 1, MPI_INT, &c[0], comm.size, MPI_INT, 0, comm);
        if (comm.rank == 0) {
            std::cerr << "Failed assumption: " << message << std::endl;
            std::cerr << "Offending processes:";
            for (int i = 0; i < comm.size; ++i)
                if (!c[i]) std::cerr << " " << i;
            std::cerr << std::endl;
        }
        MPI_Barrier(comm);
        throw std::runtime_error(message);
    }
}

} // namespace mpi
} // namespace amgcl

#endif
