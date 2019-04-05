#ifndef AMGCL_PERF_COUNTER_MPI_AGGREGATOR_HPP
#define AMGCL_PERF_COUNTER_MPI_AGGREGATOR_HPP

/*
The MIT License

Copyright (c) 2012-2019 Denis Demidov <dennis.demidov@gmail.com>
Copyright (c) 2016 Mohammad Siahatgar <siahatgar@luis.uni-hannover.de>

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
 * \file   amgcl/perf_counter/mpi_aggregator.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Aggregate performace counter over MPI.
 */

#include <functional>
#include <cmath>
#include <type_traits>
#include <amgcl/mpi/util.hpp>

namespace amgcl {
namespace perf_counter {

template <class Counter, bool SingleReaderPerNode = true>
class mpi_aggregator {
    public:
        typedef typename Counter::value_type value_type;

        mpi_aggregator(MPI_Comm comm = MPI_COMM_WORLD)
            : world(comm), dtype(amgcl::mpi::datatype<value_type>())
        {
            if (SingleReaderPerNode) {
                typedef std::integral_constant<bool, sizeof(size_t) == sizeof(int)>::type _32bit;

                char node_name[MPI_MAX_PROCESSOR_NAME];
                int node_name_len, node_master;

                MPI_Get_processor_name(node_name, &node_name_len);
                MPI_Comm_split(world, name_hash(node_name, _32bit()), world.rank, &node_comm);
                MPI_Allreduce(&world.rank, &node_master, 1, MPI_INT, MPI_MIN, node_comm);

                reader = (world.rank == node_master);
                MPI_Comm_split(world, reader, world.rank, &reader_comm);
            }
        }

        ~mpi_aggregator() {
            if (SingleReaderPerNode) {
                MPI_Comm_free(&node_comm);
                MPI_Comm_free(&reader_comm);
            }
        }

        static const char* units() {
            return Counter::units();
        }

        value_type current() {
            value_type gval;

            if (SingleReaderPerNode) {
                if (reader) {
                    value_type lval = counter.current();
                    MPI_Allreduce(&lval, &gval, 1, dtype, MPI_SUM, reader_comm);
                }

                MPI_Bcast(&gval, 1, dtype, 0, node_comm);
            } else {
                value_type lval = counter.current();
                MPI_Allreduce(&lval, &gval, 1, dtype, MPI_SUM, world);
            }

            return gval;
        }
    private:
        amgcl::mpi::communicator world;
        MPI_Datatype dtype;

        bool reader;
        MPI_Comm node_comm, reader_comm;

        Counter counter;

        int name_hash(const char *name, std::true_type) {
            return std::hash<std::string>()(name);
        }

        int name_hash(const char *name, std::false_type) {
            union {
                size_t full;
                struct {
                    int lo, hi;
                } part;
            } h;

            h.full = std::hash<std::string>()(name);
            return std::abs(h.part.lo ^ h.part.hi);
        }
};

} // namespace perf_counter
} // namespace amgcl

#endif
