#ifndef KRATOS_AMGCL_MPI_SOLVE_FUNCTIONS_H
#define KRATOS_AMGCL_MPI_SOLVE_FUNCTIONS_H

#include <boost/range/iterator_range.hpp>

#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/adapter/epetra.hpp>
#include <amgcl/adapter/ublas.hpp>
#include <amgcl/adapter/zero_copy.hpp>
#include <amgcl/adapter/block_matrix.hpp>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/value_type/static_matrix.hpp>
#include <amgcl/solver/runtime.hpp>

#include <amgcl/mpi/util.hpp>
#include <amgcl/mpi/make_solver.hpp>
#include <amgcl/mpi/amg.hpp>
#include <amgcl/mpi/coarsening/runtime.hpp>
#include <amgcl/mpi/relaxation/runtime.hpp>
#include <amgcl/mpi/direct_solver/runtime.hpp>
#include <amgcl/mpi/partition/runtime.hpp>

namespace Kratos
{

// Spacialization of AMGCLScalarSolve for distribued systems.
template <class TSparseSpaceType>
typename std::enable_if<TSparseSpaceType::IsDistributed(), void>::type
AMGCLScalarSolve(
    typename TSparseSpaceType::MatrixType& rA,
    typename TSparseSpaceType::VectorType& rX,
    typename TSparseSpaceType::VectorType& rB,
    typename TSparseSpaceType::IndexType& rIterationNumber,
    double& rResidual,
    const boost::property_tree::ptree &amgclParams,
    int verbosity_level
    )
{
    typedef amgcl::backend::builtin<double> Backend;

    typedef
        amgcl::mpi::make_solver<
            amgcl::mpi::amg<
                Backend,
                amgcl::runtime::mpi::coarsening::wrapper<Backend>,
                amgcl::runtime::mpi::relaxation::wrapper<Backend>,
                amgcl::runtime::mpi::direct::solver<double>,
                amgcl::runtime::mpi::partition::wrapper<Backend>
                >,
            amgcl::runtime::solver::wrapper
            >
        Solver;

    Solver solve(MPI_COMM_WORLD, amgcl::adapter::map(rA), amgclParams);

    std::size_t n = rA.NumMyRows();

    auto b_range = boost::make_iterator_range(rB.Values(), rB.Values() + n);
    auto x_range = boost::make_iterator_range(rX.Values(), rX.Values() + n);

    std::tie(rIterationNumber, rResidual) = solve(b_range, x_range);
}

// Spacialization of AMGCLBlockSolve for distribued systems.
template <int TBlockSize, class TSparseSpaceType>
typename std::enable_if<TSparseSpaceType::IsDistributed(), void>::type
AMGCLBlockSolve(
    typename TSparseSpaceType::MatrixType & rA,
    typename TSparseSpaceType::VectorType& rX,
    typename TSparseSpaceType::VectorType& rB,
    typename TSparseSpaceType::IndexType& rIterationNumber,
    double& rResidual,
    boost::property_tree::ptree amgclParams,
    int verbosity_level
    )
{
    amgclParams.put("precond.coarsening.aggr.block_size",1);

    typedef amgcl::static_matrix<double, TBlockSize, TBlockSize> val_type;
    typedef amgcl::static_matrix<double, TBlockSize, 1> rhs_type;
    typedef amgcl::backend::builtin<val_type> Backend;

    std::size_t n = rA.RowMap().NumMyElements();

    typedef
        amgcl::mpi::make_solver<
            amgcl::mpi::amg<
                Backend,
                amgcl::runtime::mpi::coarsening::wrapper<Backend>,
                amgcl::runtime::mpi::relaxation::wrapper<Backend>,
                amgcl::runtime::mpi::direct::solver<val_type>,
                amgcl::runtime::mpi::partition::wrapper<Backend>
                >,
            amgcl::runtime::solver::wrapper
            >
        Solver;

    Solver solve(
            MPI_COMM_WORLD,
            amgcl::adapter::block_matrix<val_type>(amgcl::adapter::map(rA)),
            amgclParams
            );

    auto b_begin = reinterpret_cast<const rhs_type*>(rB.Values());
    auto x_begin = reinterpret_cast<rhs_type*>(rX.Values());

    auto b_range = boost::make_iterator_range(b_begin, b_begin + n / TBlockSize);
    auto x_range = boost::make_iterator_range(x_begin, x_begin + n / TBlockSize);

    std::tie(rIterationNumber, rResidual) = solve(b_range, x_range);
}

} // namespace Kratos

#endif
