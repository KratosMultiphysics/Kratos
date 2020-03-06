#include <boost/range/iterator_range.hpp>
#include <boost/property_tree/ptree.hpp>

#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/adapter/epetra.hpp>
#include <amgcl/adapter/ublas.hpp>
#include <amgcl/adapter/zero_copy.hpp>
#include <amgcl/adapter/block_matrix.hpp>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/value_type/static_matrix.hpp>

#include <amgcl/mpi/util.hpp>
#include <amgcl/mpi/make_solver.hpp>
#include <amgcl/mpi/preconditioner.hpp>
#include <amgcl/mpi/solver/runtime.hpp>

#ifdef AMGCL_GPGPU
#  include <amgcl/backend/vexcl.hpp>
#  include <amgcl/backend/vexcl_static_matrix.hpp>
#endif

#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"
#include "trilinos_space.h"

namespace Kratos
{

#ifdef AMGCL_GPGPU
vex::Context& vexcl_context();

template <int TBlockSize>
void register_vexcl_static_matrix_type();
#endif

// Spacialization of AMGCLScalarSolve for distribued systems.
void AMGCLScalarSolve(
    TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>::MatrixType& rA,
    TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>::VectorType& rX,
    TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>::VectorType& rB,
    TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>::IndexType& rIterationNumber,
    double& rResidual,
    const boost::property_tree::ptree &amgclParams,
    int verbosity_level,
    bool use_gpgpu
    )
{
#ifdef AMGCL_GPGPU
    if (use_gpgpu && vexcl_context()) {
        auto &ctx = vexcl_context();

        typedef amgcl::backend::vexcl<double> Backend;

        typedef
            amgcl::mpi::make_solver<
                amgcl::runtime::mpi::preconditioner<Backend>,
                amgcl::runtime::mpi::solver::wrapper<Backend>
                >
            Solver;

        Backend::params bprm;
        bprm.q = ctx;

        Solver solve(MPI_COMM_WORLD, amgcl::adapter::map(rA), amgclParams, bprm);

        std::size_t n = rA.NumMyRows();

        vex::vector<double> b(ctx, n, rB.Values());
        vex::vector<double> x(ctx, n, rX.Values());

        std::tie(rIterationNumber, rResidual) = solve(b, x);

        vex::copy(x.begin(), x.end(), rX.Values());
    } else
#endif
    {
        typedef amgcl::backend::builtin<double> Backend;

        typedef
            amgcl::mpi::make_solver<
                amgcl::runtime::mpi::preconditioner<Backend>,
                amgcl::runtime::mpi::solver::wrapper<Backend>
                >
            Solver;

        Solver solve(MPI_COMM_WORLD, amgcl::adapter::map(rA), amgclParams);

        std::size_t n = rA.NumMyRows();

        auto b_range = boost::make_iterator_range(rB.Values(), rB.Values() + n);
        auto x_range = boost::make_iterator_range(rX.Values(), rX.Values() + n);

        std::tie(rIterationNumber, rResidual) = solve(b_range, x_range);
    }
}

// Spacialization of AMGCLBlockSolve for distribued systems.
template <int TBlockSize>
void AMGCLBlockSolve(
    TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>::MatrixType & rA,
    TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>::VectorType& rX,
    TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>::VectorType& rB,
    TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>::IndexType& rIterationNumber,
    double& rResidual,
    boost::property_tree::ptree amgclParams,
    int verbosity_level,
    bool use_gpgpu
    )
{
    if(amgclParams.get<std::string>("precond.class") != "amg")
        amgclParams.erase("precond.coarsening");
    else
        amgclParams.put("precond.coarsening.aggr.block_size",1);

    typedef amgcl::static_matrix<double, TBlockSize, TBlockSize> val_type;
    typedef amgcl::static_matrix<double, TBlockSize, 1> rhs_type;

    std::size_t n = rA.RowMap().NumMyElements();
    std::size_t nb = n / TBlockSize;

#ifdef AMGCL_GPGPU
    if (use_gpgpu && vexcl_context()) {
        auto &ctx = vexcl_context();
        register_vexcl_static_matrix_type<TBlockSize>();

        typedef amgcl::backend::vexcl<val_type> Backend;

        typedef
            amgcl::mpi::make_solver<
                amgcl::runtime::mpi::preconditioner<Backend>,
                amgcl::runtime::mpi::solver::wrapper<Backend>
                >
            Solver;

        typename Backend::params bprm;
        bprm.q = ctx;

        Solver solve(
                MPI_COMM_WORLD,
                amgcl::adapter::block_matrix<val_type>(amgcl::adapter::map(rA)),
                amgclParams, bprm
                );

        auto b_begin = reinterpret_cast<const rhs_type*>(rB.Values());
        auto x_begin = reinterpret_cast<rhs_type*>(rX.Values());

        vex::vector<rhs_type> x(ctx, nb, x_begin);
        vex::vector<rhs_type> b(ctx, nb, b_begin);

        std::tie(rIterationNumber, rResidual) = solve(b, x);

        vex::copy(x.begin(), x.end(), x_begin);
    } else
#endif
    {
        typedef amgcl::backend::builtin<val_type> Backend;

        typedef
            amgcl::mpi::make_solver<
                amgcl::runtime::mpi::preconditioner<Backend>,
                amgcl::runtime::mpi::solver::wrapper<Backend>
                >
            Solver;

        Solver solve(
                MPI_COMM_WORLD,
                amgcl::adapter::block_matrix<val_type>(amgcl::adapter::map(rA)),
                amgclParams
                );

        auto b_begin = reinterpret_cast<const rhs_type*>(rB.Values());
        auto x_begin = reinterpret_cast<rhs_type*>(rX.Values());

        auto b_range = boost::make_iterator_range(b_begin, b_begin + nb);
        auto x_range = boost::make_iterator_range(x_begin, x_begin + nb);

        std::tie(rIterationNumber, rResidual) = solve(b_range, x_range);
    }
}

void AMGCLSolve(
    int block_size,
    TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>::MatrixType& rA,
    TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>::VectorType& rX,
    TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>::VectorType& rB,
    TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>::IndexType& rIterationNumber,
    double& rResidual,
    boost::property_tree::ptree amgclParams,
    int verbosity_level,
    bool use_gpgpu
    )
{
    switch (block_size) {
        case 2:
            AMGCLBlockSolve<2>(rA, rX, rB, rIterationNumber, rResidual, amgclParams, verbosity_level, use_gpgpu);
            return;
        case 3:
            AMGCLBlockSolve<3>(rA, rX, rB, rIterationNumber, rResidual, amgclParams, verbosity_level, use_gpgpu);
            return;
        case 4:
            AMGCLBlockSolve<4>(rA, rX, rB, rIterationNumber, rResidual, amgclParams, verbosity_level, use_gpgpu);
            return;
        default:
            AMGCLScalarSolve(rA, rX, rB, rIterationNumber, rResidual, amgclParams, verbosity_level, use_gpgpu);
            return;
    }
}

} // namespace Kratos
