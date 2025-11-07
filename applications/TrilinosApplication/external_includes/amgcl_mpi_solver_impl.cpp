//  KRATOS  _____     _ _ _
//         |_   _| __(_) (_)_ __   ___  ___
//           | || '__| | | | '_ \ / _ \/ __|
//           | || |  | | | | | | | (_) \__
//           |_||_|  |_|_|_|_| |_|\___/|___/ APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Denis Demidov
//                   Riccardo Rossi
//

// The implementation of AMGCLSolver is split between
// - an implementation header ("linear_solvers/amgcl_solver_impl.hpp")
// - implementation sources (this source file and "linear_solvers/amgcl_solver_impl.cpp")
//
// The reason is twofold:
// - includes from the AMGCL library are extremely heavy, so they are
//   avoided in the class declaration ("linear_solvers/amgcl_solver.h").
//   Instead, the implementation header includes them and defines logic
//   common to any matrix/vector representations. Each source file that
//   defines an instantiation of AMGCLSolver includes the implementation
//   header.
// - Shared memory and distributed memory matrix/vector representations
//   are handled in separate source files to avoid adding a Trilinos
//   dependency to core.

// External includes
#include "amgcl/adapter/epetra.hpp"

// Project includes
#include "trilinos_space.h"
#include "amgcl_mpi_solver.h"
#include "custom_utilities/trilinos_solver_utilities.h"

#define KRATOS_AMGCL_MPI // <= avoid including mpi.h in KratosCore
#include "linear_solvers/amgcl_solver_impl.hpp"
#undef KRATOS_AMGCL_MPI

// System includes
#include <optional>



namespace Kratos {


template <>
struct AMGCLAdaptor<TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>>
{
    template <int BlockSize>
    auto MakeMatrixAdaptor(const Epetra_FECrsMatrix& rMatrix)
    {
        mAdaptor.emplace(amgcl::adapter::map(rMatrix));
        if constexpr (BlockSize == 1) {
            return mAdaptor.value();
        } else {
            using BlockType = amgcl::static_matrix<
                double,
                BlockSize,
                BlockSize
            >;
            return amgcl::adapter::block_matrix<BlockType>(mAdaptor.value());
        }
    }

    template <class TStaticMatrix>
    std::size_t BlockSystemSize(const Epetra_FECrsMatrix& rMatrix) const noexcept
    {
        return rMatrix.RowMap().NumMyElements() / AMGCLStaticVectorTraits<TStaticMatrix>::value;
    }

    auto MakeVectorIterator(const Epetra_FEVector& rVector) const
    {
        return rVector.Values();
    }

    auto MakeVectorIterator(Epetra_FEVector& rVector) const
    {
        return rVector.Values();
    }

    MPI_Comm GetCommunicator(Epetra_FECrsMatrix& rMatrix) const noexcept
    {
        return TrilinosSolverUtilities::GetMPICommFromEpetraComm(rMatrix.Comm());
    }

private:
    // amgcl::adapter::block_matrix constructs a class
    // that stores a reference to the "matrix" passed
    // into it, which in this case means the adaptor
    // defined below. We need to keep it alive until
    // the hierarchy construction finishes, hence the
    // convoluted member variable.
    // Optional is used here to represent an invalid state
    // of the matrix view, before InitializeSolutionStep is
    // called.
    std::optional<amgcl::adapter::epetra_map> mAdaptor;
};


template class KRATOS_API(TRILINOS_APPLICATION) AmgclMPISolver<
    TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>,
    UblasSpace<double, Matrix, Vector>
>;


} // namespace Kratos
