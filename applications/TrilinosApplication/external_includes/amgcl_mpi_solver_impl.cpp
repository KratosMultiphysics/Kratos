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

// External includes
#include "amgcl/adapter/epetra.hpp"

// Project includes
#include "trilinos_space.h"
#include "amgcl_mpi_solver.h"
#include "custom_utilities/trilinos_solver_utilities.h"

#define KRATOS_AMGCL_MPI
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
    std::optional<amgcl::adapter::epetra_map> mAdaptor;
};


template class KRATOS_API(TRILINOS_APPLICATION) AmgclMPISolver<
    TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>,
    UblasSpace<double, Matrix, Vector>
>;


} // namespace Kratos
