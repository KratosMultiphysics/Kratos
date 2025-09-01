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



namespace Kratos {


template <>
struct AMGCLAdaptor<TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>>
{
    template <int BlockSize>
    auto MakeMatrixAdaptor(const Epetra_FECrsMatrix& rMatrix) const
    {
        if constexpr (BlockSize == 1) {
            return amgcl::adapter::map(rMatrix);
        } else {
            using BlockType = amgcl::static_matrix<
                double,
                BlockSize,
                BlockSize
            >;
            return amgcl::adapter::block_matrix<BlockType>(amgcl::adapter::map(rMatrix));
        }
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
};


template class KRATOS_API(TRILINOS_APPLICATION) AmgclMPISolver<
    TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>,
    UblasSpace<double, Matrix, Vector>
>;


} // namespace Kratos
