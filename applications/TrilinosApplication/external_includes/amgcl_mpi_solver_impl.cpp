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
#include "linear_solvers/amgcl_solver_impl.hpp"
#include "amgcl_mpi_solver.h"


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
            using BackendMatrix = amgcl::static_matrix<
                double,
                BlockSize,
                BlockSize
            >;
            return amgcl::adapter::block_matrix<BackendMatrix>(amgcl::adapter::map(rMatrix));
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
};


template class AmgclMPISolver<
    TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>,
    UblasSpace<double, Matrix, Vector>
>;


} // namespace Kratos
