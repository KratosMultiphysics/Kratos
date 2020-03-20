//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

#include "add_trilinos_utilities_to_python.h"

// Trilinos includes
#include "Epetra_FEVector.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_MpiComm.h"

// KratosCore dependencies
#include "containers/variable.h"
#include "processes/process.h"
#include "solving_strategies/schemes/scheme.h"
#include "spaces/ublas_space.h"

// TrilinosApplication dependencies
#include "trilinos_space.h"

// FluidDynamicsApplication dependencies
#include "custom_strategies/strategies/gear_scheme.h"
#include "custom_strategies/strategies/residualbased_predictorcorrector_velocity_bossak_scheme_turbulent.h"
#include "custom_strategies/strategies/residualbased_predictorcorrector_velocity_bdf_scheme_turbulent.h"
#include "custom_strategies/strategies/residualbased_simple_steady_scheme.h"

namespace Kratos {
namespace Python {

void AddTrilinosSchemesToPython(pybind11::module& m)
{
    // namespace py = pybind11;

    // using TrilinosSparseSpace = TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>;
    // using UblasLocalSpace = UblasSpace<double, Matrix, Vector>;
}

}
}
