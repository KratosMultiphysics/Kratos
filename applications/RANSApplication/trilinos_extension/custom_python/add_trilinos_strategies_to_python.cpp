//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// Trilinos includes
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"
#include "Epetra_MpiComm.h"

// KratosCore dependencies
#include "includes/model_part.h"
#include "spaces/ublas_space.h"

// TrilinosApplication dependencies
#include "trilinos_space.h"

// RANS trilinos extensions
// schemes
#include "custom_strategies/algebraic_flux_corrected_steady_scalar_scheme.h"
#include "custom_strategies/bossak_relaxation_scalar_scheme.h"
#include "custom_strategies/steady_scalar_scheme.h"

// Include base h
#include "add_trilinos_strategies_to_python.h"

namespace Kratos
{
namespace Python
{
void AddTrilinosStrategiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    using LocalSpaceType = UblasSpace<double, Matrix, Vector>;
    using MPISparseSpaceType = TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>;
    using MPIBaseSchemeType = Scheme<MPISparseSpaceType, LocalSpaceType>;

    // add schemes
    using MPIAlgebraicFluxCorrectedSteadyScalarSchemeType = AlgebraicFluxCorrectedSteadyScalarScheme<MPISparseSpaceType, LocalSpaceType>;
    py::class_<MPIAlgebraicFluxCorrectedSteadyScalarSchemeType, typename MPIAlgebraicFluxCorrectedSteadyScalarSchemeType::Pointer, MPIBaseSchemeType>(m, "MPIAlgebraicFluxCorrectedSteadyScalarScheme")
        .def(py::init<const double, const Flags&>())
        .def(py::init<const double, const Flags&, const Variable<int>&>());

    using MPISteadyScalarSchemeType = SteadyScalarScheme<MPISparseSpaceType, LocalSpaceType>;
    py::class_<MPISteadyScalarSchemeType, typename MPISteadyScalarSchemeType::Pointer, MPIBaseSchemeType>(m, "MPISteadyScalarScheme")
        .def(py::init<const double>());

    using MPIBossakRelaxationScalarSchemeType = BossakRelaxationScalarScheme<MPISparseSpaceType, LocalSpaceType>;
    py::class_<MPIBossakRelaxationScalarSchemeType, typename MPIBossakRelaxationScalarSchemeType::Pointer, MPIBaseSchemeType>(m, "MPIBossakRelaxationScalarScheme")
        .def(py::init<const double, const double, const Variable<double>&>());

}

} // namespace Python
} // namespace Kratos
