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

// adjoint schemes
#include "custom_strategies/schemes/simple_steady_adjoint_scheme.h"
#include "custom_strategies/schemes/velocity_bossak_adjoint_scheme.h"

// sensitivity builder schemes
#include "custom_strategies/schemes/simple_steady_sensitivity_builder_scheme.h"
#include "custom_strategies/schemes/velocity_bossak_sensitivity_builder_scheme.h"

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

    using  SimpleSteadyAdjointScheme2DType = SimpleSteadyAdjointScheme<2, MPISparseSpaceType, LocalSpaceType, 5>;
    py::class_<SimpleSteadyAdjointScheme2DType, typename SimpleSteadyAdjointScheme2DType::Pointer, MPIBaseSchemeType>
        (m, "MPIRansSimpleSteadyAdjointScheme2D")
        .def(py::init<AdjointResponseFunction::Pointer>())
        ;

    using  SimpleSteadyAdjointScheme3DType = SimpleSteadyAdjointScheme<3, MPISparseSpaceType, LocalSpaceType, 6>;
    py::class_<SimpleSteadyAdjointScheme3DType, typename SimpleSteadyAdjointScheme3DType::Pointer, MPIBaseSchemeType>
        (m, "MPIRansSimpleSteadyAdjointScheme3D")
        .def(py::init<AdjointResponseFunction::Pointer>())
        ;

    using  VelocityBossakAdjointScheme2DType = VelocityBossakAdjointScheme<2, MPISparseSpaceType, LocalSpaceType, 5>;
    py::class_<VelocityBossakAdjointScheme2DType, typename VelocityBossakAdjointScheme2DType::Pointer, MPIBaseSchemeType>
        (m, "MPIRansVelocityBossakAdjointScheme2D")
        .def(py::init<Parameters, AdjointResponseFunction::Pointer>())
        ;

    using  VelocityBossakAdjointScheme3DType = VelocityBossakAdjointScheme<3, MPISparseSpaceType, LocalSpaceType, 6>;
    py::class_<VelocityBossakAdjointScheme3DType, typename VelocityBossakAdjointScheme3DType::Pointer, MPIBaseSchemeType>
        (m, "MPIRansVelocityBossakAdjointScheme3D")
        .def(py::init<Parameters, AdjointResponseFunction::Pointer>())
        ;
}

} // namespace Python
} // namespace Kratos
