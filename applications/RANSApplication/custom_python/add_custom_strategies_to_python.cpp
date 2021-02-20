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

// System includes

// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define_python.h"
#include "spaces/ublas_space.h"

// schemes
#include "custom_strategies/bossak_relaxation_scalar_scheme.h"
#include "custom_strategies/steady_scalar_scheme.h"
#include "custom_strategies/algebraic_flux_corrected_steady_scalar_scheme.h"

// adjoint schemes
#include "custom_strategies/schemes/simple_steady_adjoint_scheme.h"
#include "custom_strategies/schemes/velocity_bossak_adjoint_scheme.h"

// sensitivity builder schemes
#include "custom_strategies/schemes/simple_steady_sensitivity_builder_scheme.h"
#include "custom_strategies/schemes/velocity_bossak_sensitivity_builder_scheme.h"

// Include base h
#include "custom_python/add_custom_strategies_to_python.h"

namespace Kratos
{
namespace Python
{
void AddCustomStrategiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    using LocalSpaceType = UblasSpace<double, Matrix, Vector>;
    using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
    using BaseSchemeType = Scheme<SparseSpaceType, LocalSpaceType>;

    // add schemes
    using SteadyScalarSchemeType = SteadyScalarScheme<SparseSpaceType, LocalSpaceType>;
    py::class_<SteadyScalarSchemeType, typename SteadyScalarSchemeType::Pointer, BaseSchemeType>(m, "SteadyScalarScheme")
        .def(py::init<const double>());

    using AlgebraicFluxCorrectedSteadyScalarSchemeType = AlgebraicFluxCorrectedSteadyScalarScheme<SparseSpaceType, LocalSpaceType>;
    py::class_<AlgebraicFluxCorrectedSteadyScalarSchemeType, typename AlgebraicFluxCorrectedSteadyScalarSchemeType::Pointer, BaseSchemeType>(m, "AlgebraicFluxCorrectedSteadyScalarScheme")
        .def(py::init<const double, const Flags&>())
        .def(py::init<const double, const Flags&, const Variable<int>&>());

    using BossakRelaxationScalarSchemeType = BossakRelaxationScalarScheme<SparseSpaceType, LocalSpaceType>;
    py::class_<BossakRelaxationScalarSchemeType, typename BossakRelaxationScalarSchemeType::Pointer, BaseSchemeType>(m, "BossakRelaxationScalarScheme")
        .def(py::init<const double, const double, const Variable<double>&>());

    using  SimpleSteadyAdjointScheme2DType = SimpleSteadyAdjointScheme<2, SparseSpaceType, LocalSpaceType, 5>;
    py::class_<SimpleSteadyAdjointScheme2DType, typename SimpleSteadyAdjointScheme2DType::Pointer, BaseSchemeType>
        (m, "RansSimpleSteadyAdjointScheme2D")
        .def(py::init<AdjointResponseFunction::Pointer>())
        ;

    using  SimpleSteadyAdjointScheme3DType = SimpleSteadyAdjointScheme<3, SparseSpaceType, LocalSpaceType, 6>;
    py::class_<SimpleSteadyAdjointScheme3DType, typename SimpleSteadyAdjointScheme3DType::Pointer, BaseSchemeType>
        (m, "RansSimpleSteadyAdjointScheme3D")
        .def(py::init<AdjointResponseFunction::Pointer>())
        ;

    using  VelocityBossakAdjointScheme2DType = VelocityBossakAdjointScheme<2, SparseSpaceType, LocalSpaceType, 5>;
    py::class_<VelocityBossakAdjointScheme2DType, typename VelocityBossakAdjointScheme2DType::Pointer, BaseSchemeType>
        (m, "RansVelocityBossakAdjointScheme2D")
        .def(py::init<Parameters, AdjointResponseFunction::Pointer>())
        ;

    using  VelocityBossakAdjointScheme3DType = VelocityBossakAdjointScheme<3, SparseSpaceType, LocalSpaceType, 6>;
    py::class_<VelocityBossakAdjointScheme3DType, typename VelocityBossakAdjointScheme3DType::Pointer, BaseSchemeType>
        (m, "RansVelocityBossakAdjointScheme3D")
        .def(py::init<Parameters, AdjointResponseFunction::Pointer>())
        ;

    using SimpleSteadySensitivityBuilderScheme2DType = SimpleSteadySensitivityBuilderScheme<2, 5>;
    py::class_<SimpleSteadySensitivityBuilderScheme2DType, typename SimpleSteadySensitivityBuilderScheme2DType::Pointer, SensitivityBuilderScheme>
        (m, "RansSimpleSteadySensitivityBuilderScheme2D")
        .def(py::init())
        ;

    using SimpleSteadySensitivityBuilderScheme3DType = SimpleSteadySensitivityBuilderScheme<3, 6>;
    py::class_<SimpleSteadySensitivityBuilderScheme3DType, typename SimpleSteadySensitivityBuilderScheme3DType::Pointer, SensitivityBuilderScheme>
        (m, "RansSimpleSteadySensitivityBuilderScheme3D")
        .def(py::init())
        ;

    using VelocityBossakSensitivityBuilderScheme2DType = VelocityBossakSensitivityBuilderScheme<2, 5>;
    py::class_<VelocityBossakSensitivityBuilderScheme2DType, typename VelocityBossakSensitivityBuilderScheme2DType::Pointer, SensitivityBuilderScheme>
        (m, "RansVelocityBossakSensitivityBuilderScheme2D")
        .def(py::init<const double>())
        ;

    using VelocityBossakSensitivityBuilderScheme3DType = VelocityBossakSensitivityBuilderScheme<3, 6>;
    py::class_<VelocityBossakSensitivityBuilderScheme3DType, typename VelocityBossakSensitivityBuilderScheme3DType::Pointer, SensitivityBuilderScheme>
        (m, "RansVelocityBossakSensitivityBuilderScheme3D")
        .def(py::init<const double>())
        ;

}

} // namespace Python.
} // Namespace Kratos
