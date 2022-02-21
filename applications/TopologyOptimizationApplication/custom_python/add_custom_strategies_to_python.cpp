//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Baumgärtner Daniel, https://github.com/dbaumgaertner
//                   Octaviano Malfavón Farías
//                   Eric Gonzales
//					 Philipp Hofer
//					 Erich Wehrle
//
// ==============================================================================

// External includes

// Project includes
#include "custom_python/add_custom_strategies_to_python.h"
#include "spaces/ublas_space.h"

// Schemes
#include "custom_strategies/structure_adjoint_sensitivity_strategy.h"
#include "custom_strategies/custom_schemes/residualbased_incrementalupdate_static_simp_scheme.h"

// Linear solvers
#include "linear_solvers/linear_solver.h"


namespace Kratos
{

namespace Python
{


void  AddCustomStrategiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

    //base types
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
    typedef SolvingStrategy< SparseSpaceType, LocalSpaceType> BaseSolvingStrategyType;
    typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;

    //custom scheme types
    typedef ResidualBasedIncrementalUpdateStaticSIMPScheme< SparseSpaceType, LocalSpaceType > ResidualBasedIncrementalUpdateStaticSIMPSchemeType;
    typedef StructureAdjointSensitivityStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > StructureAdjointSensitivityStrategyType;
    // =============================================================================================================================================
    // Scheme Classes
    // =============================================================================================================================================

    // Static TIMP Scheme Type
    py::class_< ResidualBasedIncrementalUpdateStaticSIMPSchemeType, typename ResidualBasedIncrementalUpdateStaticSIMPSchemeType::Pointer, BaseSchemeType >
    (m,"ResidualBasedIncrementalUpdateStaticSIMPScheme")
    .def(py::init<>())
    .def("Initialize", &ResidualBasedIncrementalUpdateStaticSIMPScheme<SparseSpaceType, LocalSpaceType>::Initialize)
    ;

    // =============================================================================================================================================
    // Strategy Classes
    // =============================================================================================================================================

    py::class_< StructureAdjointSensitivityStrategyType, typename StructureAdjointSensitivityStrategyType::Pointer, BaseSolvingStrategyType >
    (m, "StructureAdjointSensitivityStrategy")
    .def(py::init<ModelPart&, LinearSolverType::Pointer, int>())
    .def("ComputeStrainEnergySensitivities",&StructureAdjointSensitivityStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::ComputeStrainEnergySensitivities)
    .def("ComputeVolumeFractionSensitivities",&StructureAdjointSensitivityStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::ComputeVolumeFractionSensitivities)
    ; 

	
}

}  // namespace Python.

} // Namespace Kratos
