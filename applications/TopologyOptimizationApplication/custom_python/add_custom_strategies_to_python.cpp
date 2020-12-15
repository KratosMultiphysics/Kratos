// ==============================================================================
//  KratosTopologyOptimizationApplication
//
//  License:         BSD License
//                   license: TopologyOptimizationApplication/license.txt
//
//  Main authors:    Baumgärtner Daniel, https://github.com/dbaumgaertner
//                   Octaviano Malfavón Farías
//                   Eric Gonzales
//
// ==============================================================================

// External includes


// Project includes
#include "includes/define.h"
#include "containers/flags.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "spaces/ublas_space.h"

// Schemes
#include "custom_strategies/structure_adjoint_sensitivity_strategy.h"
#include "custom_strategies/custom_schemes/residualbased_incrementalupdate_static_simp_scheme.hpp"

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
	typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;
	typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;

	//custom scheme types
	typedef ResidualBasedIncrementalUpdateStaticSIMPScheme< SparseSpaceType, LocalSpaceType > ResidualBasedIncrementalUpdateStaticSIMPSchemeType;


	// =============================================================================================================================================
	// Scheme Classes
	// =============================================================================================================================================

	// Static TIMP Scheme Type
	py::class_< ResidualBasedIncrementalUpdateStaticSIMPSchemeType, ResidualBasedIncrementalUpdateStaticSIMPSchemeType::Pointer, BaseSchemeType >
	(m,"ResidualBasedIncrementalUpdateStaticSIMPScheme")
	.def(py::init<>());

	// =============================================================================================================================================
	// Strategy Classes
	// =============================================================================================================================================
///###############################################################################################auskommentiert
	/* py::class_< StructureAdjointSensitivityStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >, BaseSolvingStrategyType >
	(m, "StructureAdjointSensitivityStrategy")
	.def(py::init<ModelPart&, LinearSolverType::Pointer, int>())
	.def("ComputeStrainEnergySensitivities",&StructureAdjointSensitivityStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::ComputeStrainEnergySensitivities)
	.def("ComputeVolumeFractionSensitivities",&StructureAdjointSensitivityStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::ComputeVolumeFractionSensitivities)
	; */

	
}

}  // namespace Python.

} // Namespace Kratos
