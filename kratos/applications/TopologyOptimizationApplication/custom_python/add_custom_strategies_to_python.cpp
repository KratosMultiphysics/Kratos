// ==============================================================================
/*
 KratosTopologyOptimizationApplication
 A library based on:
 Kratos
 A General Purpose Software for Multi-Physics Finite Element Analysis
 (Released on march 05, 2007).

 Copyright (c) 2016: Daniel Baumgaertner
                     daniel.baumgaertner@tum.de
                     Chair of Structural Analysis
                     Technische Universitaet Muenchen
                     Arcisstrasse 21 80333 Munich, Germany

 Permission is hereby granted, free  of charge, to any person obtaining
 a  copy  of this  software  and  associated  documentation files  (the
 "Software"), to  deal in  the Software without  restriction, including
 without limitation  the rights to  use, copy, modify,  merge, publish,
 distribute,  sublicense and/or  sell copies  of the  Software,  and to
 permit persons to whom the Software  is furnished to do so, subject to
 the following condition:

 Distribution of this code for  any  commercial purpose  is permissible
 ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

 The  above  copyright  notice  and  this permission  notice  shall  be
 included in all copies or substantial portions of the Software.

 THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
 EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
 CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
 TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
 SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
//==============================================================================
//
//   Project Name:        KratosTopology                        $
//   Last modified by:	  $Author:   daniel.baumgaertner@tum.de $
// 						  $Co-Author: Octaviano Malfavón Farías $
//   Date:                $Date:                    August 2016 $
//   Revision:            $Revision:                        0.0 $
//
// ==============================================================================

// External includes
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/timer.hpp>

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
using namespace boost::python;

void  AddCustomStrategiesToPython()
{
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
	class_< ResidualBasedIncrementalUpdateStaticSIMPSchemeType,bases< BaseSchemeType >, boost::noncopyable >
	("ResidualBasedIncrementalUpdateStaticSIMPScheme", init< >())
	.def("Initialize", &ResidualBasedIncrementalUpdateStaticSIMPScheme<SparseSpaceType, LocalSpaceType>::Initialize)
	;

	// =============================================================================================================================================
	// Strategy Classes
	// =============================================================================================================================================

	class_< StructureAdjointSensitivityStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >,bases< BaseSolvingStrategyType >,  boost::noncopyable >
	("StructureAdjointSensitivityStrategy",init<ModelPart&, LinearSolverType::Pointer, int>())
	.def("ComputeStrainEnergySensitivities",&StructureAdjointSensitivityStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::ComputeStrainEnergySensitivities)
	.def("ComputeVolumeFractionSensitivities",&StructureAdjointSensitivityStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::ComputeVolumeFractionSensitivities)
	;
}

}  // namespace Python.

} // Namespace Kratos
