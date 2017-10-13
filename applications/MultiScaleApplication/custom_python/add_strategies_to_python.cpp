/*
==============================================================================
KratosMultiScaleApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


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

==============================================================================
*/
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Massimo Petracca $
//   Date:                $Date: 2013-10-03 19:37:00 $
//   Revision:            $Revision: 1.00 $
//
//

// System includes

// External includes
#include <boost/python.hpp>
#include "spaces/ublas_space.h"

// Project includes

#include "linear_solvers/linear_solver.h"
#include "add_strategies_to_python.h"

#include "custom_strategies/strategies/static_general_strategy.h"
#include "custom_strategies/strategies/static_general_strategy_krylov_newton.h"
#include "custom_strategies/strategies/arc_length_strategy.h"
#include "custom_strategies/strategies/arc_length_riks_strategy.h"
#include "custom_strategies/strategies/arc_length_diss_nrg_strategy.h"
#include "custom_strategies/strategies/displacement_control_strategy.h"

#include "custom_strategies/builder_and_solvers/static_general_builder_and_solver.h"

#include "custom_strategies/schemes/static_general_scheme.h"
#include "custom_strategies/schemes/rve_static_scheme.h"

#include "custom_strategies/convergencecriterias/residual_norm_criteria.h"
#include "custom_strategies/convergencecriterias/displacement_norm_criteria.h"
#include "custom_strategies/convergencecriterias/energy_norm_criteria.h"


namespace Kratos
{
namespace Python
{

using namespace boost::python;



typedef UblasSpace<double, CompressedMatrix, Vector>	SparseSpaceType;
typedef UblasSpace<double, Matrix, Vector>				LocalSpaceType;



void MoveMesh(Scheme< SparseSpaceType, LocalSpaceType >& dummy, ModelPart::NodesContainerType& rNodes)
{
    for (ModelPart::NodeIterator i = rNodes.begin(); i != rNodes.end(); ++i)
    {
        const array_1d<double, 3 > & disp = i->FastGetSolutionStepValue(DISPLACEMENT);
        (i)->X() = (i)->X0() + disp[0];
        (i)->Y() = (i)->Y0() + disp[1];
        (i)->Z() = (i)->Z0() + disp[2];
    }
}



void AddStrategiesToPython()
{

	typedef LinearSolver<SparseSpaceType, LocalSpaceType >							LinearSolverType;
	typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >	BaseSolvingStrategyType;
	typedef Scheme< SparseSpaceType, LocalSpaceType >								BaseSchemeType;

	typedef ConvergenceCriteria< SparseSpaceType, LocalSpaceType >					TConvergenceCriteriaType;
	typedef ConvergenceCriteria< SparseSpaceType, LocalSpaceType > ::Pointer		TConvergenceCriteriaPointer;

	typedef BuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType >	BuilderAndSolverType;

	// =========================================================================================================
	//
	// Strategies
	//
	// =========================================================================================================

	class_< BaseSolvingStrategyType, boost::noncopyable > (
		"SolvingStrategyBase", 
		init < ModelPart&, bool >())
        .def("Predict", &BaseSolvingStrategyType::Predict)
        .def("Solve", &BaseSolvingStrategyType::Solve)
        .def("IsConverged", &BaseSolvingStrategyType::IsConverged)
        .def("CalculateOutputData", &BaseSolvingStrategyType::CalculateOutputData)
        .def("SetEchoLevel", &BaseSolvingStrategyType::SetEchoLevel)
        .def("GetEchoLevel", &BaseSolvingStrategyType::GetEchoLevel)
        .def("SetRebuildLevel", &BaseSolvingStrategyType::SetRebuildLevel)
        .def("GetRebuildLevel", &BaseSolvingStrategyType::GetRebuildLevel)
        .def("SetMoveMeshFlag", &BaseSolvingStrategyType::SetMoveMeshFlag)
        .def("MoveMeshFlag", &BaseSolvingStrategyType::MoveMeshFlag)
        .def("MoveMesh", &BaseSolvingStrategyType::MoveMesh)
        .def("Clear", &BaseSolvingStrategyType::Clear)
        .def("Check", &BaseSolvingStrategyType::Check)
        ;

	typedef StaticGeneralStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > StaticGeneralStrategyType;
	class_< StaticGeneralStrategyType, bases< BaseSolvingStrategyType >, boost::noncopyable >(
		"StaticGeneralStrategy",
		init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, TConvergenceCriteriaType::Pointer, int, bool, bool, bool>())
		.def(init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, TConvergenceCriteriaType::Pointer, BuilderAndSolverType::Pointer, int, bool, bool, bool >())
		.def("SetMaxIterationNumber", &StaticGeneralStrategyType::SetMaxIterationNumber)
		.def("GetMaxIterationNumber", &StaticGeneralStrategyType::GetMaxIterationNumber)
		.def("SetKeepSystemConstantDuringIterations", &StaticGeneralStrategyType::SetKeepSystemConstantDuringIterations)
		.def("GetKeepSystemConstantDuringIterations", &StaticGeneralStrategyType::GetKeepSystemConstantDuringIterations)
		;
	
	typedef StaticGeneralStrategyKrylovNewton< SparseSpaceType, LocalSpaceType, LinearSolverType > StaticGeneralStrategyKrylovNewtonType;
	class_< StaticGeneralStrategyKrylovNewtonType, bases< BaseSolvingStrategyType >, boost::noncopyable >(
		"StaticGeneralStrategyKrylovNewton",
		init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, TConvergenceCriteriaType::Pointer, int, bool, bool, bool>())
		.def(init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, TConvergenceCriteriaType::Pointer, BuilderAndSolverType::Pointer, int, bool, bool, bool >())
		.def("SetMaxIterationNumber", &StaticGeneralStrategyKrylovNewtonType::SetMaxIterationNumber)
		.def("GetMaxIterationNumber", &StaticGeneralStrategyKrylovNewtonType::GetMaxIterationNumber)
		.def("SetKeepSystemConstantDuringIterations", &StaticGeneralStrategyKrylovNewtonType::SetKeepSystemConstantDuringIterations)
		.def("GetKeepSystemConstantDuringIterations", &StaticGeneralStrategyKrylovNewtonType::GetKeepSystemConstantDuringIterations)
		;

	typedef ArcLengthStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > ArcLengthStrategyType;
	class_< ArcLengthStrategyType, bases< BaseSolvingStrategyType >, boost::noncopyable >(
		"ArcLengthStrategy",
		init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, TConvergenceCriteriaType::Pointer, double, double, double, double, int, bool, bool, bool>())
		.def("SetMaxIterationNumber", &ArcLengthStrategyType::SetMaxIterationNumber)
		.def("GetMaxIterationNumber", &ArcLengthStrategyType::GetMaxIterationNumber)
		.def("SetKeepSystemConstantDuringIterations", &ArcLengthStrategyType::SetKeepSystemConstantDuringIterations)
		.def("GetKeepSystemConstantDuringIterations", &ArcLengthStrategyType::GetKeepSystemConstantDuringIterations)
		.def("PrepareFirstSolve", &ArcLengthStrategyType::PrepareFirstSolve)
		;

	typedef ArcLengthRiksStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > ArcLengthRiksStrategyType;
	class_< ArcLengthRiksStrategyType, bases< BaseSolvingStrategyType >, boost::noncopyable >(
		"ArcLengthRiksStrategy",
		init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, TConvergenceCriteriaType::Pointer, int, bool, bool, bool>())
		.def(init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, TConvergenceCriteriaType::Pointer, BuilderAndSolverType::Pointer, int, bool, bool, bool >())
		.def("SetMaxIterationNumber", &ArcLengthRiksStrategyType::SetMaxIterationNumber)
		.def("GetMaxIterationNumber", &ArcLengthRiksStrategyType::GetMaxIterationNumber)
		.def("SetKeepSystemConstantDuringIterations", &ArcLengthRiksStrategyType::SetKeepSystemConstantDuringIterations)
		.def("GetKeepSystemConstantDuringIterations", &ArcLengthRiksStrategyType::GetKeepSystemConstantDuringIterations)
		.def("PrepareFirstSolve", &ArcLengthRiksStrategyType::PrepareFirstSolve)
		.def("SetLoadFactors", &ArcLengthRiksStrategyType::SetLoadFactors)
		;
	
	typedef ArcLengthDissNrgStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > ArcLengthDissNrgStrategyType;
	class_< ArcLengthDissNrgStrategyType, bases< BaseSolvingStrategyType >, boost::noncopyable >(
		"ArcLengthDissNrgStrategy",
		init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, TConvergenceCriteriaType::Pointer, int, bool, bool, bool>())
		.def(init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, TConvergenceCriteriaType::Pointer, BuilderAndSolverType::Pointer, int, bool, bool, bool >())
		.def("SetMaxIterationNumber", &ArcLengthDissNrgStrategyType::SetMaxIterationNumber)
		.def("GetMaxIterationNumber", &ArcLengthDissNrgStrategyType::GetMaxIterationNumber)
		.def("SetKeepSystemConstantDuringIterations", &ArcLengthDissNrgStrategyType::SetKeepSystemConstantDuringIterations)
		.def("GetKeepSystemConstantDuringIterations", &ArcLengthDissNrgStrategyType::GetKeepSystemConstantDuringIterations)
		.def("PrepareFirstSolve", &ArcLengthDissNrgStrategyType::PrepareFirstSolve)
		.def("SetLoadFactors", &ArcLengthDissNrgStrategyType::SetLoadFactors)
		;

	typedef DisplacementControlStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > DisplacementControlStrategyType;
	class_< DisplacementControlStrategyType, bases< BaseSolvingStrategyType >, boost::noncopyable >(
		"DisplacementControlStrategy",
		init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, TConvergenceCriteriaType::Pointer, int, bool, bool, bool, unsigned int, unsigned int, double>())
		.def(init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, TConvergenceCriteriaType::Pointer, BuilderAndSolverType::Pointer, int, bool, bool, bool, unsigned int, unsigned int, double >())
		.def("SetMaxIterationNumber", &DisplacementControlStrategyType::SetMaxIterationNumber)
		.def("GetMaxIterationNumber", &DisplacementControlStrategyType::GetMaxIterationNumber)
		.def("SetKeepSystemConstantDuringIterations", &DisplacementControlStrategyType::SetKeepSystemConstantDuringIterations)
		.def("GetKeepSystemConstantDuringIterations", &DisplacementControlStrategyType::GetKeepSystemConstantDuringIterations)
		.def("PrepareFirstSolve", &DisplacementControlStrategyType::PrepareFirstSolve)
		.def("SetLoadFactors", &DisplacementControlStrategyType::SetLoadFactors)
		;

	// =========================================================================================================
	//
	// BuilderAndSolvers
	//
	// =========================================================================================================

	typedef StaticGeneralBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > StaticGeneralBuilderAndSolverType;
	class_<StaticGeneralBuilderAndSolverType, bases< BuilderAndSolverType >, boost::noncopyable >(
		"StaticGeneralBuilderAndSolver",
		init<LinearSolverType::Pointer>())
		;

	// =========================================================================================================
	//
	// Schemes
	//
	// =========================================================================================================

	class_< BaseSchemeType, boost::noncopyable >
            ("Scheme", init< >())
            .def("Initialize", &BaseSchemeType::Initialize)
            .def("SchemeIsInitialized", &BaseSchemeType::SchemeIsInitialized)
            .def("ElementsAreInitialized", &BaseSchemeType::ElementsAreInitialized)
            .def("ConditionsAreInitialized", &BaseSchemeType::ConditionsAreInitialized)
            .def("InitializeElements", &BaseSchemeType::InitializeElements)
            .def("InitializeConditions", &BaseSchemeType::InitializeConditions)
            .def("InitializeSolutionStep", &BaseSchemeType::InitializeSolutionStep)
            .def("FinalizeSolutionStep", &BaseSchemeType::FinalizeSolutionStep)
            .def("InitializeNonLinIteration", &BaseSchemeType::InitializeNonLinIteration)
            .def("FinalizeNonLinIteration", &BaseSchemeType::FinalizeNonLinIteration)
            .def("Predict", &BaseSchemeType::Predict)
            .def("Update", &BaseSchemeType::Update)
            .def("CalculateOutputData", &BaseSchemeType::CalculateOutputData)
            .def("Clean", &BaseSchemeType::Clean)
            .def("Clear",&BaseSchemeType::Clear)
            .def("MoveMesh", MoveMesh)
            .def("Check", &BaseSchemeType::Check)
            ;

	typedef StaticGeneralScheme< SparseSpaceType, LocalSpaceType> StaticGeneralSchemeType;
    class_< StaticGeneralSchemeType, bases< BaseSchemeType >, boost::noncopyable >(
            "StaticGeneralScheme", 
			init< >())
			;

	typedef RveStaticScheme< SparseSpaceType, LocalSpaceType> RveStaticSchemeType;
    class_< RveStaticSchemeType, bases< BaseSchemeType >, boost::noncopyable >(
            "RveStaticScheme", 
			init< >())
			;

	// =========================================================================================================
	//
	// Convergence criteria
	//
	// =========================================================================================================

	typedef ConvergenceCriteria< SparseSpaceType, LocalSpaceType > ConvergenceCriteriaBaseType;
	class_< ConvergenceCriteriaBaseType, boost::noncopyable > (
		"ConvergenceCriteriaBase", 
		init<>())
		.def("SetActualizeRHSFlag", &ConvergenceCriteria<SparseSpaceType, LocalSpaceType >::SetActualizeRHSFlag)
		.def("GetActualizeRHSflag", &ConvergenceCriteria<SparseSpaceType, LocalSpaceType >::GetActualizeRHSflag)
		.def("PreCriteria", &ConvergenceCriteria<SparseSpaceType, LocalSpaceType >::PreCriteria)
		.def("PostCriteria", &ConvergenceCriteria<SparseSpaceType, LocalSpaceType >::PostCriteria)
		.def("Initialize", &ConvergenceCriteria<SparseSpaceType, LocalSpaceType >::Initialize)
		.def("InitializeSolutionStep", &ConvergenceCriteria<SparseSpaceType, LocalSpaceType >::InitializeSolutionStep)
		.def("FinalizeSolutionStep", &ConvergenceCriteria<SparseSpaceType, LocalSpaceType >::FinalizeSolutionStep)
		.def("Check", &ConvergenceCriteria<SparseSpaceType, LocalSpaceType >::Check)
		;

	typedef ResidualNormCriteria<SparseSpaceType, LocalSpaceType > ResidualNormCriteriaType;
    class_< ResidualNormCriteriaType, bases<ConvergenceCriteriaBaseType>, boost::noncopyable >(
		"ResidualNormCriteria", 
		init< double, double>())
		.def(init< double, double, bool >())
		;

	typedef DisplacementNormCriteria<SparseSpaceType, LocalSpaceType > DisplacementNormCriteriaType;
    class_< DisplacementNormCriteriaType, bases<ConvergenceCriteriaBaseType>, boost::noncopyable >(
		"DisplacementNormCriteria", 
		init< double, double>())
		.def(init< double, double, bool >())
		;

	typedef EnergyNormCriteria<SparseSpaceType, LocalSpaceType > EnergyNormCriteriaType;
    class_< EnergyNormCriteriaType, bases<ConvergenceCriteriaBaseType>, boost::noncopyable >(
		"EnergyNormCriteria", 
		init< double, double>())
		.def(init< double, double, bool >())
		;

}


}

}