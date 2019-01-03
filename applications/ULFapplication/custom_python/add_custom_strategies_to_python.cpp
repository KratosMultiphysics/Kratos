/*
==============================================================================
KratosULFApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).
Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
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
//   Last modified by:    $Author: anonymous $
//   Date:                $Date: 2008-04-07 09:50:09 $
//   Revision:            $Revision: 1.3 $
//
//


// System includes


// External includes


// Project includes
#include <pybind11/pybind11.h>
#include "includes/define.h"
#include "includes/define_python.h"

#include "custom_python/add_custom_strategies_to_python.h"
#include "processes/process.h"
#include "spaces/ublas_space.h"
#include <boost/timer.hpp>


//builder_and_solvers
#include "custom_strategies/builder_and_solvers/residualbased_elimination_quasiincompresible_builder_and_solver.h"
#include "custom_strategies/strategies/modified_linear_strategy.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "custom_strategies/strategies/runge_kutta_fracstep_GLS_strategy.h"
//schemes
#include "solving_strategies/schemes/scheme.h"
#include "custom_strategies/schemes/residualbased_predictorcorrector_bossak_scheme.h"

//linear solvers
#include "linear_solvers/linear_solver.h"

namespace Kratos
{

namespace Python
{
namespace py = pybind11;

void  AddCustomStrategiesToPython(pybind11::module& m)
{
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;


    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
    typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;

    typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;

    typedef BuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType> BuilderAndSolverType;
    typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;
//			typedef ResidualBasedEliminationBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType> ResidualBasedEliminationBuilderAndSolverType;

    //typedef ResidualBasedPredictorCorrectorBossakScheme< SparseSpaceType, LocalSpaceType >   ResidualBasedPredictorCorrectorBossakSchemeType;
    //********************************************************************
    //********************************************************************
    //typedef ResidualBasedEliminationQuasiIncompressibleBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType, 2> //ResidualBasedEliminationQuasiIncompressibleBuilderAndSolverType2D;

    typedef ResidualBasedEliminationQuasiIncompressibleBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType, 2> ResidualBasedIncompressibleBuilderType2D;



    //py::class_< ResidualBasedEliminationQuasiIncompressibleBuilderAndSolverType2D, boost::noncopyable>
//("ResidualBasedEliminationQuasiIncompressibleBuilderAndSolver2D", init< LinearSolverType::Pointer>() )
//
//    //********************************************************************
    //********************************************************************


py::class_< ResidualBasedPredictorCorrectorBossakScheme< SparseSpaceType, LocalSpaceType>,
                    typename ResidualBasedPredictorCorrectorBossakScheme< SparseSpaceType, LocalSpaceType>::Pointer,
                    BaseSchemeType >
                    (m, "ResidualBasedPredictorCorrectorBossakScheme")
                    .def(py::init< double >()
                    );


    py::class_< ResidualBasedIncompressibleBuilderType2D, ResidualBasedIncompressibleBuilderType2D::Pointer, BuilderAndSolverType > (m, "ResidualBasedIncompressibleBuilder2D")
    .def(py::init< LinearSolverType::Pointer>() )
    .def("AssembleLHS", &ResidualBasedIncompressibleBuilderType2D::AssembleLHS )
    .def("AssembleRHS", &ResidualBasedIncompressibleBuilderType2D::AssembleRHS )
    .def("BuildAndSolve", &ResidualBasedIncompressibleBuilderType2D::BuildAndSolve)
    .def("SetUpDofSet", &ResidualBasedIncompressibleBuilderType2D::SetUpDofSet)
    .def("SetUpSystem", &ResidualBasedIncompressibleBuilderType2D::SetUpSystem)
    .def("ResizeAndInitializeVectors", &ResidualBasedIncompressibleBuilderType2D::ResizeAndInitializeVectors)
    .def("Build", &ResidualBasedIncompressibleBuilderType2D::Build)
    .def("Solve", &ResidualBasedIncompressibleBuilderType2D::SystemSolve)
    .def("ConstructMatrixStructure", &ResidualBasedIncompressibleBuilderType2D::ConstructMatrixStructure)
    .def("ConstructMatrixStructure_Mconsistent", &ResidualBasedIncompressibleBuilderType2D::ConstructMatrixStructure_Mconsistent)
    .def("ConstructMatrixStructure_DivergenceMatrixD", &ResidualBasedIncompressibleBuilderType2D::ConstructMatrixStructure_DivergenceMatrixD)
    .def("BuildAuxiliaries", &ResidualBasedIncompressibleBuilderType2D::BuildAuxiliaries)
    .def("AssembleMassMatrices", &ResidualBasedIncompressibleBuilderType2D::AssembleMassMatrices)
    .def("calc_GMinvD_prod", &ResidualBasedIncompressibleBuilderType2D::calc_GMinvD_prod)
    .def("CalculatePreconditionerDiagonalMatrix",  &ResidualBasedIncompressibleBuilderType2D::CalculatePreconditionerDiagonalMatrix)
    .def("calc_prod_precond_vec", &ResidualBasedIncompressibleBuilderType2D::calc_prod_precond_vec)
    .def("ModifyForDirichlet", &ResidualBasedIncompressibleBuilderType2D::ModifyForDirichlet)
    .def("UpdatePressures", &ResidualBasedIncompressibleBuilderType2D::UpdatePressures)
    .def("ReturnDx", &ResidualBasedIncompressibleBuilderType2D::ReturnDx)
    .def("UpdatePressuresNew", &ResidualBasedIncompressibleBuilderType2D::UpdatePressuresNew)
    .def("UpdateAfterProjection", &ResidualBasedIncompressibleBuilderType2D::UpdateAfterProjection)
    .def("ComputePressureAtFreeSurface", &ResidualBasedIncompressibleBuilderType2D::ComputePressureAtFreeSurface)
    .def("SavePressureIteration", &ResidualBasedIncompressibleBuilderType2D::SavePressureIteration)
    .def("FractionalStepProjection", &ResidualBasedIncompressibleBuilderType2D::FractionalStepProjection)
    .def("CalculateNodalPressureForce", &ResidualBasedIncompressibleBuilderType2D::CalculateNodalPressureForce )
    .def("ConvergenceCheck", &ResidualBasedIncompressibleBuilderType2D::ConvergenceCheck)
    //.def("BuildAuxiliariesFSI", &ResidualBasedIncompressibleBuilderType2D::BuildAuxiliariesFSI)
    //.def("ConstructMatrixStructure_Fluid_DivergenceMatrixD",&ResidualBasedIncompressibleBuilderType2D::ConstructMatrixStructure_Fluid_DivergenceMatrixD)

    ;

    typedef ResidualBasedEliminationQuasiIncompressibleBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType, 3> ResidualBasedEliminationQuasiIncompressibleBuilderAndSolverType3D;


    py::class_< ResidualBasedEliminationQuasiIncompressibleBuilderAndSolverType3D, ResidualBasedEliminationQuasiIncompressibleBuilderAndSolverType3D::Pointer, BuilderAndSolverType > (m,"ResidualBasedEliminationQuasiIncompressibleBuilderAndSolver3D")
    .def(py::init< LinearSolverType::Pointer>() )
    .def("AssembleLHS", &ResidualBasedEliminationQuasiIncompressibleBuilderAndSolverType3D::AssembleLHS )
    .def("AssembleRHS", &ResidualBasedEliminationQuasiIncompressibleBuilderAndSolverType3D::AssembleRHS )
    .def("BuildAndSolve", &ResidualBasedEliminationQuasiIncompressibleBuilderAndSolverType3D::BuildAndSolve)
    .def("SetUpDofSet", &ResidualBasedEliminationQuasiIncompressibleBuilderAndSolverType3D::SetUpDofSet)
    .def("SetUpSystem", &ResidualBasedEliminationQuasiIncompressibleBuilderAndSolverType3D::SetUpSystem)
    .def("ResizeAndInitializeVectors", &ResidualBasedEliminationQuasiIncompressibleBuilderAndSolverType3D::ResizeAndInitializeVectors)
    .def("Build", &ResidualBasedEliminationQuasiIncompressibleBuilderAndSolverType3D::Build)
    .def("Solve", &ResidualBasedEliminationQuasiIncompressibleBuilderAndSolverType3D::SystemSolve)
    .def("ConstructMatrixStructure", &ResidualBasedEliminationQuasiIncompressibleBuilderAndSolverType3D::ConstructMatrixStructure)
    .def("ConstructMatrixStructure_Mconsistent", &ResidualBasedEliminationQuasiIncompressibleBuilderAndSolverType3D::ConstructMatrixStructure_Mconsistent)
    .def("ConstructMatrixStructure_DivergenceMatrixD", &ResidualBasedEliminationQuasiIncompressibleBuilderAndSolverType3D::ConstructMatrixStructure_DivergenceMatrixD)
    .def("BuildAuxiliaries", &ResidualBasedEliminationQuasiIncompressibleBuilderAndSolverType3D::BuildAuxiliaries)
    .def("AssembleMassMatrices", &ResidualBasedEliminationQuasiIncompressibleBuilderAndSolverType3D::AssembleMassMatrices)
    .def("calc_GMinvD_prod", &ResidualBasedEliminationQuasiIncompressibleBuilderAndSolverType3D::calc_GMinvD_prod)
    .def("CalculatePreconditionerDiagonalMatrix",  &ResidualBasedEliminationQuasiIncompressibleBuilderAndSolverType3D::CalculatePreconditionerDiagonalMatrix)
    .def("calc_prod_precond_vec", &ResidualBasedEliminationQuasiIncompressibleBuilderAndSolverType3D::calc_prod_precond_vec)
    .def("ModifyForDirichlet", &ResidualBasedEliminationQuasiIncompressibleBuilderAndSolverType3D::ModifyForDirichlet)
    .def("UpdatePressures", &ResidualBasedEliminationQuasiIncompressibleBuilderAndSolverType3D::UpdatePressures)
    .def("ReturnDx", &ResidualBasedEliminationQuasiIncompressibleBuilderAndSolverType3D::ReturnDx)
    .def("UpdatePressuresNew", &ResidualBasedEliminationQuasiIncompressibleBuilderAndSolverType3D::UpdatePressuresNew)
    .def("CalculateNodalPressureForce", &ResidualBasedEliminationQuasiIncompressibleBuilderAndSolverType3D::CalculateNodalPressureForce )
    .def("ConvergenceCheck", &ResidualBasedEliminationQuasiIncompressibleBuilderAndSolverType3D::ConvergenceCheck)
    .def("UpdateAfterProjection", &ResidualBasedEliminationQuasiIncompressibleBuilderAndSolverType3D::UpdateAfterProjection)
    .def("ComputePressureAtFreeSurface", &ResidualBasedEliminationQuasiIncompressibleBuilderAndSolverType3D::ComputePressureAtFreeSurface)
    .def("SavePressureIteration", &ResidualBasedEliminationQuasiIncompressibleBuilderAndSolverType3D::SavePressureIteration)
    ;


    py::class_< RungeKuttaFracStepStrategy < 2, SparseSpaceType, LocalSpaceType, LinearSolverType >, RungeKuttaFracStepStrategy < 2, SparseSpaceType, LocalSpaceType, LinearSolverType >::Pointer, BaseSolvingStrategyType > (m,"RungeKuttaFracStepStrategy2D")
     .def(py::init < ModelPart&, LinearSolverType::Pointer, bool, bool, bool >())
     //.def("SolveStep1", &RungeKuttaFracStepStrategy < 2, SparseSpaceType, LocalSpaceType, LinearSolverType >::SolveStep1)
     .def("SolveStep2", &RungeKuttaFracStepStrategy < 2, SparseSpaceType, LocalSpaceType, LinearSolverType >::SolveStep2)
     .def("SolveStep3", &RungeKuttaFracStepStrategy < 2, SparseSpaceType, LocalSpaceType, LinearSolverType >::SolveStep3)
     //.def("SolveStep_ForwardEuler", &RungeKuttaFracStepStrategy < 2, SparseSpaceType, LocalSpaceType, LinearSolverType >::SolveStep_ForwardEuler)
     .def("SolveStep1", &RungeKuttaFracStepStrategy < 2, SparseSpaceType, LocalSpaceType, LinearSolverType >::SolveStep1)
     .def("Clear", &RungeKuttaFracStepStrategy < 2, SparseSpaceType, LocalSpaceType, LinearSolverType >::Clear)
     ;


//    py::class_< ResidualBasedSemiEulerianConvectionDiffusionStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >,
//            ResidualBasedSemiEulerianConvectionDiffusionStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::Pointer,
//            BaseSolvingStrategyType >
//            (m,"ResidualBasedSemiEulerianConvectionDiffusionStrategy")
//            .def( init<	ModelPart&, LinearSolverType::Pointer,	bool, int	>() )
//            .def("Clear",&ResidualBasedSemiEulerianConvectionDiffusionStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::Clear)
//            ;



    //********************************************************************
    //********************************************************************
    typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;
//strategy base class
    py::class_< BaseSolvingStrategyType, BaseSolvingStrategyType::Pointer >(m,"SolvingStrategy", py::module_local())
    .def(py::init< ModelPart&, bool >() )
    .def("Predict", &BaseSolvingStrategyType::Predict )
    .def("Solve", &BaseSolvingStrategyType::Solve )
    .def("IsConverged", &BaseSolvingStrategyType::IsConverged )
    .def("CalculateOutputData", &BaseSolvingStrategyType::CalculateOutputData )
    .def("SetEchoLevel", &BaseSolvingStrategyType::SetEchoLevel )
    .def("GetEchoLevel", &BaseSolvingStrategyType::GetEchoLevel )
    .def("SetRebuildLevel", &BaseSolvingStrategyType::SetRebuildLevel )
    .def("GetRebuildLevel", &BaseSolvingStrategyType::GetRebuildLevel )
    .def("SetMoveMeshFlag", &BaseSolvingStrategyType::SetMoveMeshFlag )
    .def("MoveMeshFlag", &BaseSolvingStrategyType::MoveMeshFlag )
    .def("MoveMesh", &BaseSolvingStrategyType::MoveMesh )
    .def("Clear", &BaseSolvingStrategyType::Clear )
    //.def("GetModelPart", &BaseSolvingStrategyType::GetModelPart )
    ;
    typedef LapModifiedLinearStrategy< 2, SparseSpaceType, LocalSpaceType, LinearSolverType> LapModifiedLinearStrategy2D;

    py::class_< LapModifiedLinearStrategy2D, LapModifiedLinearStrategy2D::Pointer, BaseSolvingStrategyType >
    (m,"LapModifiedLinearStrategy2D")
     .def(py::init<ModelPart&,BaseSchemeType::Pointer, LinearSolverType::Pointer, bool, bool, bool, bool	>() )
    .def("Solve", &LapModifiedLinearStrategy< 2, SparseSpaceType, LocalSpaceType, LinearSolverType >::Solve )
    ;
    typedef LapModifiedLinearStrategy< 3, SparseSpaceType, LocalSpaceType, LinearSolverType> LapModifiedLinearStrategy3D;

    py::class_< LapModifiedLinearStrategy3D, LapModifiedLinearStrategy3D::Pointer, BaseSolvingStrategyType >
    (m,"LapModifiedLinearStrategy3D")
     .def(py::init<ModelPart&,BaseSchemeType::Pointer, LinearSolverType::Pointer, bool, bool, bool, bool	>() )
    .def("Solve", &LapModifiedLinearStrategy< 3, SparseSpaceType, LocalSpaceType, LinearSolverType >::Solve )
    ;

}

}  // namespace Python.

} // Namespace Kratos
