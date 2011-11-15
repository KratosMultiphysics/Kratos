/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

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
//   Last modified by:    $Author: rrossi $
//   Date:                $Date: 2008-11-10 14:23:33 $
//   Revision:            $Revision: 1.9 $
//
//


// System includes 


// External includes 
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/timer.hpp> 


// Project includes
#include "includes/define.h"
#include "python/add_strategies_to_python.h"
#include "includes/model_part.h"
#include "spaces/ublas_space.h"

//strategies
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"
#include "solving_strategies/strategies/adaptive_residualbased_newton_raphson_strategy.h"
//#include "solving_strategies/strategies/residualbased_arc_lenght_strategy.h"

//schemes
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"

//convergence criterias
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "solving_strategies/convergencecriterias/displacement_criteria.h"
#include "solving_strategies/convergencecriterias/incremental_displacement_criteria.h"
#include "solving_strategies/convergencecriterias/residual_criteria.h"
#include "solving_strategies/convergencecriterias/and_criteria.h"
#include "solving_strategies/convergencecriterias/or_criteria.h"
//#include "solving_strategies/convergencecriterias/and_criteria.h"

//Builder And Solver
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_deactivation.h"


//linear solvers
#include "linear_solvers/linear_solver.h"

//utilities
#include "python/pointer_vector_set_python_interface.h"


namespace Kratos
{

    namespace Python
    {
        using namespace boost::python;



        typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
        typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

        //ADDED BY PAOLO (next two)

        double Dot(SparseSpaceType& dummy, SparseSpaceType::VectorType& rX, SparseSpaceType::VectorType& rY)
        {
            return dummy.Dot(rX, rY);
        }

        void ScaleAndAdd(SparseSpaceType& dummy, const double A, const SparseSpaceType::VectorType& rX, const double B, SparseSpaceType::VectorType& rY)
        //(const double A,const  VectorType& rX, const double B, VectorType& rY) // rY = (A * rX) + (B * rY)
        {
            dummy.ScaleAndAdd(A, rX, B, rY);
        }

        void Mult(SparseSpaceType& dummy, SparseSpaceType::MatrixType& rA, SparseSpaceType::VectorType& rX, SparseSpaceType::VectorType& rY)
        //rY=A*rX (the product is stored inside the rY)
        {
            dummy.Mult(rA, rX, rY);
        }

        void TransposeMult(SparseSpaceType& dummy, SparseSpaceType::MatrixType& rA, SparseSpaceType::VectorType& rX, SparseSpaceType::VectorType& rY)
        //rY=A*rX (the product is stored inside the rY)
        {
            dummy.TransposeMult(rA, rX, rY);
        }

        SparseSpaceType::IndexType Size(SparseSpaceType& dummy, SparseSpaceType::VectorType const& rV)
        {
            return rV.size();
        }

        void ResizeMatrix(SparseSpaceType& dummy, SparseSpaceType::MatrixType& A, unsigned int i1, unsigned int i2)
        {
            dummy.Resize(A, i1, i2);
        }

        void ResizeVector(SparseSpaceType& dummy, SparseSpaceType::VectorType& x, unsigned int i1)
        {
            dummy.Resize(x, i1);
        }

        void SetToZeroMatrix(SparseSpaceType& dummy, SparseSpaceType::MatrixType& A)
        {
            dummy.SetToZero(A);
        }

        void SetToZeroVector(SparseSpaceType& dummy, SparseSpaceType::VectorType& x)
        {
            dummy.SetToZero(x);
        }

        void ClearMatrix(SparseSpaceType& dummy, SparseSpaceType::MatrixPointerType& A)
        {
            dummy.Clear(A);
        }

        void ClearVector(SparseSpaceType& dummy, SparseSpaceType::VectorPointerType& x)
        {
            dummy.Clear(x);
        }

        double TwoNorm(SparseSpaceType& dummy, SparseSpaceType::VectorType& x)
        {
            return dummy.TwoNorm(x);
        }

        void UnaliasedAdd(SparseSpaceType& dummy, SparseSpaceType::VectorType& x, const double A, const SparseSpaceType::VectorType& rY) // x+= a*Y
        {
            dummy.UnaliasedAdd(x, A, rY);
        }

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

        SparseSpaceType::MatrixPointerType CreateEmptyMatrixPointer(SparseSpaceType& dummy)
        {
            return dummy.CreateEmptyMatrixPointer();
        }

        SparseSpaceType::VectorPointerType CreateEmptyVectorPointer(SparseSpaceType& dummy)
        {
            return dummy.CreateEmptyVectorPointer();
        }

        // 	boost::shared_ptr< CompressedMatrix > CreateEmptyMatrixPointer()
        // 	{
        // 		boost::shared_ptr<CompressedMatrix> pNewMat = boost::shared_ptr<CompressedMatrix>(new CompressedMatrix() );
        // 		return pNewMat;
        // 	}
        //
        // 	boost::shared_ptr< Vector > CreateEmptyVectorPointer()
        // 	{
        // 		boost::shared_ptr<Vector > pNewVec = boost::shared_ptr<Vector >(new Vector() );
        // 		return pNewVec;
        // 	}

        CompressedMatrix& GetMatRef(boost::shared_ptr<CompressedMatrix>& dummy)
        {
            return *dummy;
        }

        Vector& GetVecRef(boost::shared_ptr<Vector>& dummy)
        {
            return *dummy;
        }

        void AddStrategiesToPython()
        {
            typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
            typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

            // 			def("CreateEmptyMatrixPointer",CreateEmptyMatrixPointer);
            // 			def("CreateEmptyVectorPointer",CreateEmptyVectorPointer);

            class_< boost::shared_ptr<CompressedMatrix> >("CompressedMatrixPointer", init<boost::shared_ptr<CompressedMatrix> >())
                    .def("GetReference", GetMatRef, return_value_policy<reference_existing_object > ())
                    //    				.def("GetReference", GetRef, return_internal_reference<1>() )
                    ;

            // // // 			class_< CompressedMatrix , boost::noncopyable >("CompressedMatrix", init< >() );


            class_< boost::shared_ptr<Vector> >("VectorPointer", init< boost::shared_ptr<Vector> >())
                    .def("GetReference", GetVecRef, return_value_policy<reference_existing_object > ())
                    ;
            // // // 			class_< Vector , boost::noncopyable >("Vector", init< >() );

            typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
            typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;
            typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;
            typedef ConvergenceCriteria< SparseSpaceType, LocalSpaceType > ::Pointer TConvergenceCriteriaPointer;

            //********************************************************************
            //********************************************************************
            //strategy base class
            class_< BaseSolvingStrategyType, boost::noncopyable > ("SolvingStrategy", init < ModelPart&, bool >())
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
                    //.def("GetModelPart", &BaseSolvingStrategyType::GetModelPart )
                    ;


            class_< ResidualBasedLinearStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >, bases< BaseSolvingStrategyType >, boost::noncopyable >
                    ("ResidualBasedLinearStrategy",
                    init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, bool, bool, bool, bool >())
                    .def("GetResidualNorm", &ResidualBasedLinearStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::GetResidualNorm)
                    ;

            typedef ConvergenceCriteria< SparseSpaceType, LocalSpaceType > TConvergenceCriteriaType;
            typedef BuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > BuilderAndSolverType;

            class_< ResidualBasedNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >, bases< BaseSolvingStrategyType >, boost::noncopyable >
                    ("ResidualBasedNewtonRaphsonStrategy",
                    init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, TConvergenceCriteriaType::Pointer, int, bool, bool, bool
                    >())
                    .def(init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, TConvergenceCriteriaType::Pointer, BuilderAndSolverType::Pointer, int, bool, bool, bool >())
                    .def("SetMaxIterationNumber", &ResidualBasedNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SetMaxIterationNumber)
                    .def("GetMaxIterationNumber", &ResidualBasedNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::GetMaxIterationNumber)
                    .def("SetKeepSystemConstantDuringIterations", &ResidualBasedNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SetKeepSystemConstantDuringIterations)
                    .def("GetKeepSystemConstantDuringIterations", &ResidualBasedNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::GetKeepSystemConstantDuringIterations)
                    ;

            class_< AdaptiveResidualBasedNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >, bases< BaseSolvingStrategyType >, boost::noncopyable >
                    ("AdaptiveResidualBasedNewtonRaphsonStrategy",
                    init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, TConvergenceCriteriaType::Pointer, int, int, bool, bool, bool, double, double, int
                    >())
                    ;

            //********************************************************************
            //********************************************************************
            class_< BaseSchemeType, boost::noncopyable >
                    ("Scheme", init< >())
                    .def("Initialize", &BaseSchemeType::Initialize)
                    .def("SchemeIsInitialized", &BaseSchemeType::SchemeIsInitialized)
                    .def("ElementsAreInitialized", &BaseSchemeType::ElementsAreInitialized)
                    .def("InitializeElements", &BaseSchemeType::InitializeElements)
                    .def("InitializeSolutionStep", &BaseSchemeType::InitializeSolutionStep)
                    .def("FinalizeSolutionStep", &BaseSchemeType::FinalizeSolutionStep)
                    .def("InitializeNonLinIteration", &BaseSchemeType::InitializeNonLinIteration)
                    .def("FinalizeNonLinIteration", &BaseSchemeType::FinalizeNonLinIteration)
                    .def("Predict", &BaseSchemeType::Predict)
                    .def("Update", &BaseSchemeType::Update)
                    .def("CalculateOutputData", &BaseSchemeType::CalculateOutputData)
                    .def("Clean", &BaseSchemeType::Clean)
                    .def("MoveMesh", MoveMesh)
                    .def("Check", &BaseSchemeType::Check)
                    ;

            class_< ResidualBasedIncrementalUpdateStaticScheme< SparseSpaceType, LocalSpaceType>,
                    bases< BaseSchemeType >, boost::noncopyable >
                    (
                    "ResidualBasedIncrementalUpdateStaticScheme", init< >()
                    );


            //********************************************************************
            //********************************************************************
            //convergence criteria base class
            class_< ConvergenceCriteria< SparseSpaceType, LocalSpaceType >, boost::noncopyable > ("ConvergenceCriteria", init<>())
                    .def("SetActualizeRHSFlag", &ConvergenceCriteria<SparseSpaceType, LocalSpaceType >::SetActualizeRHSFlag)
                    .def("GetActualizeRHSflag", &ConvergenceCriteria<SparseSpaceType, LocalSpaceType >::GetActualizeRHSflag)
                    .def("PreCriteria", &ConvergenceCriteria<SparseSpaceType, LocalSpaceType >::PreCriteria)
                    .def("PostCriteria", &ConvergenceCriteria<SparseSpaceType, LocalSpaceType >::PostCriteria)
                    .def("Initialize", &ConvergenceCriteria<SparseSpaceType, LocalSpaceType >::Initialize)
                    .def("InitializeSolutionStep", &ConvergenceCriteria<SparseSpaceType, LocalSpaceType >::InitializeSolutionStep)
                    .def("FinalizeSolutionStep", &ConvergenceCriteria<SparseSpaceType, LocalSpaceType >::FinalizeSolutionStep)
                    .def("Check", &ConvergenceCriteria<SparseSpaceType, LocalSpaceType >::Check)
                    ;

            class_< DisplacementCriteria<SparseSpaceType, LocalSpaceType >,
                    bases<ConvergenceCriteria< SparseSpaceType, LocalSpaceType > >,
                    boost::noncopyable >
                    ("DisplacementCriteria", init< double, double>());
		    
	    class_< IncrementalDisplacementCriteria<SparseSpaceType, LocalSpaceType >,
                    bases<ConvergenceCriteria< SparseSpaceType, LocalSpaceType > >,
                    boost::noncopyable >
                    ("IncrementalDisplacementCriteria", init< double, double>());

            class_<ResidualCriteria<SparseSpaceType, LocalSpaceType >,
                    bases<ConvergenceCriteria< SparseSpaceType, LocalSpaceType > >,
                    boost::noncopyable >
                    ("ResidualCriteria", init< double, double>());

            /*			class_< ResidualCriteria< SparseSpaceType >,
                                             bases<ConvergenceCriteria< SparseSpaceType > >,
                                             boost::noncopyable >
                                            ("ResidualCriteria", init<Model::Pointer, double >() );
			
                                    class_< AndCriteria< SparseSpaceType >,
                                             bases<ConvergenceCriteria< SparseSpaceType > >,
                                             boost::noncopyable >
                                            ("AndCriteria", init<Model::Pointer, ConvergenceCriteria< SparseSpaceType >::Pointer, ConvergenceCriteria< SparseSpaceType >::Pointer >()*/
            //);


            class_<And_Criteria<SparseSpaceType, LocalSpaceType >,
                    bases<ConvergenceCriteria< SparseSpaceType, LocalSpaceType > >,
                    boost::noncopyable >
                    ("AndCriteria", init<TConvergenceCriteriaPointer, TConvergenceCriteriaPointer > ());

            class_<Or_Criteria<SparseSpaceType, LocalSpaceType >,
                    bases<ConvergenceCriteria< SparseSpaceType, LocalSpaceType > >,
                    boost::noncopyable >
                    ("OrCriteria", init<TConvergenceCriteriaPointer, TConvergenceCriteriaPointer > ());

            //********************************************************************
            //********************************************************************
            //

            //Builder and Solver

            class_< BuilderAndSolverType::DofsArrayType, boost::noncopyable > ("DofsArrayType", init<>());

            class_< BuilderAndSolverType, boost::noncopyable > ("BuilderAndSolver", init<LinearSolverType::Pointer > ())
                    .def("SetCalculateReactionsFlag", &BuilderAndSolverType::SetCalculateReactionsFlag)
                    .def("GetCalculateReactionsFlag", &BuilderAndSolverType::GetCalculateReactionsFlag)
                    .def("SetDofSetIsInitializedFlag", &BuilderAndSolverType::SetDofSetIsInitializedFlag)
                    .def("GetDofSetIsInitializedFlag", &BuilderAndSolverType::GetDofSetIsInitializedFlag)
                    .def("SetReshapeMatrixFlag", &BuilderAndSolverType::SetReshapeMatrixFlag)
                    .def("GetReshapeMatrixFlag", &BuilderAndSolverType::GetReshapeMatrixFlag)
                    .def("GetEquationSystemSize", &BuilderAndSolverType::GetEquationSystemSize)
                    .def("BuildLHS", &BuilderAndSolverType::BuildLHS)
                    .def("BuildRHS", &BuilderAndSolverType::BuildRHS)
                    .def("Build", &BuilderAndSolverType::Build)
                    .def("SystemSolve", &BuilderAndSolverType::SystemSolve)
                    .def("BuildAndSolve", &BuilderAndSolverType::BuildAndSolve)
                    .def("BuildRHSAndSolve", &BuilderAndSolverType::BuildRHSAndSolve)
                    .def("ApplyDirichletConditions", &BuilderAndSolverType::ApplyDirichletConditions)
                    .def("SetUpDofSet", &BuilderAndSolverType::SetUpDofSet)
                    .def("GetDofSet", &BuilderAndSolverType::GetDofSet, return_internal_reference<>())
                    .def("SetUpSystem", &BuilderAndSolverType::SetUpSystem)
                    .def("ResizeAndInitializeVectors", &BuilderAndSolverType::ResizeAndInitializeVectors)
                    .def("InitializeSolutionStep", &BuilderAndSolverType::InitializeSolutionStep)
                    .def("FinalizeSolutionStep", &BuilderAndSolverType::FinalizeSolutionStep)
                    .def("CalculateReactions", &BuilderAndSolverType::CalculateReactions)
                    .def("Clear", &BuilderAndSolverType::Clear)
                    .def("Check", &BuilderAndSolverType::Check)
                    .def("SetEchoLevel", &BuilderAndSolverType::SetEchoLevel)
                    .def("GetEchoLevel", &BuilderAndSolverType::GetEchoLevel)
                    ;

            typedef ResidualBasedEliminationBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedEliminationBuilderAndSolverType;

            class_< ResidualBasedEliminationBuilderAndSolverType, bases<BuilderAndSolverType>, boost::noncopyable > ("ResidualBasedEliminationBuilderAndSolver", init< LinearSolverType::Pointer > ());

            typedef ResidualBasedEliminationBuilderAndSolverDeactivation< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedEliminationBuilderAndSolverDeactivationType;

            class_< ResidualBasedEliminationBuilderAndSolverDeactivationType, bases<BuilderAndSolverType>, boost::noncopyable > ("ResidualBasedEliminationBuilderAndSolverDeactivation", init< LinearSolverType::Pointer > ());


            //********************************************************************
            //********************************************************************

            class_< SparseSpaceType, boost::noncopyable > ("UblasSparseSpace", init<>())
                    .def("ClearMatrix", ClearMatrix)
                    .def("ClearVector", ClearVector)
                    .def("ResizeMatrix", ResizeMatrix)
                    .def("ResizeVector", ResizeVector)
                    .def("SetToZeroMatrix", SetToZeroMatrix)
                    .def("SetToZeroVector", SetToZeroVector)
                    .def("TwoNorm", TwoNorm)
                    //the dot product of two vectors
                    .def("Dot", Dot)
                    //the matrix-vector multiplication
                    .def("Mult", Mult)
                    .def("TransposeMult", TransposeMult)
                    .def("Size", Size)
                    .def("UnaliasedAdd", UnaliasedAdd)
                    .def("ScaleAndAdd", ScaleAndAdd)
                    .def("CreateEmptyMatrixPointer", CreateEmptyMatrixPointer)
                    .def("CreateEmptyVectorPointer", CreateEmptyVectorPointer)
                    ;
        }

    } // namespace Python.

} // Namespace Kratos

