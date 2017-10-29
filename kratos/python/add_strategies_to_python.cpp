//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand, Riccardo Rossi
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

// Strategies
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"
#include "solving_strategies/strategies/adaptive_residualbased_newton_raphson_strategy.h"
#include "solving_strategies/strategies/line_search_strategy.h"
#include "solving_strategies/strategies/explicit_strategy.h"
//#include "solving_strategies/strategies/residualbased_arc_lenght_strategy.h"

// Schemes
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme_slip.h"
#include "solving_strategies/schemes/residual_based_bossak_displacement_scheme.hpp"
#include "solving_strategies/schemes/residual_based_newmark_displacement_scheme.hpp"

// Convergence criterias
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "solving_strategies/convergencecriterias/displacement_criteria.h"
#include "solving_strategies/convergencecriterias/residual_criteria.h"
#include "solving_strategies/convergencecriterias/and_criteria.h"
#include "solving_strategies/convergencecriterias/or_criteria.h"

// Builder And Solver
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"

// Linear solvers
#include "linear_solvers/linear_solver.h"

// Utilities
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

        SparseSpaceType::IndexType Size1(SparseSpaceType& dummy, SparseSpaceType::MatrixType const& rM)
        {
            return rM.size1();
        }

        SparseSpaceType::IndexType Size2(SparseSpaceType& dummy, SparseSpaceType::MatrixType const& rM)
        {
            return rM.size2();
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
            int numNodes = static_cast<int>(rNodes.size());
            
            #pragma omp parallel for
            for(int i = 0; i < numNodes; i++)  
            {
                auto itNode = rNodes.begin() + i;

                noalias(itNode->Coordinates()) = itNode->GetInitialPosition().Coordinates();
                noalias(itNode->Coordinates()) += itNode->FastGetSolutionStepValue(DISPLACEMENT);
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
	  //typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType; //already done up in this file
	  //typedef UblasSpace<double, Matrix, Vector> LocalSpaceType; //already done up in this file

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
                    .def("Initialize", &BaseSolvingStrategyType::Initialize)
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
					.def("InitializeSolutionStep", &BaseSolvingStrategyType::InitializeSolutionStep)
					.def("FinalizeSolutionStep", &BaseSolvingStrategyType::FinalizeSolutionStep)
					.def("SolveSolutionStep", &BaseSolvingStrategyType::SolveSolutionStep)
                    //.def("GetModelPart", &BaseSolvingStrategyType::GetModelPart )
                    ;

            typedef ConvergenceCriteria< SparseSpaceType, LocalSpaceType > TConvergenceCriteriaType;
            typedef BuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > BuilderAndSolverType;

            class_< ResidualBasedLinearStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >, bases< BaseSolvingStrategyType >, boost::noncopyable >
                    ("ResidualBasedLinearStrategy",init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, bool, bool, bool, bool >())
                    .def(init < ModelPart& ,  BaseSchemeType::Pointer, LinearSolverType::Pointer, BuilderAndSolverType::Pointer, bool, bool, bool,  bool  >())
                    .def("GetResidualNorm", &ResidualBasedLinearStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::GetResidualNorm)
                    .def("SetBuilderAndSolver", &ResidualBasedLinearStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SetBuilderAndSolver)
                    .def("GetSystemMatrix", &ResidualBasedLinearStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::GetSystemMatrix, return_internal_reference<>())
                    ;

            typedef ResidualBasedNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedNewtonRaphsonStrategyType;

            class_< ResidualBasedNewtonRaphsonStrategyType, bases< BaseSolvingStrategyType >, boost::noncopyable >
                    ("ResidualBasedNewtonRaphsonStrategy", init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, TConvergenceCriteriaType::Pointer, int, bool, bool, bool >())
                    .def(init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, TConvergenceCriteriaType::Pointer, BuilderAndSolverType::Pointer, int, bool, bool, bool >())
                    .def("SetMaxIterationNumber", &ResidualBasedNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SetMaxIterationNumber)
                    .def("GetMaxIterationNumber", &ResidualBasedNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::GetMaxIterationNumber)
                    .def("SetKeepSystemConstantDuringIterations", &ResidualBasedNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SetKeepSystemConstantDuringIterations)
                    .def("GetKeepSystemConstantDuringIterations", &ResidualBasedNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::GetKeepSystemConstantDuringIterations)
                    .def("SetInitializePerformedFlag", &ResidualBasedNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SetInitializePerformedFlag)
                    .def("GetInitializePerformedFlag", &ResidualBasedNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::GetInitializePerformedFlag)
                    .def("GetSystemMatrix", &ResidualBasedNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::GetSystemMatrix, return_internal_reference<>())
                    ;

            class_< AdaptiveResidualBasedNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >, bases< BaseSolvingStrategyType >, boost::noncopyable >
                    ("AdaptiveResidualBasedNewtonRaphsonStrategy",
                    init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, TConvergenceCriteriaType::Pointer, int, int, bool, bool, bool, double, double, int
                    >())
                    ;

            class_< LineSearchStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >, bases< ResidualBasedNewtonRaphsonStrategyType >, boost::noncopyable >
                    ("LineSearchStrategy", init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, TConvergenceCriteriaType::Pointer, int, bool, bool, bool >())
                    .def(init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, TConvergenceCriteriaType::Pointer, BuilderAndSolverType::Pointer, int, bool, bool, bool >())
                    ;

            class_< ExplicitStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >,
                    bases< BaseSolvingStrategyType >,  boost::noncopyable >
                    ("Explicit_Strategy",
                    init<ModelPart&, int, bool >() )
                    //AssembleLoop loops the elements calling AddExplicitContribution. Using processinfo the element is the one who "decides" which variable to modify.
                    .def("AssembleLoop",&ExplicitStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::AssembleLoop)
                    //once the assembleloop has been performed, the variable must be normalized. (for example with the nodal mass or the nodal area). Loop on nodes.
                    .def("NormalizeVariable",&ExplicitStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::NormalizeVariable)
                    //ExplicitUpdateLoop modifies a vectorial variable by adding another variable (the RHS, PRESS_PROJ,etc) multiplied by a user-given factor (ie delta_time)
                    .def("ExplicitUpdateLoop",&ExplicitStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::ExplicitUpdateLoop)
                    //initialize and finalize.
                    .def("InitializeSolutionStep",&ExplicitStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::InitializeSolutionStep)
                    .def("FinalizeSolutionStep",&ExplicitStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::FinalizeSolutionStep)
                    ;

            //********************************************************************
            //********************************************************************

	    typedef ResidualBasedBossakDisplacementScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedBossakDisplacementSchemeType;
	    typedef ResidualBasedNewmarkDisplacementScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedNewmarkDisplacementSchemeType;

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

            class_< ResidualBasedIncrementalUpdateStaticScheme< SparseSpaceType, LocalSpaceType>,
                    bases< BaseSchemeType >, boost::noncopyable >
                    (
                    "ResidualBasedIncrementalUpdateStaticScheme", init< >()
                    );

            class_< ResidualBasedIncrementalUpdateStaticSchemeSlip< SparseSpaceType, LocalSpaceType>,
                    bases< ResidualBasedIncrementalUpdateStaticScheme< SparseSpaceType, LocalSpaceType> >,
                    boost::noncopyable >
                    ("ResidualBasedIncrementalUpdateStaticSchemeSlip", init<unsigned int, unsigned int>());

	    // Residual Based Bossak Scheme Type
	    class_< ResidualBasedBossakDisplacementSchemeType,
            bases< BaseSchemeType >,  boost::noncopyable >
            (
                "ResidualBasedBossakDisplacementScheme", init< double >() )
            .def("Initialize", &ResidualBasedBossakDisplacementScheme<SparseSpaceType, LocalSpaceType>::Initialize)
            ;

	    // Residual Based Newmark Scheme Type
	    class_< ResidualBasedNewmarkDisplacementSchemeType,
            bases< BaseSchemeType >,  boost::noncopyable >
            (
                "ResidualBasedNewmarkDisplacementScheme", init< >() )
            .def("Initialize", &ResidualBasedNewmarkDisplacementScheme<SparseSpaceType, LocalSpaceType>::Initialize)
            ;

            //********************************************************************
            //********************************************************************
            // Convergence criteria base class
            class_< ConvergenceCriteria< SparseSpaceType, LocalSpaceType >, boost::noncopyable > ("ConvergenceCriteria", init<>())
                    .def("SetActualizeRHSFlag", &ConvergenceCriteria<SparseSpaceType, LocalSpaceType >::SetActualizeRHSFlag)
                    .def("GetActualizeRHSflag", &ConvergenceCriteria<SparseSpaceType, LocalSpaceType >::GetActualizeRHSflag)
                    .def("PreCriteria", &ConvergenceCriteria<SparseSpaceType, LocalSpaceType >::PreCriteria)
                    .def("PostCriteria", &ConvergenceCriteria<SparseSpaceType, LocalSpaceType >::PostCriteria)
                    .def("Initialize", &ConvergenceCriteria<SparseSpaceType, LocalSpaceType >::Initialize)
                    .def("InitializeNonLinearIteration", &ConvergenceCriteria<SparseSpaceType, LocalSpaceType >::InitializeNonLinearIteration)
                    .def("InitializeSolutionStep", &ConvergenceCriteria<SparseSpaceType, LocalSpaceType >::InitializeSolutionStep)
                    .def("FinalizeNonLinearIteration", &ConvergenceCriteria<SparseSpaceType, LocalSpaceType >::FinalizeNonLinearIteration)
                    .def("FinalizeSolutionStep", &ConvergenceCriteria<SparseSpaceType, LocalSpaceType >::FinalizeSolutionStep)
                    .def("Check", &ConvergenceCriteria<SparseSpaceType, LocalSpaceType >::Check)
                    .def("SetEchoLevel", &ConvergenceCriteria<SparseSpaceType, LocalSpaceType >::SetEchoLevel)
                    ;

            class_< DisplacementCriteria<SparseSpaceType, LocalSpaceType >,
                    bases<ConvergenceCriteria< SparseSpaceType, LocalSpaceType > >,
                    boost::noncopyable >
                    ("DisplacementCriteria", init< double, double>())
                    .def("SetEchoLevel",&ResidualCriteria<SparseSpaceType, LocalSpaceType >::SetEchoLevel)
                    .def("SetActualizeRHSFlag",&ResidualCriteria<SparseSpaceType, LocalSpaceType >::SetActualizeRHSFlag)
                    ;

            class_<ResidualCriteria<SparseSpaceType, LocalSpaceType >,
                    bases<ConvergenceCriteria< SparseSpaceType, LocalSpaceType > >,
                    boost::noncopyable >
                    ("ResidualCriteria", init< double, double>())
					.def("SetEchoLevel",&ResidualCriteria<SparseSpaceType, LocalSpaceType >::SetEchoLevel)
					.def("SetActualizeRHSFlag",&ResidualCriteria<SparseSpaceType, LocalSpaceType >::SetActualizeRHSFlag)
					;

            class_<And_Criteria<SparseSpaceType, LocalSpaceType >,
                    bases<ConvergenceCriteria< SparseSpaceType, LocalSpaceType > >,
                    boost::noncopyable >
                    ("AndCriteria", init<TConvergenceCriteriaPointer, TConvergenceCriteriaPointer > ())
                    .def("SetEchoLevel",&ResidualCriteria<SparseSpaceType, LocalSpaceType >::SetEchoLevel)
                    .def("SetActualizeRHSFlag",&ResidualCriteria<SparseSpaceType, LocalSpaceType >::SetActualizeRHSFlag)
                    ;

            class_<Or_Criteria<SparseSpaceType, LocalSpaceType >,
                    bases<ConvergenceCriteria< SparseSpaceType, LocalSpaceType > >,
                    boost::noncopyable >
                    ("OrCriteria", init<TConvergenceCriteriaPointer, TConvergenceCriteriaPointer > ())
                    .def("SetEchoLevel",&ResidualCriteria<SparseSpaceType, LocalSpaceType >::SetEchoLevel)
                    .def("SetActualizeRHSFlag",&ResidualCriteria<SparseSpaceType, LocalSpaceType >::SetActualizeRHSFlag)
                    ;

            //********************************************************************
            //********************************************************************

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

            typedef ResidualBasedBlockBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedBlockBuilderAndSolverType;
            class_< ResidualBasedBlockBuilderAndSolverType, bases<BuilderAndSolverType>, boost::noncopyable > ("ResidualBasedBlockBuilderAndSolver", init< LinearSolverType::Pointer > ());

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
                    .def("Size1", Size1)
                    .def("Size2", Size2)
                    .def("UnaliasedAdd", UnaliasedAdd)
                    .def("ScaleAndAdd", ScaleAndAdd)
                    .def("CreateEmptyMatrixPointer", CreateEmptyMatrixPointer)
                    .def("CreateEmptyVectorPointer", CreateEmptyVectorPointer)
                    ;
        }

    } // namespace Python.

} // Namespace Kratos
