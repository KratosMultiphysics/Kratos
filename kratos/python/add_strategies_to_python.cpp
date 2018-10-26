//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand, Riccardo Rossi
//

// System includes


// External includes
#ifdef KRATOS_USE_AMATRIX   // This macro definition is for the migration period and to be removed afterward please do not use it 
#include "boost/numeric/ublas/matrix.hpp" // for the sparse space dense vector
#else
#endif // KRATOS_USE_AMATRIX

// Project includes
#include "includes/define_python.h"
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
#include "solving_strategies/schemes/residual_based_pseudo_static_displacement_scheme.h"
#include "solving_strategies/schemes/residual_based_bdf_displacement_scheme.h"
#include "solving_strategies/schemes/residual_based_bdf_custom_scheme.h"

// Convergence criterias
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "solving_strategies/convergencecriterias/displacement_criteria.h"
#include "solving_strategies/convergencecriterias/residual_criteria.h"
#include "solving_strategies/convergencecriterias/and_criteria.h"
#include "solving_strategies/convergencecriterias/or_criteria.h"

// Builder And Solver
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver_with_constraints.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"

// Linear solvers
#include "linear_solvers/linear_solver.h"


namespace Kratos
{
    namespace Python
    {
        using namespace pybind11;



        typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>> SparseSpaceType;
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

        CompressedMatrix& GetMatRef(Kratos::shared_ptr<CompressedMatrix>& dummy)
        {
            return *dummy;
        }

        Vector& GetVecRef(Kratos::shared_ptr<Vector>& dummy)
        {
            return *dummy;
        }

        void AddStrategiesToPython(pybind11::module& m)
        {

//             class_< Kratos::shared_ptr<CompressedMatrix> >(m,"CompressedMatrixPointer")
//             .def(init<Kratos::shared_ptr<CompressedMatrix> >())
//                     .def("GetReference", GetMatRef, return_value_policy::reference_internal)
//                     ;
//
//             class_< Kratos::shared_ptr<Vector> >(m,"VectorPointer")
//             .def(init< Kratos::shared_ptr<Vector> >())
//                     .def("GetReference", GetVecRef, return_value_policy::reference_internal)
//                     ;

            typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;





            //********************************************************************
            //********************************************************************
            typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;
            typedef ResidualBasedBossakDisplacementScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedBossakDisplacementSchemeType;
            typedef ResidualBasedNewmarkDisplacementScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedNewmarkDisplacementSchemeType;
            typedef ResidualBasedPseudoStaticDisplacementScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedPseudoStaticDisplacementSchemeType;
            typedef ResidualBasedBDFDisplacementScheme< SparseSpaceType, LocalSpaceType > ResidualBasedBDFDisplacementSchemeType;
            typedef ResidualBasedBDFCustomScheme< SparseSpaceType, LocalSpaceType > ResidualBasedBDFCustomSchemeType;

            class_< BaseSchemeType, typename BaseSchemeType::Pointer >(m,"Scheme")
            .def(init< >())
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
                    typename ResidualBasedIncrementalUpdateStaticScheme< SparseSpaceType, LocalSpaceType>::Pointer,
                    BaseSchemeType >
                    (m, "ResidualBasedIncrementalUpdateStaticScheme")
                    .def(init< >()
                    );


            typedef typename ResidualBasedIncrementalUpdateStaticSchemeSlip< SparseSpaceType, LocalSpaceType>::RotationToolPointerType RotationToolPointerType;



            class_< ResidualBasedIncrementalUpdateStaticSchemeSlip< SparseSpaceType, LocalSpaceType>,
                    typename ResidualBasedIncrementalUpdateStaticSchemeSlip< SparseSpaceType, LocalSpaceType>::Pointer,
                    ResidualBasedIncrementalUpdateStaticScheme< SparseSpaceType, LocalSpaceType> >
                    (m,"ResidualBasedIncrementalUpdateStaticSchemeSlip")
                    .def(init<unsigned int, unsigned int>())
                    .def(init<RotationToolPointerType>());


	         // Residual Based Bossak Scheme Type
	         class_< ResidualBasedBossakDisplacementSchemeType,
                    typename ResidualBasedBossakDisplacementSchemeType::Pointer,
                    BaseSchemeType  >
                    (m,"ResidualBasedBossakDisplacementScheme")
                    .def(init< double >() )
                    ;

	         // Residual Based Newmark Scheme Type
	         class_< ResidualBasedNewmarkDisplacementSchemeType,
                   typename ResidualBasedNewmarkDisplacementSchemeType::Pointer,
                   BaseSchemeType >(m,"ResidualBasedNewmarkDisplacementScheme")
                   .def(init< >() )
                   ;

	         // Residual Based Pseudo-Static Scheme Type
	         class_< ResidualBasedPseudoStaticDisplacementSchemeType,
                   typename ResidualBasedPseudoStaticDisplacementSchemeType::Pointer,
                   BaseSchemeType >(m,"ResidualBasedPseudoStaticDisplacementScheme")
                   .def(init< const Variable<double> >() )
                   ;

            // Residual Based BDF displacement Scheme Type
            class_< ResidualBasedBDFDisplacementSchemeType,typename ResidualBasedBDFDisplacementSchemeType::Pointer, BaseSchemeType  >(m,"ResidualBasedBDFDisplacementScheme")
                .def(init<  >() )
                .def(init <const std::size_t>())
            ;

            // Residual Based BDF custom Scheme Type
            class_< ResidualBasedBDFCustomSchemeType, typename ResidualBasedBDFCustomSchemeType::Pointer, BaseSchemeType  >(m,"ResidualBasedBDFCustomScheme")
                .def(init<  >() )
                .def(init <const std::size_t>())
                .def(init <const std::size_t, Parameters>())
            ;

            //********************************************************************
            //********************************************************************
            //********************************************************************
            //********************************************************************
            //********************************************************************
            //********************************************************************
            typedef ConvergenceCriteria<SparseSpaceType, LocalSpaceType> ConvergenceCriteriaType;
            typedef typename ConvergenceCriteriaType::Pointer ConvergenceCriteriaPointerType;

            // Convergence criteria base class
            class_< ConvergenceCriteriaType,
                    ConvergenceCriteriaPointerType >(m,"ConvergenceCriteria")
                    .def(init<>())
                    .def("SetActualizeRHSFlag", &ConvergenceCriteriaType::SetActualizeRHSFlag)
                    .def("GetActualizeRHSflag", &ConvergenceCriteriaType::GetActualizeRHSflag)
                    .def("PreCriteria", &ConvergenceCriteriaType::PreCriteria)
                    .def("PostCriteria", &ConvergenceCriteriaType::PostCriteria)
                    .def("Initialize", &ConvergenceCriteriaType::Initialize)
                    .def("InitializeNonLinearIteration", &ConvergenceCriteriaType::InitializeNonLinearIteration)
                    .def("InitializeSolutionStep", &ConvergenceCriteriaType::InitializeSolutionStep)
                    .def("FinalizeNonLinearIteration", &ConvergenceCriteriaType::FinalizeNonLinearIteration)
                    .def("FinalizeSolutionStep", &ConvergenceCriteriaType::FinalizeSolutionStep)
                    .def("Check", &ConvergenceCriteriaType::Check)
                    .def("SetEchoLevel", &ConvergenceCriteriaType::SetEchoLevel)
                    ;

            class_< DisplacementCriteria<SparseSpaceType, LocalSpaceType >,
                    typename DisplacementCriteria< SparseSpaceType, LocalSpaceType >::Pointer,
                    ConvergenceCriteriaType >
                    (m,"DisplacementCriteria")
                    .def(init< double, double>())
                    ;

            class_<ResidualCriteria<SparseSpaceType, LocalSpaceType >,
                    typename ResidualCriteria< SparseSpaceType, LocalSpaceType >::Pointer,
                    ConvergenceCriteriaType >
                    (m,"ResidualCriteria")
                    .def(init< double, double>())
                    ;

            class_<And_Criteria<SparseSpaceType, LocalSpaceType >,
                    typename And_Criteria< SparseSpaceType, LocalSpaceType >::Pointer,
                    ConvergenceCriteriaType >
                    (m,"AndCriteria")
                    .def(init<ConvergenceCriteriaPointerType, ConvergenceCriteriaPointerType > ())
                    ;

            class_<Or_Criteria<SparseSpaceType, LocalSpaceType >,
                    typename Or_Criteria< SparseSpaceType, LocalSpaceType >::Pointer,
                    ConvergenceCriteriaType >
                    (m,"OrCriteria")
                    .def(init<ConvergenceCriteriaPointerType, ConvergenceCriteriaPointerType > ())
                    ;

            //********************************************************************
            //********************************************************************
            //********************************************************************
            //********************************************************************
            //********************************************************************
            //********************************************************************
            //Builder and Solver
            typedef BuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > BuilderAndSolverType;


            class_< BuilderAndSolverType::DofsArrayType, BuilderAndSolverType::DofsArrayType::Pointer>(m,"DofsArrayType")
            .def(init<>());

            class_< BuilderAndSolverType, typename BuilderAndSolverType::Pointer>(m,"BuilderAndSolver")
            .def(init<LinearSolverType::Pointer > ())
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
                    .def("GetDofSet", &BuilderAndSolverType::GetDofSet, return_value_policy::reference_internal)
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
            class_< ResidualBasedEliminationBuilderAndSolverType, ResidualBasedEliminationBuilderAndSolverType::Pointer, BuilderAndSolverType>(m,"ResidualBasedEliminationBuilderAndSolver")
            .def(init< LinearSolverType::Pointer > ());

            typedef ResidualBasedBlockBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedBlockBuilderAndSolverType;
            class_< ResidualBasedBlockBuilderAndSolverType, ResidualBasedBlockBuilderAndSolverType::Pointer,BuilderAndSolverType>(m,"ResidualBasedBlockBuilderAndSolver")
            .def(init< LinearSolverType::Pointer > ());

            typedef ResidualBasedBlockBuilderAndSolverWithConstraints< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedBlockBuilderAndSolverWithConstraintsType;
            class_< ResidualBasedBlockBuilderAndSolverWithConstraintsType, ResidualBasedBlockBuilderAndSolverWithConstraintsType::Pointer,ResidualBasedBlockBuilderAndSolverType>(m,"ResidualBasedBlockBuilderAndSolverWithConstraints")
            .def(init< LinearSolverType::Pointer > ());

            //********************************************************************
            //********************************************************************
            //********************************************************************
            //********************************************************************
            //********************************************************************
            //********************************************************************

            class_< SparseSpaceType>(m,"UblasSparseSpace")
                    .def(init<>())
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

            //********************************************************************
            //********************************************************************
            //********************************************************************
            //********************************************************************
            //********************************************************************
            //********************************************************************
            //strategy base class
            typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;

            class_< BaseSolvingStrategyType, typename BaseSolvingStrategyType::Pointer >(m,"SolvingStrategy")
            .def(init < ModelPart&, bool >())
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


            typedef ResidualBasedLinearStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedLinearStrategyType;

            class_< ResidualBasedLinearStrategyType, typename ResidualBasedLinearStrategyType::Pointer,BaseSolvingStrategyType >
                    (m,"ResidualBasedLinearStrategy")
                    .def(init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, bool, bool, bool, bool >())
                    .def(init < ModelPart& ,  BaseSchemeType::Pointer, LinearSolverType::Pointer, BuilderAndSolverType::Pointer, bool, bool, bool,  bool  >())
                    .def("GetResidualNorm", &ResidualBasedLinearStrategyType::GetResidualNorm)
                    .def("SetBuilderAndSolver", &ResidualBasedLinearStrategyType::SetBuilderAndSolver)
                    .def("GetSystemMatrix", &ResidualBasedLinearStrategyType::GetSystemMatrix, return_value_policy::reference_internal)
                    .def("GetSystemVector", &ResidualBasedLinearStrategyType::GetSystemVector, return_value_policy::reference_internal)
                    .def("GetSolutionVector", &ResidualBasedLinearStrategyType::GetSolutionVector, return_value_policy::reference_internal)
                    ;

            typedef ResidualBasedNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedNewtonRaphsonStrategyType;

            class_< ResidualBasedNewtonRaphsonStrategyType, typename ResidualBasedNewtonRaphsonStrategyType::Pointer, BaseSolvingStrategyType >
                    (m,"ResidualBasedNewtonRaphsonStrategy")
                    .def(init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, int, bool, bool, bool >())
                    .def(init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, BuilderAndSolverType::Pointer, int, bool, bool, bool >())
                    .def("SetMaxIterationNumber", &ResidualBasedNewtonRaphsonStrategyType::SetMaxIterationNumber)
                    .def("GetMaxIterationNumber", &ResidualBasedNewtonRaphsonStrategyType::GetMaxIterationNumber)
                    .def("SetKeepSystemConstantDuringIterations", &ResidualBasedNewtonRaphsonStrategyType::SetKeepSystemConstantDuringIterations)
                    .def("GetKeepSystemConstantDuringIterations", &ResidualBasedNewtonRaphsonStrategyType::GetKeepSystemConstantDuringIterations)
                    .def("SetInitializePerformedFlag", &ResidualBasedNewtonRaphsonStrategyType::SetInitializePerformedFlag)
                    .def("GetInitializePerformedFlag", &ResidualBasedNewtonRaphsonStrategyType::GetInitializePerformedFlag)
                    .def("GetSystemMatrix", &ResidualBasedNewtonRaphsonStrategyType::GetSystemMatrix, return_value_policy::reference_internal)
                    .def("GetSystemVector", &ResidualBasedNewtonRaphsonStrategyType::GetSystemVector, return_value_policy::reference_internal)
                    .def("GetSolutionVector", &ResidualBasedNewtonRaphsonStrategyType::GetSolutionVector, return_value_policy::reference_internal)
                    ;

            class_< AdaptiveResidualBasedNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >,
                    typename AdaptiveResidualBasedNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::Pointer,
                    BaseSolvingStrategyType >
                    (m,"AdaptiveResidualBasedNewtonRaphsonStrategy")
                    .def(init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, int, int, bool, bool, bool, double, double, int
                    >())
                    ;

            class_< LineSearchStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >,
                    typename LineSearchStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::Pointer,
                    ResidualBasedNewtonRaphsonStrategyType  >
                    (m,"LineSearchStrategy")
                    .def(init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, int, bool, bool, bool >())
                    .def(init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, BuilderAndSolverType::Pointer, int, bool, bool, bool >())
                    ;

            class_< ExplicitStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >,
                    typename ExplicitStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::Pointer,
                    BaseSolvingStrategyType >(m,"Explicit_Strategy")
                    .def(init<ModelPart&, int, bool >() )
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

        }

    } // namespace Python.

} // Namespace Kratos
