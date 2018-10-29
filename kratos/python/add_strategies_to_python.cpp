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
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver_with_constraints.h"

// Linear solvers
#include "linear_solvers/linear_solver.h"


namespace Kratos
{
    namespace Python
    {
        namespace py = pybind11;



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

//             py::class_< Kratos::shared_ptr<CompressedMatrix> >(m,"CompressedMatrixPointer")
//             .def(py::init<Kratos::shared_ptr<CompressedMatrix> >())
//                     .def("GetReference", GetMatRef, py::return_value_policy::reference_internal)
//                     ;
//
//             py::class_< Kratos::shared_ptr<Vector> >(m,"VectorPointer")
//             .def(py::init< Kratos::shared_ptr<Vector> >())
//                     .def("GetReference", GetVecRef, py::return_value_policy::reference_internal)
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

            py::class_< BaseSchemeType, typename BaseSchemeType::Pointer >(m,"Scheme")
                .def(py::init< >())
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

            py::class_< ResidualBasedIncrementalUpdateStaticScheme< SparseSpaceType, LocalSpaceType>,
                typename ResidualBasedIncrementalUpdateStaticScheme< SparseSpaceType, LocalSpaceType>::Pointer,
                BaseSchemeType >
                (m, "ResidualBasedIncrementalUpdateStaticScheme")
                .def(py::init< >()
                );

            typedef typename ResidualBasedIncrementalUpdateStaticSchemeSlip< SparseSpaceType, LocalSpaceType>::RotationToolPointerType RotationToolPointerType;

            py::class_< ResidualBasedIncrementalUpdateStaticSchemeSlip< SparseSpaceType, LocalSpaceType>,
                typename ResidualBasedIncrementalUpdateStaticSchemeSlip< SparseSpaceType, LocalSpaceType>::Pointer,
                ResidualBasedIncrementalUpdateStaticScheme< SparseSpaceType, LocalSpaceType> >
                (m,"ResidualBasedIncrementalUpdateStaticSchemeSlip")
                .def(py::init<unsigned int, unsigned int>())
                .def(py::init<RotationToolPointerType>());

	         // Residual Based Bossak Scheme Type
	         py::class_< ResidualBasedBossakDisplacementSchemeType,
                typename ResidualBasedBossakDisplacementSchemeType::Pointer,
                BaseSchemeType  >
                (m,"ResidualBasedBossakDisplacementScheme")
                .def(py::init< double >() )
                ;

	         // Residual Based Newmark Scheme Type
	         py::class_< ResidualBasedNewmarkDisplacementSchemeType,
                typename ResidualBasedNewmarkDisplacementSchemeType::Pointer,
                BaseSchemeType >(m,"ResidualBasedNewmarkDisplacementScheme")
                .def(py::init< >() )
                ;

	         // Residual Based Pseudo-Static Scheme Type
	         py::class_< ResidualBasedPseudoStaticDisplacementSchemeType,
                typename ResidualBasedPseudoStaticDisplacementSchemeType::Pointer,
                BaseSchemeType >(m,"ResidualBasedPseudoStaticDisplacementScheme")
                .def(py::init< const Variable<double> >() )
                ;

            // Residual Based BDF displacement Scheme Type
            py::class_< ResidualBasedBDFDisplacementSchemeType,typename ResidualBasedBDFDisplacementSchemeType::Pointer, BaseSchemeType  >(m,"ResidualBasedBDFDisplacementScheme")
                .def(py::init<  >() )
                .def(py::init <const std::size_t>())
                ;

            // Residual Based BDF custom Scheme Type
            py::class_< ResidualBasedBDFCustomSchemeType, typename ResidualBasedBDFCustomSchemeType::Pointer, BaseSchemeType  >(m,"ResidualBasedBDFCustomScheme")
                .def(py::init<  >() )
                .def(py::init <const std::size_t>())
                .def(py::init <const std::size_t, Parameters>())
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
            py::class_< ConvergenceCriteriaType,
                ConvergenceCriteriaPointerType >(m,"ConvergenceCriteria")
                .def(py::init<>())
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

            py::class_< DisplacementCriteria<SparseSpaceType, LocalSpaceType >,
                typename DisplacementCriteria< SparseSpaceType, LocalSpaceType >::Pointer,
                ConvergenceCriteriaType >
                (m,"DisplacementCriteria")
                .def(py::init< double, double>())
                ;

            py::class_<ResidualCriteria<SparseSpaceType, LocalSpaceType >,
                typename ResidualCriteria< SparseSpaceType, LocalSpaceType >::Pointer,
                ConvergenceCriteriaType >
                (m,"ResidualCriteria")
                .def(py::init< double, double>())
                ;

            py::class_<And_Criteria<SparseSpaceType, LocalSpaceType >,
                typename And_Criteria< SparseSpaceType, LocalSpaceType >::Pointer,
                ConvergenceCriteriaType >
                (m,"AndCriteria")
                .def(py::init<ConvergenceCriteriaPointerType, ConvergenceCriteriaPointerType > ())
                ;

            py::class_<Or_Criteria<SparseSpaceType, LocalSpaceType >,
                typename Or_Criteria< SparseSpaceType, LocalSpaceType >::Pointer,
                ConvergenceCriteriaType >
                (m,"OrCriteria")
                .def(py::init<ConvergenceCriteriaPointerType, ConvergenceCriteriaPointerType > ())
                ;

            //********************************************************************
            //********************************************************************
            //********************************************************************
            //********************************************************************
            //********************************************************************
            //********************************************************************
            //Builder and Solver
            typedef BuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > BuilderAndSolverType;


            py::class_< BuilderAndSolverType::DofsArrayType, BuilderAndSolverType::DofsArrayType::Pointer>(m,"DofsArrayType")
                .def(py::init<>());

            py::class_< BuilderAndSolverType, typename BuilderAndSolverType::Pointer>(m,"BuilderAndSolver")
            .def(py::init<LinearSolverType::Pointer > ())
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
                .def("GetDofSet", &BuilderAndSolverType::GetDofSet, py::return_value_policy::reference_internal)
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
            py::class_< ResidualBasedEliminationBuilderAndSolverType, ResidualBasedEliminationBuilderAndSolverType::Pointer, BuilderAndSolverType>(m,"ResidualBasedEliminationBuilderAndSolver")
                .def(py::init< LinearSolverType::Pointer > ());

            typedef ResidualBasedBlockBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedBlockBuilderAndSolverType;
            py::class_< ResidualBasedBlockBuilderAndSolverType, ResidualBasedBlockBuilderAndSolverType::Pointer,BuilderAndSolverType>(m,"ResidualBasedBlockBuilderAndSolver")
                .def(py::init< LinearSolverType::Pointer > ());

            typedef ResidualBasedBlockBuilderAndSolverWithConstraints< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedBlockBuilderAndSolverWithConstraintsType;
            py::class_< ResidualBasedBlockBuilderAndSolverWithConstraintsType, ResidualBasedBlockBuilderAndSolverWithConstraintsType::Pointer,ResidualBasedBlockBuilderAndSolverType>(m,"ResidualBasedBlockBuilderAndSolverWithConstraints")
                .def(py::init< LinearSolverType::Pointer > ());

            //********************************************************************
            //********************************************************************
            //********************************************************************
            //********************************************************************
            //********************************************************************
            //********************************************************************

            py::class_< SparseSpaceType>(m,"UblasSparseSpace")
                .def(py::init<>())
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

            py::class_< BaseSolvingStrategyType, typename BaseSolvingStrategyType::Pointer >(m,"SolvingStrategy")
                .def(py::init < ModelPart&, bool >())
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

            py::class_< ResidualBasedLinearStrategyType, typename ResidualBasedLinearStrategyType::Pointer,BaseSolvingStrategyType >
                (m,"ResidualBasedLinearStrategy")
                .def(py::init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, bool, bool, bool, bool >())
                .def(py::init < ModelPart& ,  BaseSchemeType::Pointer, LinearSolverType::Pointer, BuilderAndSolverType::Pointer, bool, bool, bool,  bool  >())
                .def("GetResidualNorm", &ResidualBasedLinearStrategyType::GetResidualNorm)
                .def("SetBuilderAndSolver", &ResidualBasedLinearStrategyType::SetBuilderAndSolver)
                .def("GetSystemMatrix", &ResidualBasedLinearStrategyType::GetSystemMatrix, py::return_value_policy::reference_internal)
                .def("GetSystemVector", &ResidualBasedLinearStrategyType::GetSystemVector, py::return_value_policy::reference_internal)
                .def("GetSolutionVector", &ResidualBasedLinearStrategyType::GetSolutionVector, py::return_value_policy::reference_internal)
                ;

            typedef ResidualBasedNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedNewtonRaphsonStrategyType;

            py::class_< ResidualBasedNewtonRaphsonStrategyType, typename ResidualBasedNewtonRaphsonStrategyType::Pointer, BaseSolvingStrategyType >
                (m,"ResidualBasedNewtonRaphsonStrategy")
                .def(py::init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, int, bool, bool, bool >())
                .def(py::init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, BuilderAndSolverType::Pointer, int, bool, bool, bool >())
                .def("SetMaxIterationNumber", &ResidualBasedNewtonRaphsonStrategyType::SetMaxIterationNumber)
                .def("GetMaxIterationNumber", &ResidualBasedNewtonRaphsonStrategyType::GetMaxIterationNumber)
                .def("SetKeepSystemConstantDuringIterations", &ResidualBasedNewtonRaphsonStrategyType::SetKeepSystemConstantDuringIterations)
                .def("GetKeepSystemConstantDuringIterations", &ResidualBasedNewtonRaphsonStrategyType::GetKeepSystemConstantDuringIterations)
                .def("SetInitializePerformedFlag", &ResidualBasedNewtonRaphsonStrategyType::SetInitializePerformedFlag)
                .def("GetInitializePerformedFlag", &ResidualBasedNewtonRaphsonStrategyType::GetInitializePerformedFlag)
                .def("GetSystemMatrix", &ResidualBasedNewtonRaphsonStrategyType::GetSystemMatrix, py::return_value_policy::reference_internal)
                .def("GetSystemVector", &ResidualBasedNewtonRaphsonStrategyType::GetSystemVector, py::return_value_policy::reference_internal)
                .def("GetSolutionVector", &ResidualBasedNewtonRaphsonStrategyType::GetSolutionVector, py::return_value_policy::reference_internal)
                ;

            py::class_< AdaptiveResidualBasedNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >,
                typename AdaptiveResidualBasedNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::Pointer,
                BaseSolvingStrategyType >
                (m,"AdaptiveResidualBasedNewtonRaphsonStrategy")
                .def(py::init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, int, int, bool, bool, bool, double, double, int
                >())
                ;

            py::class_< LineSearchStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >,
                typename LineSearchStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::Pointer,
                ResidualBasedNewtonRaphsonStrategyType  >
                (m,"LineSearchStrategy")
                .def(py::init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, int, bool, bool, bool >())
                .def(py::init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, BuilderAndSolverType::Pointer, int, bool, bool, bool >())
                ;

            py::class_< ExplicitStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >,
                typename ExplicitStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::Pointer,
                BaseSolvingStrategyType >(m,"Explicit_Strategy")
                .def(py::init<ModelPart&, int, bool >() )
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
