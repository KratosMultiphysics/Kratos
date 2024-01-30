//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand, Riccardo Rossi
//

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "python/add_strategies_to_python.h"
#include "includes/model_part.h"
#include "spaces/ublas_space.h"
#include "includes/ublas_complex_interface.h"
#include "utilities/variable_utils.h"

// Strategies
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/implicit_solving_strategy.h"
#include "solving_strategies/strategies/explicit_solving_strategy.h"
#include "solving_strategies/strategies/explicit_solving_strategy_runge_kutta.h"
#include "solving_strategies/strategies/explicit_solving_strategy_bfecc.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"
#include "solving_strategies/strategies/adaptive_residualbased_newton_raphson_strategy.h"
#include "solving_strategies/strategies/line_search_strategy.h"
#include "solving_strategies/strategies/arc_length_strategy.h"

// Schemes
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme_slip.h"
#include "solving_strategies/schemes/residual_based_bossak_displacement_scheme.hpp"
#include "solving_strategies/schemes/residual_based_newmark_displacement_scheme.hpp"
#include "solving_strategies/schemes/residual_based_pseudo_static_displacement_scheme.h"
#include "solving_strategies/schemes/residual_based_bdf_displacement_scheme.h"
#include "solving_strategies/schemes/residual_based_bdf_custom_scheme.h"
#include "solving_strategies/schemes/residual_based_adjoint_static_scheme.h"
#include "solving_strategies/schemes/residual_based_adjoint_steady_scheme.h"
#include "solving_strategies/schemes/residual_based_adjoint_bossak_scheme.h"

// sensitivity builder schemes
#include "solving_strategies/schemes/sensitivity_builder_scheme.h"

// Convergence criterias
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "solving_strategies/convergencecriterias/displacement_criteria.h"
#include "solving_strategies/convergencecriterias/residual_criteria.h"
#include "solving_strategies/convergencecriterias/and_criteria.h"
#include "solving_strategies/convergencecriterias/or_criteria.h"
#include "solving_strategies/convergencecriterias/mixed_generic_criteria.h"

// Builder And Solver
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "solving_strategies/builder_and_solvers/explicit_builder.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver_with_lagrange_multiplier.h"

// Linear solvers
#include "linear_solvers/linear_solver.h"

namespace Kratos:: Python
{
    namespace py = pybind11;

    typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

    typedef UblasSpace<std::complex<double>, ComplexCompressedMatrix, boost::numeric::ublas::vector<std::complex<double>>> ComplexSparseSpaceType;
    typedef UblasSpace<std::complex<double>, ComplexMatrix, ComplexVector> ComplexLocalSpaceType;

    //ADDED BY PAOLO (next two)

    double Dot(SparseSpaceType& dummy, SparseSpaceType::VectorType& rX, SparseSpaceType::VectorType& rY)
    {
        return dummy.Dot(rX, rY);
    }

    template< typename TSpaceType >
    void ScaleAndAdd(TSpaceType& dummy, const double A, const typename TSpaceType::VectorType& rX, const double B, typename TSpaceType::VectorType& rY)
    //(const double A,const  VectorType& rX, const double B, VectorType& rY) // rY = (A * rX) + (B * rY)
    {
        dummy.ScaleAndAdd(A, rX, B, rY);
    }

    template< typename TSpaceType >
    void Mult(TSpaceType& dummy, typename TSpaceType::MatrixType& rA, typename TSpaceType::VectorType& rX, typename TSpaceType::VectorType& rY)
    //rY=A*rX (the product is stored inside the rY)
    {
        dummy.Mult(rA, rX, rY);
    }

    void TransposeMult(SparseSpaceType& dummy, SparseSpaceType::MatrixType& rA, SparseSpaceType::VectorType& rX, SparseSpaceType::VectorType& rY)
    //rY=A*rX (the product is stored inside the rY)
    {
        dummy.TransposeMult(rA, rX, rY);
    }

    template< typename TSpaceType >
    typename TSpaceType::IndexType Size(TSpaceType& dummy, typename TSpaceType::VectorType const& rV)
    {
        return rV.size();
    }

    template< typename TSpaceType >
    typename TSpaceType::IndexType Size1(TSpaceType& dummy, typename TSpaceType::MatrixType const& rM)
    {
        return rM.size1();
    }

    template< typename TSpaceType >
    typename TSpaceType::IndexType Size2(TSpaceType& dummy, typename TSpaceType::MatrixType const& rM)
    {
        return rM.size2();
    }

    template< typename TSpaceType >
    void ResizeMatrix(TSpaceType& dummy, typename TSpaceType::MatrixType& A, unsigned int i1, unsigned int i2)
    {
        dummy.Resize(A, i1, i2);
    }

    template< typename TSpaceType >
    void ResizeVector(TSpaceType& dummy, typename TSpaceType::VectorType& x, unsigned int i1)
    {
        dummy.Resize(x, i1);
    }

    template< typename TSpaceType >
    void SetToZeroMatrix(TSpaceType& dummy, typename TSpaceType::MatrixType& A)
    {
        dummy.SetToZero(A);
    }

    template< typename TSpaceType >
    void SetToZeroVector(TSpaceType& dummy, typename TSpaceType::VectorType& x)
    {
        dummy.SetToZero(x);
    }

    template< typename TSpaceType >
    void ClearMatrix(TSpaceType& dummy, typename TSpaceType::MatrixPointerType& A)
    {
        dummy.Clear(A);
    }

    template< typename TSpaceType >
    void ClearVector(TSpaceType& dummy, typename TSpaceType::VectorPointerType& x)
    {
        dummy.Clear(x);
    }

    double TwoNorm(SparseSpaceType& dummy, SparseSpaceType::VectorType& x)
    {
        return dummy.TwoNorm(x);
    }

    template< typename TSpaceType >
    void UnaliasedAdd(TSpaceType& dummy, typename TSpaceType::VectorType& x, const double A, const typename TSpaceType::VectorType& rY) // x+= a*Y
    {
        dummy.UnaliasedAdd(x, A, rY);
    }

    void MoveMesh(Scheme< SparseSpaceType, LocalSpaceType >& dummy, ModelPart::NodesContainerType& rNodes)
    {
        VariableUtils().UpdateCurrentPosition(rNodes, DISPLACEMENT);
    }

    template< typename TSpaceType >
    typename TSpaceType::MatrixPointerType CreateEmptyMatrixPointer(TSpaceType& dummy)
    {
        return dummy.CreateEmptyMatrixPointer();
    }

    template< typename TSpaceType >
    typename TSpaceType::VectorPointerType CreateEmptyVectorPointer(TSpaceType& dummy)
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

    template< typename TSpaceType >
    py::class_< TSpaceType > CreateSpaceInterface(pybind11::module& m, std::string Name)
    {
        py::class_< TSpaceType > binder(m,Name.c_str());
        binder.def(py::init<>());

        binder.def("ClearMatrix", ClearMatrix<TSpaceType>);
        binder.def("ClearVector", ClearVector<TSpaceType>);
        binder.def("ResizeMatrix", ResizeMatrix<TSpaceType>);
        binder.def("ResizeVector", ResizeVector<TSpaceType>);
        binder.def("SetToZeroMatrix", SetToZeroMatrix<TSpaceType>);
        binder.def("SetToZeroVector", SetToZeroVector<TSpaceType>);
        binder.def("ScaleAndAdd", ScaleAndAdd<TSpaceType>);
        //the matrix-vector multiplication
        binder.def("Mult", Mult<TSpaceType>);
        binder.def("Size", Size<TSpaceType>);
        binder.def("Size1", Size1<TSpaceType>);
        binder.def("Size2", Size2<TSpaceType>);
        binder.def("UnaliasedAdd", UnaliasedAdd<TSpaceType>);
        binder.def("CreateEmptyMatrixPointer", CreateEmptyMatrixPointer<TSpaceType>);
        binder.def("CreateEmptyVectorPointer", CreateEmptyVectorPointer<TSpaceType>);

        return binder;
    }

    void AddStrategiesToPython(pybind11::module& m)
    {

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
            .def("Create", &BaseSchemeType::Create)
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
            .def("Check", [](const BaseSchemeType& self, const ModelPart& rModelPart){ return self.Check(rModelPart); })
            .def("GetDefaultParameters",&BaseSchemeType::GetDefaultParameters)
            .def("Info", &BaseSchemeType::Info)
            ;

        py::class_< ResidualBasedIncrementalUpdateStaticScheme< SparseSpaceType, LocalSpaceType>,
            typename ResidualBasedIncrementalUpdateStaticScheme< SparseSpaceType, LocalSpaceType>::Pointer,
            BaseSchemeType >
            (m, "ResidualBasedIncrementalUpdateStaticScheme")
            .def(py::init<Parameters >() )
            .def(py::init< >()
            );

        typedef typename ResidualBasedIncrementalUpdateStaticSchemeSlip< SparseSpaceType, LocalSpaceType>::RotationToolPointerType RotationToolPointerType;

        py::class_< ResidualBasedIncrementalUpdateStaticSchemeSlip< SparseSpaceType, LocalSpaceType>,
            typename ResidualBasedIncrementalUpdateStaticSchemeSlip< SparseSpaceType, LocalSpaceType>::Pointer,
            ResidualBasedIncrementalUpdateStaticScheme< SparseSpaceType, LocalSpaceType> >
            (m,"ResidualBasedIncrementalUpdateStaticSchemeSlip")
            .def(py::init<Parameters >() )
            .def(py::init<unsigned int, unsigned int>())
            .def(py::init<RotationToolPointerType>());

            // Residual Based Bossak Scheme Type
            py::class_< ResidualBasedBossakDisplacementSchemeType,
            typename ResidualBasedBossakDisplacementSchemeType::Pointer,
            BaseSchemeType  >
            (m,"ResidualBasedBossakDisplacementScheme")
            .def(py::init< double >())
            .def(py::init< double, double >())
            ;

            // Residual Based Newmark Scheme Type
            py::class_< ResidualBasedNewmarkDisplacementSchemeType,
            typename ResidualBasedNewmarkDisplacementSchemeType::Pointer,
            BaseSchemeType >(m,"ResidualBasedNewmarkDisplacementScheme")
            .def(py::init< >() )
            .def(py::init<Parameters>() )
            ;

            // Residual Based Pseudo-Static Scheme Type
            py::class_< ResidualBasedPseudoStaticDisplacementSchemeType,
            typename ResidualBasedPseudoStaticDisplacementSchemeType::Pointer,
            BaseSchemeType >(m,"ResidualBasedPseudoStaticDisplacementScheme")
            .def(py::init< const Variable<double>& >() )
            .def(py::init<Parameters>() )
            ;

        // Residual Based BDF displacement Scheme Type
        py::class_< ResidualBasedBDFDisplacementSchemeType,typename ResidualBasedBDFDisplacementSchemeType::Pointer, BaseSchemeType  >(m,"ResidualBasedBDFDisplacementScheme")
            .def(py::init<Parameters >() )
            .def(py::init<  >() )
            .def(py::init <const std::size_t>())
            ;

        // Residual Based BDF custom Scheme Type
        py::class_< ResidualBasedBDFCustomSchemeType, typename ResidualBasedBDFCustomSchemeType::Pointer, BaseSchemeType  >(m,"ResidualBasedBDFCustomScheme")
            .def(py::init<Parameters >() )
            .def(py::init<  >() )
            .def(py::init <const std::size_t>())
            .def(py::init <const std::size_t, Parameters>())
            ;

        // Residual Based Adjoint Static Scheme Type
        typedef ResidualBasedAdjointStaticScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedAdjointStaticSchemeType;
        py::class_<ResidualBasedAdjointStaticSchemeType, typename ResidualBasedAdjointStaticSchemeType::Pointer, BaseSchemeType>
        (m, "ResidualBasedAdjointStaticScheme")
        .def(py::init<AdjointResponseFunction::Pointer>())
        .def("SetResponseFunction", &ResidualBasedAdjointStaticSchemeType::SetResponseFunction, py::arg("new_response_function"))
        ;

        // Residual Based Adjoint Steady Scheme Type
        typedef ResidualBasedAdjointSteadyScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedAdjointSteadySchemeType;
        py::class_<ResidualBasedAdjointSteadySchemeType, typename ResidualBasedAdjointSteadySchemeType::Pointer, ResidualBasedAdjointStaticSchemeType>
        (m, "ResidualBasedAdjointSteadyScheme")
        .def(py::init<AdjointResponseFunction::Pointer>())
        ;

        // Residual Based Adjoint Bossak Scheme Type
        typedef ResidualBasedAdjointBossakScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedAdjointBossakSchemeType;
        py::class_<ResidualBasedAdjointBossakSchemeType, typename ResidualBasedAdjointBossakSchemeType::Pointer, BaseSchemeType>
        (m, "ResidualBasedAdjointBossakScheme")
        .def(py::init<Kratos::Parameters, AdjointResponseFunction::Pointer>())
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
                .def("Create", &ConvergenceCriteriaType::Create)
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
                .def("GetDefaultParameters",&ConvergenceCriteriaType::GetDefaultParameters)
                .def("SetEchoLevel", &ConvergenceCriteriaType::SetEchoLevel)
                .def("Info", &ConvergenceCriteriaType::Info)
                ;

        py::class_< DisplacementCriteria<SparseSpaceType, LocalSpaceType >,
            typename DisplacementCriteria< SparseSpaceType, LocalSpaceType >::Pointer,
            ConvergenceCriteriaType >
            (m,"DisplacementCriteria")
            .def(py::init<Parameters >() )
            .def(py::init< double, double>())
            ;

        py::class_<ResidualCriteria<SparseSpaceType, LocalSpaceType >,
            typename ResidualCriteria< SparseSpaceType, LocalSpaceType >::Pointer,
            ConvergenceCriteriaType >
            (m,"ResidualCriteria")
            .def(py::init<Parameters >() )
            .def(py::init< double, double>())
            ;

        py::class_<And_Criteria<SparseSpaceType, LocalSpaceType >,
            typename And_Criteria< SparseSpaceType, LocalSpaceType >::Pointer,
            ConvergenceCriteriaType >
            (m,"AndCriteria")
            .def(py::init<Parameters >() )
            .def(py::init<ConvergenceCriteriaPointerType, ConvergenceCriteriaPointerType > ())
            ;

        py::class_<Or_Criteria<SparseSpaceType, LocalSpaceType >,
            typename Or_Criteria< SparseSpaceType, LocalSpaceType >::Pointer,
            ConvergenceCriteriaType >
            (m,"OrCriteria")
            .def(py::init<Parameters >() )
            .def(py::init<ConvergenceCriteriaPointerType, ConvergenceCriteriaPointerType > ())
            ;

        typedef typename MixedGenericCriteria<SparseSpaceType, LocalSpaceType>::ConvergenceVariableListType ConvergenceVariableListType;
        py::class_<MixedGenericCriteria<SparseSpaceType, LocalSpaceType >,
            typename MixedGenericCriteria< SparseSpaceType, LocalSpaceType >::Pointer,
            ConvergenceCriteriaType >
            (m,"MixedGenericCriteria")
            .def(py::init< const ConvergenceVariableListType& > ())
            ;

        //********************************************************************
        //********************************************************************
        //********************************************************************
        //********************************************************************
        //********************************************************************
        //********************************************************************
        //Builder and Solver
        typedef BuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > BuilderAndSolverType;
        typedef typename ModelPart::DofsArrayType DofsArrayType;

        py::class_< BuilderAndSolverType, typename BuilderAndSolverType::Pointer>(m,"BuilderAndSolver")
        .def(py::init<LinearSolverType::Pointer > ())
            .def(py::init<LinearSolverType::Pointer, Parameters >() )
            .def(py::init<LinearSolverType::Pointer > ())
            .def("Create", &BuilderAndSolverType::Create)
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
            .def("BuildAndSolveLinearizedOnPreviousIteration", &BuilderAndSolverType::BuildAndSolveLinearizedOnPreviousIteration)
            .def("BuildRHSAndSolve", &BuilderAndSolverType::BuildRHSAndSolve)
            .def("ApplyDirichletConditions", &BuilderAndSolverType::ApplyDirichletConditions)
            .def("ApplyConstraints", &BuilderAndSolverType::ApplyConstraints)
            .def("SetUpDofSet", &BuilderAndSolverType::SetUpDofSet)
            .def("GetDofSet",  [](BuilderAndSolverType& self) -> DofsArrayType& {return self.GetDofSet();}, py::return_value_policy::reference_internal)
            .def("SetUpSystem", &BuilderAndSolverType::SetUpSystem)
            .def("ResizeAndInitializeVectors", &BuilderAndSolverType::ResizeAndInitializeVectors)
            .def("InitializeSolutionStep", &BuilderAndSolverType::InitializeSolutionStep)
            .def("FinalizeSolutionStep", &BuilderAndSolverType::FinalizeSolutionStep)
            .def("CalculateReactions", &BuilderAndSolverType::CalculateReactions)
            .def("Clear", &BuilderAndSolverType::Clear)
            .def("Check", &BuilderAndSolverType::Check)
            .def("GetDefaultParameters",&BuilderAndSolverType::GetDefaultParameters)
            .def("SetEchoLevel", &BuilderAndSolverType::SetEchoLevel)
            .def("GetEchoLevel", &BuilderAndSolverType::GetEchoLevel)
            .def("Info", &BuilderAndSolverType::Info)
            ;

        // Explicit builder
        typedef ExplicitBuilder< SparseSpaceType, LocalSpaceType > ExplicitBuilderType;

        py::class_<ExplicitBuilderType, typename ExplicitBuilderType::Pointer>(m, "ExplicitBuilder")
            .def(py::init<>())
            .def("SetCalculateReactionsFlag", &ExplicitBuilderType::SetCalculateReactionsFlag)
            .def("GetCalculateReactionsFlag", &ExplicitBuilderType::GetCalculateReactionsFlag)
            .def("SetDofSetIsInitializedFlag", &ExplicitBuilderType::SetDofSetIsInitializedFlag)
            .def("GetDofSetIsInitializedFlag", &ExplicitBuilderType::GetDofSetIsInitializedFlag)
            .def("SetResetDofSetFlag", &ExplicitBuilderType::SetResetDofSetFlag)
            .def("GetResetDofSetFlag", &ExplicitBuilderType::GetResetDofSetFlag)
            .def("SetResetLumpedMassVectorFlag", &ExplicitBuilderType::SetResetLumpedMassVectorFlag)
            .def("GetResetLumpedMassVectorFlag", &ExplicitBuilderType::GetResetLumpedMassVectorFlag)
            .def("GetEquationSystemSize", &ExplicitBuilderType::GetEquationSystemSize)
            .def("BuildRHS", &ExplicitBuilderType::BuildRHS)
            .def("ApplyConstraints", &ExplicitBuilderType::ApplyConstraints)
            .def("GetDofSet", [](ExplicitBuilderType& self) -> DofsArrayType& {return self.GetDofSet();}, py::return_value_policy::reference_internal)
            .def("GetLumpedMassMatrixVector", &ExplicitBuilderType::GetLumpedMassMatrixVector, py::return_value_policy::reference_internal)
            .def("Initialize", &ExplicitBuilderType::Initialize)
            .def("InitializeSolutionStep", &ExplicitBuilderType::InitializeSolutionStep)
            .def("FinalizeSolutionStep", &ExplicitBuilderType::FinalizeSolutionStep)
            .def("Clear", &ExplicitBuilderType::Clear)
            .def("Check", &ExplicitBuilderType::Check)
            .def("GetDefaultParameters",&ExplicitBuilderType::GetDefaultParameters)
            .def("SetEchoLevel", &ExplicitBuilderType::SetEchoLevel)
            .def("GetEchoLevel", &ExplicitBuilderType::GetEchoLevel)
            .def("Info", &ExplicitBuilderType::Info)
            ;

        typedef ResidualBasedEliminationBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedEliminationBuilderAndSolverType;
        py::class_< ResidualBasedEliminationBuilderAndSolverType, ResidualBasedEliminationBuilderAndSolverType::Pointer, BuilderAndSolverType>(m,"ResidualBasedEliminationBuilderAndSolver")
        .def(py::init<LinearSolverType::Pointer, Parameters >() )
        .def(py::init< LinearSolverType::Pointer > ())
        ;

        typedef ResidualBasedEliminationBuilderAndSolverWithConstraints< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedEliminationBuilderAndSolverWithConstraintsType;
        py::class_< ResidualBasedEliminationBuilderAndSolverWithConstraintsType, ResidualBasedEliminationBuilderAndSolverWithConstraintsType::Pointer, BuilderAndSolverType>(m,"ResidualBasedEliminationBuilderAndSolverWithConstraints")
        .def(py::init< LinearSolverType::Pointer > ())
        .def(py::init< LinearSolverType::Pointer, bool > ())
        .def(py::init< LinearSolverType::Pointer, bool, bool > ())
        .def(py::init< LinearSolverType::Pointer, Parameters > ())
        ;

        typedef ResidualBasedBlockBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedBlockBuilderAndSolverType;
        py::class_< ResidualBasedBlockBuilderAndSolverType, ResidualBasedBlockBuilderAndSolverType::Pointer,BuilderAndSolverType>(m,"ResidualBasedBlockBuilderAndSolver")
        .def(py::init< LinearSolverType::Pointer > ())
        .def(py::init< LinearSolverType::Pointer, Parameters > ())
        ;

        typedef ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplierType;
        py::class_< ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplierType, ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplierType::Pointer,BuilderAndSolverType>(m,"ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier")
        .def(py::init< LinearSolverType::Pointer > ())
        .def(py::init< LinearSolverType::Pointer, Parameters > ())
        ;

        //********************************************************************
        //********************************************************************
        //********************************************************************
        //********************************************************************
        //********************************************************************
        //********************************************************************

        auto sparse_space_binder = CreateSpaceInterface< SparseSpaceType >(m,"UblasSparseSpace");
        sparse_space_binder.def("TwoNorm", TwoNorm);
        //the dot product of two vectors
        sparse_space_binder.def("Dot", Dot);
        sparse_space_binder.def("TransposeMult", TransposeMult);

        auto cplx_sparse_space_binder = CreateSpaceInterface< ComplexSparseSpaceType >(m,"UblasComplexSparseSpace");

        //********************************************************************
        //********************************************************************
        //********************************************************************
        //********************************************************************
        //********************************************************************
        //********************************************************************

        //strategy base class
        typedef SolvingStrategy< SparseSpaceType, LocalSpaceType > BaseSolvingStrategyType;
        py::class_< BaseSolvingStrategyType, typename BaseSolvingStrategyType::Pointer >(m,"BaseSolvingStrategy")
            .def(py::init<ModelPart&, Parameters >() )
            .def(py::init < ModelPart&, bool >())
            .def("Create", &BaseSolvingStrategyType::Create)
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
            .def("MoveMeshFlag", &BaseSolvingStrategyType::GetMoveMeshFlag)
            .def("GetMoveMeshFlag", &BaseSolvingStrategyType::GetMoveMeshFlag)
            .def("MoveMesh", &BaseSolvingStrategyType::MoveMesh)
            .def("Clear", &BaseSolvingStrategyType::Clear)
            .def("Check", &BaseSolvingStrategyType::Check)
            .def("GetDefaultParameters",&BaseSolvingStrategyType::GetDefaultParameters)
            .def("GetSystemMatrix", &BaseSolvingStrategyType::GetSystemMatrix, py::return_value_policy::reference_internal)
            .def("GetSystemVector", &BaseSolvingStrategyType::GetSystemVector, py::return_value_policy::reference_internal)
            .def("GetSolutionVector", &BaseSolvingStrategyType::GetSolutionVector, py::return_value_policy::reference_internal)
            .def("InitializeSolutionStep", &BaseSolvingStrategyType::InitializeSolutionStep)
            .def("FinalizeSolutionStep", &BaseSolvingStrategyType::FinalizeSolutionStep)
            .def("SolveSolutionStep", &BaseSolvingStrategyType::SolveSolutionStep)
            .def("GetModelPart", [](BaseSolvingStrategyType& self) -> ModelPart& { return self.GetModelPart(); })
            .def("Info", &BaseSolvingStrategyType::Info)
            ;

        //implicit strategy base class
        typedef ImplicitSolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > ImplicitSolvingStrategyType;
        py::class_<ImplicitSolvingStrategyType, typename ImplicitSolvingStrategyType::Pointer, BaseSolvingStrategyType>(m,"ImplicitSolvingStrategy")
            .def(py::init<ModelPart&, Parameters >() )
            .def(py::init < ModelPart&, bool >())
            .def("SetStiffnessMatrixIsBuilt", &ImplicitSolvingStrategyType::SetStiffnessMatrixIsBuilt)
            .def("GetStiffnessMatrixIsBuilt", &ImplicitSolvingStrategyType::GetStiffnessMatrixIsBuilt)
            ;

        typedef ExplicitSolvingStrategy< SparseSpaceType, LocalSpaceType > BaseExplicitSolvingStrategyType;
        py::class_<BaseExplicitSolvingStrategyType, typename BaseExplicitSolvingStrategyType::Pointer, BaseSolvingStrategyType>(m, "ExplicitSolvingStrategy")
            .def(py::init<ModelPart &, bool, int>())
            .def(py::init<ModelPart&, typename ExplicitBuilderType::Pointer, bool, int>())
            ;

        typedef ExplicitSolvingStrategyRungeKutta4< SparseSpaceType, LocalSpaceType > ExplicitSolvingStrategyRungeKutta4Type;
        py::class_<ExplicitSolvingStrategyRungeKutta4Type, typename ExplicitSolvingStrategyRungeKutta4Type::Pointer, BaseExplicitSolvingStrategyType>(m, "ExplicitSolvingStrategyRungeKutta4")
            .def(py::init<ModelPart&, bool, int>())
            .def(py::init<ModelPart&, Parameters>())
            .def(py::init<ModelPart&, typename ExplicitBuilderType::Pointer, bool, int>())
            ;

        typedef ExplicitSolvingStrategyRungeKutta3TVD< SparseSpaceType, LocalSpaceType > ExplicitSolvingStrategyRungeKutta3TVDType;
        py::class_<ExplicitSolvingStrategyRungeKutta3TVDType, typename ExplicitSolvingStrategyRungeKutta3TVDType::Pointer, BaseExplicitSolvingStrategyType>(m, "ExplicitSolvingStrategyRungeKutta3TVD")
            .def(py::init<ModelPart&, bool, int>())
            .def(py::init<ModelPart&, Parameters>())
            .def(py::init<ModelPart&, typename ExplicitBuilderType::Pointer, bool, int>())
            ;

        typedef ExplicitSolvingStrategyRungeKutta1< SparseSpaceType, LocalSpaceType > ExplicitSolvingStrategyRungeKutta1Type;
        py::class_<ExplicitSolvingStrategyRungeKutta1Type, typename ExplicitSolvingStrategyRungeKutta1Type::Pointer, BaseExplicitSolvingStrategyType>(m, "ExplicitSolvingStrategyRungeKutta1")
            .def(py::init<ModelPart&, bool, int>())
            .def(py::init<ModelPart&, Parameters>())
            .def(py::init<ModelPart&, typename ExplicitBuilderType::Pointer, bool, int>())
            ;

        typedef ExplicitSolvingStrategyBFECC< SparseSpaceType, LocalSpaceType > ExplicitSolvingStrategyBFECCType;
        py::class_<ExplicitSolvingStrategyBFECCType, typename ExplicitSolvingStrategyBFECCType::Pointer, BaseExplicitSolvingStrategyType>(m, "ExplicitSolvingStrategyBFECC")
            .def(py::init<ModelPart&, bool, int>())
            .def(py::init<ModelPart&, Parameters>())
            .def(py::init<ModelPart&, typename ExplicitBuilderType::Pointer, bool, int>())
            ;

        typedef ResidualBasedLinearStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedLinearStrategyType;

        py::class_< ResidualBasedLinearStrategyType, typename ResidualBasedLinearStrategyType::Pointer,ImplicitSolvingStrategyType >
            (m,"ResidualBasedLinearStrategy")
            .def(py::init<ModelPart&, Parameters >() )
            .def(py::init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, bool, bool, bool, bool >())
            .def(py::init < ModelPart& ,  BaseSchemeType::Pointer, BuilderAndSolverType::Pointer, bool, bool, bool,  bool  >())
            .def(py::init([](ModelPart& rModelPart, BaseSchemeType::Pointer pScheme, LinearSolverType::Pointer pLinearSolver, BuilderAndSolverType::Pointer pBuilderAndSolver, bool CalculateReactionFlag, bool ReformDofSetAtEachStep, bool CalculateNormDxFlag, bool MoveMeshFlag) {
                KRATOS_WARNING("ResidualBasedLinearStrategy") << "Using deprecated constructor. Please use constructor without linear solver.";
                return std::shared_ptr<ResidualBasedLinearStrategyType>(new ResidualBasedLinearStrategyType(rModelPart, pScheme, pBuilderAndSolver, CalculateReactionFlag, ReformDofSetAtEachStep, CalculateNormDxFlag, MoveMeshFlag));
            }))
            .def("GetScheme", &ResidualBasedLinearStrategyType::GetScheme)
            .def("GetResidualNorm", &ResidualBasedLinearStrategyType::GetResidualNorm)
            .def("SetBuilderAndSolver", &ResidualBasedLinearStrategyType::SetBuilderAndSolver)
            .def("SetReformDofSetAtEachStepFlag", &ResidualBasedLinearStrategyType::SetReformDofSetAtEachStepFlag)
            .def("GetReformDofSetAtEachStepFlag", &ResidualBasedLinearStrategyType::GetReformDofSetAtEachStepFlag)
            ;

        typedef ResidualBasedNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedNewtonRaphsonStrategyType;

        py::class_< ResidualBasedNewtonRaphsonStrategyType, typename ResidualBasedNewtonRaphsonStrategyType::Pointer, ImplicitSolvingStrategyType >
            (m,"ResidualBasedNewtonRaphsonStrategy")
            .def(py::init<ModelPart&, Parameters >() )
            .def(py::init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, int, bool, bool, bool >())
            .def(py::init < ModelPart&, BaseSchemeType::Pointer, ConvergenceCriteriaType::Pointer, BuilderAndSolverType::Pointer, int, bool, bool, bool >())
            .def(py::init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, Parameters>())
            .def(py::init < ModelPart&, BaseSchemeType::Pointer, ConvergenceCriteriaType::Pointer, BuilderAndSolverType::Pointer, Parameters>())
            .def(py::init([](ModelPart& rModelPart, BaseSchemeType::Pointer pScheme, LinearSolverType::Pointer pLinearSolver, ConvergenceCriteriaType::Pointer pConvergenceCriteria, BuilderAndSolverType::Pointer pBuilderAndSolver, int MaxIterations, bool CalculateReactions, bool ReformDofSetAtEachStep, bool MoveMeshFlag) {
                    KRATOS_WARNING("ResidualBasedNewtonRaphsonStrategy") << "Using deprecated constructor. Please use constructor without linear solver.";
                    return std::shared_ptr<ResidualBasedNewtonRaphsonStrategyType>(new ResidualBasedNewtonRaphsonStrategyType(rModelPart, pScheme, pConvergenceCriteria, pBuilderAndSolver, MaxIterations, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag));
                }))
            .def(py::init([](ModelPart& rModelPart, BaseSchemeType::Pointer pScheme, LinearSolverType::Pointer pLinearSolver, ConvergenceCriteriaType::Pointer pConvergenceCriteria, BuilderAndSolverType::Pointer pBuilderAndSolver, Parameters Settings) {
                    KRATOS_WARNING("ResidualBasedNewtonRaphsonStrategy") << "Using deprecated constructor. Please use constructor without linear solver.";
                    return std::shared_ptr<ResidualBasedNewtonRaphsonStrategyType>(new ResidualBasedNewtonRaphsonStrategyType(rModelPart, pScheme, pConvergenceCriteria, pBuilderAndSolver, Settings));
                }))
            .def("SetMaxIterationNumber", &ResidualBasedNewtonRaphsonStrategyType::SetMaxIterationNumber)
            .def("GetMaxIterationNumber", &ResidualBasedNewtonRaphsonStrategyType::GetMaxIterationNumber)
            .def("SetKeepSystemConstantDuringIterations", &ResidualBasedNewtonRaphsonStrategyType::SetKeepSystemConstantDuringIterations)
            .def("GetKeepSystemConstantDuringIterations", &ResidualBasedNewtonRaphsonStrategyType::GetKeepSystemConstantDuringIterations)
            .def("SetInitializePerformedFlag", &ResidualBasedNewtonRaphsonStrategyType::SetInitializePerformedFlag)
            .def("GetInitializePerformedFlag", &ResidualBasedNewtonRaphsonStrategyType::GetInitializePerformedFlag)
            .def("SetUseOldStiffnessInFirstIterationFlag", &ResidualBasedNewtonRaphsonStrategyType::SetUseOldStiffnessInFirstIterationFlag)
            .def("GetUseOldStiffnessInFirstIterationFlag", &ResidualBasedNewtonRaphsonStrategyType::GetUseOldStiffnessInFirstIterationFlag)
            .def("SetReformDofSetAtEachStepFlag", &ResidualBasedNewtonRaphsonStrategyType::SetReformDofSetAtEachStepFlag)
            .def("GetReformDofSetAtEachStepFlag", &ResidualBasedNewtonRaphsonStrategyType::GetReformDofSetAtEachStepFlag)
            ;

        // ARC-LENGTH
        typedef ArcLengthStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > ArcLengthStrategyStrategyType;
        py::class_< ArcLengthStrategyStrategyType, typename ArcLengthStrategyStrategyType::Pointer, ResidualBasedNewtonRaphsonStrategyType >
            (m,"ArcLengthStrategy")
            .def(py::init < ModelPart&, BaseSchemeType::Pointer, ConvergenceCriteriaType::Pointer, BuilderAndSolverType::Pointer, Parameters >())
            ;

        py::class_< AdaptiveResidualBasedNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >,
            typename AdaptiveResidualBasedNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::Pointer,
            ImplicitSolvingStrategyType >
            (m,"AdaptiveResidualBasedNewtonRaphsonStrategy")
            .def(py::init<ModelPart&, Parameters >() )
            .def(py::init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, int, int, bool, bool, bool, double, double, int
            >())
            ;

        typedef LineSearchStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > LineSearchStrategyType;
        py::class_< LineSearchStrategyType, typename LineSearchStrategyType::Pointer, ResidualBasedNewtonRaphsonStrategyType >
            (m,"LineSearchStrategy")
            .def(py::init<ModelPart&, Parameters >() )
            .def(py::init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, int, bool, bool, bool >())
            .def(py::init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, Parameters >())
            .def(py::init < ModelPart&, BaseSchemeType::Pointer, ConvergenceCriteriaType::Pointer, BuilderAndSolverType::Pointer, int, bool, bool, bool >())
            .def(py::init([](ModelPart& rModelPart, BaseSchemeType::Pointer pScheme, LinearSolverType::Pointer pLinearSolver, ConvergenceCriteriaType::Pointer pConvergenceCriteria, BuilderAndSolverType::Pointer pBuilderAndSolver, int MaxIterations, bool CalculateReactions, bool ReformDofSetAtEachStep, bool MoveMeshFlag) {
                    KRATOS_WARNING("LineSearchStrategy") << "Using deprecated constructor. Please use constructor without linear solver.";
                    return std::shared_ptr<LineSearchStrategyType>(new LineSearchStrategyType(rModelPart, pScheme, pConvergenceCriteria, pBuilderAndSolver, MaxIterations, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag));
                }))
            .def(py::init < ModelPart&, BaseSchemeType::Pointer, ConvergenceCriteriaType::Pointer, BuilderAndSolverType::Pointer, Parameters >())
            ;

        py::class_<SensitivityBuilderScheme, typename SensitivityBuilderScheme::Pointer>
        (m,"SensitivityBuilderScheme")
            .def(py::init< >())
            .def("Initialize", &SensitivityBuilderScheme::Initialize)
            .def("InitializeSolutionStep", &SensitivityBuilderScheme::InitializeSolutionStep)
            .def("FinalizeSolutionStep", &SensitivityBuilderScheme::FinalizeSolutionStep)
            .def("Finalize", &SensitivityBuilderScheme::Finalize)
            .def("Update", &SensitivityBuilderScheme::Update)
            .def("Clear", &SensitivityBuilderScheme::Clear)
            .def("Check", &SensitivityBuilderScheme::Check)
            .def("Info", &SensitivityBuilderScheme::Info)
            ;

    }

} // namespace Kratos::Python.
