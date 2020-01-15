//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#if !defined(MOR_SECOND_ORDER_KRYLOV_STRATEGY)
#define MOR_SECOND_ORDER_KRYLOV_STRATEGY

// System includes

// External includes

// Project includes
#include "includes/define.h"
// #include "solving_strategies/strategies/solving_strategy.h"
// #include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "utilities/builtin_timer.h"
#include "utilities/qr_utility.h"
#include "custom_strategies/custom_strategies/mor_offline_second_order_strategy.hpp"
#include "custom_utilities/eigen_qr_utility.hpp"

//default builder and solver
#include "custom_strategies/custom_builder_and_solvers/system_matrix_builder_and_solver.hpp"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}

///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class MorSecondOrderKrylovStrategy
 * @ingroup KratosCore
 * @brief This is the linear MOR matrix output strategy
 * @details This strategy builds the K and M matrices and outputs them
 * @author Aditya Ghantasala
 */
template <class TSparseSpace,
          class TDenseSpace,  // = DenseSpace<double>,
          class TLinearSolver, //= LinearSolver<TSparseSpace,TDenseSpace>
          class TReducedSparseSpace,
          class TReducedDenseSpace
          >
class MorSecondOrderKrylovStrategy
    // : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
    : public MorOfflineSecondOrderStrategy< TSparseSpace, TDenseSpace, TLinearSolver, TReducedSparseSpace, TReducedDenseSpace >
{
  public:
    ///@name Type Definitions
    ///@{
    // Counted pointer of ClassName
    KRATOS_CLASS_POINTER_DEFINITION(MorSecondOrderKrylovStrategy);

    // typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
    typedef MorOfflineSecondOrderStrategy< TSparseSpace, TDenseSpace, TLinearSolver, TReducedSparseSpace, TReducedDenseSpace > BaseType;

    typedef SystemMatrixBuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver > TBuilderAndSolverType;

    typedef typename BaseType::TDataType TDataType;

    typedef TSparseSpace SparseSpaceType;

    typedef TDenseSpace DenseSpaceType;

    typedef typename TDenseSpace::MatrixType TDenseMatrixType;

    typedef typename TDenseSpace::MatrixPointerType TDenseMatrixPointerType;

    typedef typename BaseType::TSchemeType TSchemeType;

    //typedef typename BaseType::DofSetType DofSetType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;

    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;

    typedef typename TReducedSparseSpace::MatrixType TReducedSparseMatrixType;

    typedef typename TReducedSparseSpace::VectorType TReducedSparseVectorType;

    typedef typename TReducedDenseSpace::MatrixType TReducedDenseMatrixType;


    ///@}
    ///@name Life Cycle

    ///@{

    /**
     * Default constructor for the damped case
     * @param rModelPart The model part of the problem
     * @param pScheme The integration schemed
     * @param MoveMeshFlag The flag that allows to move the mesh
     */
    MorSecondOrderKrylovStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TBuilderAndSolverType::Pointer pBuilderAndSolver,
        typename LinearSolver< TReducedSparseSpace, TReducedDenseSpace >::Pointer pNewLinearSolver,
        vector< double > samplingPoints,
        bool MoveMeshFlag = false)
        : BaseType(rModelPart, pScheme, pBuilderAndSolver, pNewLinearSolver, MoveMeshFlag)
    {
        KRATOS_TRY;

        // Saving the scheme
        // this->SetScheme(pScheme);

        // this->mpBuilderAndSolver = pBuilderAndSolver;
        // // Setting up the default builder and solver
        // // this->SetBuilderAndSolver(typename TBuilderAndSolverType::Pointer(
        // //     new TBuilderAndSolverType(pNewLinearSolver)));

        // // Saving the linear solver
        // this->SetLinearSolver(pNewLinearSolver);

        // // Set flags to start correcty the calculations
        // mSolutionStepIsInitialized = false;
        // mInitializeWasPerformed = false;

        // // Tells to the builder and solver if the reactions have to be Calculated or not
        // this->GetBuilderAndSolver()->SetCalculateReactionsFlag(false);

        // // Tells to the Builder And Solver if the system matrix and vectors need to
        // // be reshaped at each step or not
        // this->GetBuilderAndSolver()->SetReshapeMatrixFlag(false);

        // // Set EchoLevel to the default value (only time is displayed)
        // this->SetEchoLevel(1);

        // // By default the matrices are rebuilt at each iteration
        // this->SetRebuildLevel(0);

        std::cout << "subclass constructor\n";
        // KRATOS_WATCH(mSolutionStepIsInitialized)
        KRATOS_WATCH(this->mSolutionStepIsInitialized)

        this->mUseDamping = true;

        mSamplingPoints = samplingPoints;

        KRATOS_CATCH("");
    }

    /**
     * @brief Destructor.
     * @details In trilinos third party library, the linear solver's preconditioner should be freed before the system matrix. We control the deallocation order with Clear().
     */
    ~MorSecondOrderKrylovStrategy() override
    {
        this->Clear();
    }


    /**
     * @brief Performs all the required operations that should be done (for each step) before solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     */
    virtual void InitializeSolutionStep() override
    {
        KRATOS_TRY;
        std::cout << "initialize????\n";
        KRATOS_WATCH(this->mSolutionStepIsInitialized)

        if (this->mSolutionStepIsInitialized == false)
        {
            BaseType::InitializeSolutionStep();
            std::cout << "subclass action!\n";
            this->mSolutionStepIsInitialized = true;
        }

        KRATOS_CATCH("");

    }


    //*********************************************************************************
    /**OPERATIONS ACCESSIBLE FROM THE INPUT: **/

    /**
     * @brief Solves the current step. This function returns true if a solution has been found, false otherwise.
     */
    bool SolveSolutionStep() override
    {
        KRATOS_TRY;
        std::cout << "hello! this is where the second order krylov MOR magic happens" << std::endl;
        typename TSchemeType::Pointer p_scheme = this->GetScheme();
        typename TBuilderAndSolverType::Pointer p_builder_and_solver = this->GetBuilderAndSolver();

        TSystemMatrixType& r_K = this->GetSystemMatrix();
        TSystemMatrixType& r_M = this->GetMassMatrix();
        TSystemMatrixType& r_D = this->GetDampingMatrix();
        TSystemVectorType& r_RHS = this->GetSystemVector();

        TSystemVectorType tmp(r_K.size1(), 0.0);

        p_builder_and_solver->BuildRHS(p_scheme, BaseType::GetModelPart(), r_RHS);

        p_builder_and_solver->BuildStiffnessMatrix(p_scheme, BaseType::GetModelPart(), r_K, tmp);
        p_builder_and_solver->ApplyDirichletConditions(p_scheme, BaseType::GetModelPart(), r_K, tmp, r_RHS);

        p_builder_and_solver->BuildMassMatrix(p_scheme, BaseType::GetModelPart(), r_M, tmp);
        p_builder_and_solver->ApplyDirichletConditionsForMassMatrix(p_scheme, BaseType::GetModelPart(), r_M);

        p_builder_and_solver->BuildDampingMatrix(p_scheme, BaseType::GetModelPart(), r_D, tmp);
        p_builder_and_solver->ApplyDirichletConditionsForDampingMatrix(p_scheme, BaseType::GetModelPart(), r_D);

        // EchoInfo(0);
        const unsigned int system_size = p_builder_and_solver->GetEquationSystemSize();
        //sampling points
        KRATOS_WATCH(mSamplingPoints)
        const std::size_t n_sampling_points = mSamplingPoints.size();
        const std::size_t reduced_system_size = 3 * n_sampling_points;

        //initialize sb, As, AAs vectors
        auto s = TReducedSparseSpace::CreateEmptyVectorPointer();
        auto& rs = *s;
        TReducedSparseSpace::Resize(rs,system_size);
        // TReducedSparseSpace::Set(rs,0.0);
        auto As = TReducedSparseSpace::CreateEmptyVectorPointer();
        auto& rAs = *As;
        TReducedSparseSpace::Resize(rAs,system_size);
        // TReducedSparseSpace::Set(rAs,0.0);
        auto AAs = TReducedSparseSpace::CreateEmptyVectorPointer();
        auto& rAAs = *AAs;
        TReducedSparseSpace::Resize(rAAs,system_size);
        // TReducedSparseSpace::Set(rAAs,0.0);

        auto kdyn = TReducedSparseSpace::CreateEmptyMatrixPointer();
        auto& r_kdyn = *kdyn;
        TReducedSparseSpace::Resize(r_kdyn, system_size, system_size);

        // auto tmp_basis = TReducedDenseSpace::CreateEmptyMatrixPointer();
        // auto& r_tmp_basis = *tmp_basis;
        // TReducedDenseSpace::Resize(r_tmp_basis, system_size, reduced_system_size);

        auto& r_basis = this->GetBasis();
        // SparseSpaceType::Resize(r_basis, system_size, reduced_system_size);
        r_basis.resize(system_size,reduced_system_size,false);


        TReducedSparseVectorType aux(r_K.size1(), std::complex<double>(0.0,0.0));
        TReducedSparseMatrixType tmp_C; 
        tmp_C = TReducedSparseMatrixType(r_D);
        TReducedSparseVectorType tmp_rhs;
        tmp_rhs = TReducedSparseVectorType(r_RHS);
        // LocalSystemMatrixType& r_Kr = *mpKr;
        // LocalSystemMatrixType& r_Dr = *mpDr;
        // LocalSystemMatrixType& r_Mr = *mpMr;
        // ComplexDenseVectorType& r_rhs = *mpRHSr;

        for( size_t i = 0; i < n_sampling_points; ++i )
        {
            KRATOS_WATCH( mSamplingPoints(i) )
            std::complex<double> c_sampling_point(0.0, mSamplingPoints(i));
            // r_kdyn = r_K - ( std::pow( mSamplingPoints(i), 2.0 ) * r_M );    // Without Damping
            r_kdyn = r_K - ( std::pow( mSamplingPoints(i), 2.0 ) * r_M ) + c_sampling_point*tmp_C; // With Damping

            this->mpLinearSolver->Solve( r_kdyn, rs, tmp_rhs );
            // p_builder_and_solver->GetLinearSystemSolver()->Solve( r_kdyn, rs, &tmp_rhs );
            aux = prod( r_M, rs );
        
            // p_builder_and_solver->GetLinearSystemSolver()->Solve( r_kdyn, rAs, aux );
            this->mpLinearSolver->Solve( r_kdyn, rAs, aux );
            aux = prod( r_M, rAs );

            // p_builder_and_solver->GetLinearSystemSolver()->Solve( r_kdyn, rAAs, aux );
            this->mpLinearSolver->Solve( r_kdyn, rAAs, aux );

            column( r_basis, (i*3) ) = rs;
            column( r_basis, (i*3)+1 ) = rAs;
            column( r_basis, (i*3)+2 ) = rAAs;

            // std::stringstream matrix_market_name;
            // matrix_market_name << "kdyn_complex.mm";
            // TSparseSpace::WriteMatrixMarketMatrix((char *)(matrix_market_name.str()).c_str(), r_kdyn, false);
            // std::stringstream matrix_market_name2;
            // matrix_market_name2 << "M_real.mm";
            // TSparseSpace::WriteMatrixMarketMatrix((char *)(matrix_market_name2.str()).c_str(), r_M, false);
        }

        // KRATOS_WATCH(r_tmp_basis)
        // r_basis = r_tmp_basis;
        KRATOS_WATCH(r_basis(1,2))

        EigenQrUtility<TReducedDenseSpace> bla;// = new EigenQrUtility();
        bla.MatrixQ(r_basis);
        KRATOS_WATCH(r_basis(1,2))

        //orthogonalize the basis -> basis_r
        // mQR_decomposition.compute( system_size, 3*n_sampling_points, &(r_basis)(0,0) );
        // mQR_decomposition.compute( system_size, 3*n_sampling_points, &(std::real(r_tmp_basis))(0,0) );
        // mQR_decomposition.compute_q();

        // for( size_t i = 0; i < system_size; ++i )
        // {
        //     for( size_t j = 0; j < (3*n_sampling_points); ++j )
        //     {
        //         r_basis(i,j) = mQR_decomposition.Q(i,j);

        //     }
        // }

        // project the system matrices onto the Krylov subspace
        KRATOS_WATCH(this->GetRHSr())
        auto& r_force_vector_reduced = this->GetRHSr();
        auto& r_stiffness_matrix_reduced = this->GetKr();
        auto& r_mass_matrix_reduced = this->GetMr();
        auto& r_damping_matrix_reduced = this->GetDr();
        r_force_vector_reduced.resize( reduced_system_size, false);

        r_force_vector_reduced = prod( r_RHS, r_basis );

        TReducedDenseMatrixType T = prod( herm( r_basis ), r_K );
        r_stiffness_matrix_reduced = prod( T, r_basis );

        T = prod( herm( r_basis ), r_M );
        r_mass_matrix_reduced = prod( T, r_basis );

        T = prod( herm( r_basis ), r_D );
        r_damping_matrix_reduced = prod( T, r_basis );

        KRATOS_WATCH(r_force_vector_reduced)
        KRATOS_WATCH(r_stiffness_matrix_reduced)
        KRATOS_WATCH(r_mass_matrix_reduced)
        KRATOS_WATCH(r_damping_matrix_reduced)

        std::cout << "MOR offline solve finished" << std::endl;

		return true;

        KRATOS_CATCH("");
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const override
    {
        return "MorSecondOrderKrylovStrategy";
    }

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}

  protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

  private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    bool mUseComplexFlag;
    vector< double > mSamplingPoints;
    QR<double, row_major> mQR_decomposition;

    /**
     * @brief Flag telling if it is needed to reform the DofSet at each
    solution step or if it is possible to form it just once
    * @details Default = false
        - true  : Reforme at each time step
        - false : Form just one (more efficient)
     */
    // bool mReformDofSetAtEachStep;

    // bool mSolutionStepIsInitialized; /// Flag to set as initialized the solution step

    // bool mInitializeWasPerformed; /// Flag to set as initialized the strategy

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /**
     * Copy constructor.
     */

    MorSecondOrderKrylovStrategy(const MorSecondOrderKrylovStrategy &Other){};

    ///@}

}; /* Class MorSecondOrderKrylovStrategy */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos. */

#endif /* MOR_SECOND_ORDER_KRYLOV_STRATEGY  defined */