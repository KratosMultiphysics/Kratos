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
#include "utilities/builtin_timer.h"
#include "custom_strategies/custom_strategies/mor_offline_second_order_strategy.hpp"
#include "custom_utilities/orthogonalization_utility.hpp"

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

    // typedef SystemMatrixBuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver > TBuilderAndSolverType;

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
        typename BaseType::TBuilderAndSolverType::Pointer pBuilderAndSolver,
        typename LinearSolver< TReducedSparseSpace, TReducedDenseSpace >::Pointer pNewLinearSolver,
        vector< double > samplingPoints,
        bool MoveMeshFlag = false)
        : BaseType(rModelPart, pScheme, pBuilderAndSolver, pNewLinearSolver, true, MoveMeshFlag),
            mSamplingPoints(samplingPoints)
    {
        KRATOS_TRY;

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

        if (this->mSolutionStepIsInitialized == false)
        {
            BaseType::InitializeSolutionStep();

            const std::size_t reduced_system_size = 3 * mSamplingPoints.size();
            const unsigned int system_size = this->GetBuilderAndSolver()->GetEquationSystemSize();

            TReducedDenseSpace::Resize(this->GetKr(), reduced_system_size, reduced_system_size);
            TReducedDenseSpace::Resize(this->GetDr(), reduced_system_size, reduced_system_size);
            TReducedDenseSpace::Resize(this->GetMr(), reduced_system_size, reduced_system_size);
            TReducedDenseSpace::Resize(this->GetBasis(), system_size, reduced_system_size);


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
        typename BaseType::TBuilderAndSolverType::Pointer p_builder_and_solver = this->GetBuilderAndSolver();
        const int rank = BaseType::GetModelPart().GetCommunicator().MyPID();

        TSystemMatrixType& r_K = this->GetSystemMatrix();
        TSystemMatrixType& r_M = this->GetMassMatrix();
        TSystemMatrixType& r_D = this->GetDampingMatrix();
        TSystemVectorType& r_RHS = this->GetSystemVector();
        TReducedDenseMatrixType& r_basis = this->GetBasis();

        const size_t system_size = p_builder_and_solver->GetEquationSystemSize();
        const size_t n_sampling_points = mSamplingPoints.size();

        //initialize s, As, AAs vectors
        auto s = TReducedSparseSpace::CreateEmptyVectorPointer();
        auto& rs = *s;
        TReducedSparseSpace::Resize(rs,system_size);
        auto As = TReducedSparseSpace::CreateEmptyVectorPointer();
        auto& rAs = *As;
        TReducedSparseSpace::Resize(rAs,system_size);
        auto AAs = TReducedSparseSpace::CreateEmptyVectorPointer();
        auto& rAAs = *AAs;
        TReducedSparseSpace::Resize(rAAs,system_size);

        auto kdyn = TReducedSparseSpace::CreateEmptyMatrixPointer();
        auto& r_kdyn = *kdyn;
        TReducedSparseSpace::Resize(r_kdyn, system_size, system_size);


        TReducedSparseVectorType aux(system_size, std::complex<double>(0.0,0.0));
        TReducedSparseVectorType tmp_rhs;
        tmp_rhs = TReducedSparseVectorType(r_RHS);
        
        BuiltinTimer basis_construction_time;

        for( size_t i = 0; i < n_sampling_points; ++i )
        {
            // build dynamic stiffness matrix
            r_kdyn = r_D;
            r_kdyn *= std::complex<double>(0.0, mSamplingPoints(i));
            r_kdyn += r_K;
            r_kdyn -= std::pow(mSamplingPoints(i), 2.0) * r_M;

            this->mpLinearSolver->Solve( r_kdyn, rs, tmp_rhs );
            noalias(aux) = prod( r_M, rs );
        
            this->mpLinearSolver->Solve( r_kdyn, rAs, aux );
            noalias(aux) = prod( r_M, rAs );

            this->mpLinearSolver->Solve( r_kdyn, rAAs, aux );

            column( r_basis, (i*3) ) = rs;
            column( r_basis, (i*3)+1 ) = rAs;
            column( r_basis, (i*3)+2 ) = rAAs;
        }

        KRATOS_INFO_IF("Basis Construction Time", BaseType::GetEchoLevel() > 0 && rank == 0)
            << basis_construction_time.ElapsedSeconds() << std::endl;

        // orthogonalization step
        BuiltinTimer basis_orthogonalization_time;
        
        OrthogonalizationUtility::OrthogonalizeQR<TReducedDenseSpace>(r_basis);

        KRATOS_INFO_IF("Basis Orthogonalization Time", BaseType::GetEchoLevel() > 0 && rank == 0)
            << basis_orthogonalization_time.ElapsedSeconds() << std::endl;

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