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

#if !defined(MOR_ONLINE_SECOND_ORDER_STRATEGY)
#define MOR_ONLINE_SECOND_ORDER_STRATEGY

// System includes

// External includes

// Project includes
#include "solving_strategies/strategies/solving_strategy.h"
#include "utilities/builtin_timer.h"

// Linear solvers
#include "linear_solvers/linear_solver.h"

//default builder and solver
#include "custom_strategies/custom_builder_and_solvers/system_matrix_builder_and_solver.hpp"
#include "custom_strategies/custom_strategies/mor_offline_second_order_strategy.hpp"

// Application includes
#include "mor_application_variables.h"

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
 * @class MorOnlineSecondOrderStrategy
 * @ingroup KratosCore
 * @brief This is the MOR online strategy
 * @details The strategy calls an offline MOR strategy and solves the reduced system
 * @author Quirin Aumann
 */
template <class TSparseSpace,
          class TDenseSpace,  // = DenseSpace<double>,
          class TLinearSolver, //= LinearSolver<TSparseSpace,TDenseSpace>
          class OfflineStrategyType
          >
class MorOnlineSecondOrderStrategy
    : public SolvingStrategy<TSparseSpace, TDenseSpace, LinearSolver<TSparseSpace,TDenseSpace>>
{
  public:
    ///@name Type Definitions
    ///@{
    // Counted pointer of ClassName
    KRATOS_CLASS_POINTER_DEFINITION(MorOnlineSecondOrderStrategy);

    typedef SolvingStrategy<TSparseSpace, TDenseSpace, LinearSolver<TSparseSpace,TDenseSpace>> BaseType;

    typedef TUblasSparseSpace< std::complex< double > > ComplexSparseSpaceType;
    typedef TUblasDenseSpace< std::complex< double > > ComplexDenseSpaceType;

    typedef typename ComplexDenseSpaceType::MatrixType ComplexDenseMatrixType;
    typedef typename ComplexDenseSpaceType::VectorType ComplexDenseVectorType;
    typedef typename ComplexDenseSpaceType::MatrixPointerType ComplexDenseMatrixPointerType;
    typedef typename ComplexDenseSpaceType::VectorPointerType ComplexDenseVectorPointerType;

    typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;

    typedef typename BaseType::TDataType TDataType;

    typedef TSparseSpace SparseSpaceType;

    typedef typename BaseType::TSchemeType TSchemeType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;

    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;

    typedef typename TDenseSpace::MatrixPointerType LocalSystemMatrixPointerType;

    typedef typename TDenseSpace::VectorPointerType LocalSystemVectorPointerType;

    typedef LinearSolver<ComplexSparseSpaceType, ComplexDenseSpaceType> ComplexLinearSolverType;


    ///@}
    ///@name Life Cycle

    ///@{

    /**
     * Default constructor
     * TODO: rename to MorOnlineSecondOrderStrategy
     * @param rModelPart The model part of the problem
     * @param pLinearSolver The linear solver
     * @param pOfflineStrategy The offline MOR strategy
     */
    MorOnlineSecondOrderStrategy(
        ModelPart& rModelPart,
        typename TLinearSolver::Pointer pLinearSolver,
        typename OfflineStrategyType::Pointer pOfflineStrategy,
        bool expandSolution = true)
        : SolvingStrategy<TSparseSpace, TDenseSpace, LinearSolver<TSparseSpace,TDenseSpace>>(rModelPart, false),
            mExpandSolution(expandSolution)
    {
        KRATOS_TRY;

        // Saving the linear solver
        mpLinearSolver = pLinearSolver;

        // Setting the offline strategy
        mpOfflineStrategy = pOfflineStrategy;

        // Set flags to correcty start the calculations
        mSolutionStepIsInitialized = false;
        mInitializeWasPerformed = false;

        // Set EchoLevel to the default value (only time is displayed)
        SetEchoLevel(0);

        // Do not rebuild the stiffness matrix
        this->SetRebuildLevel(0);

        KRATOS_CATCH("");
    }

    /**
     * @brief Destructor.
     * @details In trilinos third party library, the linear solver's preconditioner should be freed before the system matrix. We control the deallocation order with Clear().
     */
    ~MorOnlineSecondOrderStrategy() override
    {
        Clear();
    }

    void SetLinearSolver(typename TLinearSolver::Pointer pNewLinearSolver)
    {
        mpLinearSolver = pNewLinearSolver;
    }

    typename TLinearSolver::Pointer GetLinearSolver()
    {
        return mpLinearSolver;
    }

    void SetOfflineStrategy(typename OfflineStrategyType::Pointer pNewOfflineStrategy)
    {
        mpOfflineStrategy = pNewOfflineStrategy;
    };

    typename OfflineStrategyType::Pointer GetOfflineStrategy()
    {
        return mpOfflineStrategy;
    };

    /**
     * @brief This method sets the flag mInitializeWasPerformed
     * @param InitializePerformedFlag The flag that tells if the initialize has been computed
     */
    void SetInitializePerformedFlag(bool InitializePerformedFlag = true)
    {
        mInitializeWasPerformed = InitializePerformedFlag;
    }

    /**
     * @brief This method gets the flag mInitializeWasPerformed
     * @return mInitializeWasPerformed: The flag that tells if the initialize has been computed
     */
    bool GetInitializePerformedFlag()
    {
        return mInitializeWasPerformed;
    }

    std::complex<double> GetScalarResult()
    {
        return mScalarResult;
    }

    /**
     * @brief It sets the level of echo for the solving strategy
     * @param Level The level to set
     * @details The different levels of echo are:
     * - 0: Mute... no echo at all
     * - 1: Printing time and basic informations
     * - 2: Printing linear solver data
     * - 3: Print of debug informations: Echo of stiffness matrix, Dx, b...
     */

    void SetEchoLevel(int Level) override
    {
        BaseType::mEchoLevel = Level;
    }

    //*********************************************************************************
    /**OPERATIONS ACCESSIBLE FROM THE INPUT: **/

    /**
     * @brief Initialization and solve of the offline strategy
     */
    void Initialize() override
    {
        KRATOS_TRY;

        if (mInitializeWasPerformed == false)
        {
            std::cout << "INITIALIZE online\n";
            typename OfflineStrategyType::Pointer p_offline_strategy = GetOfflineStrategy();
            p_offline_strategy->Initialize();
            p_offline_strategy->Solve();

            // cast rhs to complex
            LocalSystemVectorPointerType p_tmp_rhs = p_offline_strategy->pGetRHSr();
            mpRHS = ComplexDenseSpaceType::CreateEmptyVectorPointer();
            *mpRHS = ComplexDenseVectorType( *p_tmp_rhs );

            mSystemSizeR = p_offline_strategy->pGetKr()->size1();

            mpA = ComplexDenseSpaceType::CreateEmptyMatrixPointer();
            mpA->resize( mSystemSizeR, mSystemSizeR, false );
            mpDx = ComplexDenseSpaceType::CreateEmptyVectorPointer();
            mpDx->resize( mSystemSizeR, false );

            if (mExpandSolution)
            {
                mSystemSize = p_offline_strategy->pGetBasis()->size1();
                mpResult = ComplexDenseSpaceType::CreateEmptyVectorPointer();
                mpResult->resize( mSystemSize, false );
            }

            mInitializeWasPerformed = true;
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief Clears the internal storage
     */
    void Clear() override
    {
        KRATOS_TRY;

        mpLinearSolver->Clear();

        if (mpA != nullptr)
            ComplexDenseSpaceType::Clear(mpA);
        if (mpRHS != nullptr)
            ComplexDenseSpaceType::Clear(mpRHS);
        if (mpDx != nullptr)
            ComplexDenseSpaceType::Clear(mpDx);
        if (mpResult != nullptr)
            ComplexDenseSpaceType::Clear(mpResult);

        mInitializeWasPerformed = false;

        KRATOS_CATCH("");
    }

    void InitializeSolutionStep() override
    {
        KRATOS_TRY;
        ModelPart& r_model_part = BaseType::GetModelPart();
        const int rank = r_model_part.GetCommunicator().MyPID();

        typename OfflineStrategyType::Pointer p_offline_strategy = GetOfflineStrategy();

        auto& r_process_info = r_model_part.GetProcessInfo();
        const double excitation_frequency = r_process_info[FREQUENCY];

        LocalSystemMatrixType& r_Kr = *(p_offline_strategy->pGetKr());
        LocalSystemMatrixType& r_Dr = *(p_offline_strategy->pGetDr());
        LocalSystemMatrixType& r_Mr = *(p_offline_strategy->pGetMr());
        ComplexDenseMatrixType& r_A = *mpA;

        BuiltinTimer system_construction_time;
        r_A = r_Dr;
        r_A *= std::complex<double>(0,excitation_frequency);
        r_A += r_Kr;
        r_A -= std::pow( excitation_frequency, 2.0 ) * r_Mr;

        KRATOS_INFO_IF("System Construction Time", BaseType::GetEchoLevel() > 0 && rank == 0)
            << system_construction_time.ElapsedSeconds() << std::endl;

        KRATOS_CATCH("");
    }

    /**
     * @brief Solves the current step. This function returns true if a solution has been found, false otherwise.
     */
    bool SolveSolutionStep() override
    {
        KRATOS_TRY;
        ModelPart& r_model_part = BaseType::GetModelPart();
        const int rank = r_model_part.GetCommunicator().MyPID();
        
        ComplexDenseVectorType& r_rhs = *mpRHS;
        ComplexDenseMatrixType& r_A = *mpA;
        ComplexDenseVectorType& r_dx = *mpDx;

        BuiltinTimer system_solve_time;

        mpLinearSolver->Solve( r_A, r_dx, r_rhs );

        KRATOS_INFO_IF("System Solve Time", BaseType::GetEchoLevel() > 0 && rank == 0)
            << system_solve_time.ElapsedSeconds() << std::endl;

        return true;
        KRATOS_CATCH("")
    }

    void FinalizeSolutionStep() override
    {
        KRATOS_TRY;
        ModelPart& r_model_part = BaseType::GetModelPart();
        const int rank = r_model_part.GetCommunicator().MyPID();

        if (mExpandSolution)
        {
            BuiltinTimer solution_expansion_time;
            ComplexDenseVectorType& r_result = *mpResult;

            r_result = prod( GetOfflineStrategy()->GetBasis(), *mpDx );
            AssignVariables(r_result);

            KRATOS_INFO_IF("Solution Expansion Time", BaseType::GetEchoLevel() > 0 && rank == 0)
                << solution_expansion_time.ElapsedSeconds() << std::endl;
        }
        else
        {
            auto& r_output_vector = mpOfflineStrategy->GetOVr();
            mScalarResult = inner_prod( r_output_vector, *mpDx );
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief Function to perform expensive checks.
     * @details It is designed to be called ONCE to verify that the input is correct.
     */
    int Check() override
    {
        KRATOS_TRY

        BaseType::Check();

        return 0;

        KRATOS_CATCH("")
    }

    /**
     * @brief This method returns the components of the system of equations depending of the echo level
     */
    virtual void EchoInfo()
    {
        // TSystemMatrixType& rA  = *mpKr;
        // TSystemMatrixType& rM = *mpMr;
        // TSystemVectorType& rRHS  = *mpRHS;
        // TSystemMatrixType& rS  = *mpDr;

        // if (this->GetEchoLevel() == 2) //if it is needed to print the debug info
        // {
        //     KRATOS_INFO("RHS") << "RHS  = " << rRHS << std::endl;
        // }
        // else if (this->GetEchoLevel() == 3) //if it is needed to print the debug info
        // {
        //     KRATOS_INFO("LHS") << "SystemMatrix = " << rA << std::endl;
        //     KRATOS_INFO("Dx")  << "Mass Matrix = " << mpMr << std::endl;
        //     KRATOS_INFO("Sx") << "Damping Matrix = " << rS << std::endl;
        //     KRATOS_INFO("RHS") << "RHS  = " << rRHS << std::endl;
        // }
        // std::stringstream matrix_market_name;
        // matrix_market_name << "A_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << "_" << IterationNumber << ".mm";
        // TSparseSpace::WriteMatrixMarketMatrix((char *)(matrix_market_name.str()).c_str(), rA, false);

        // std::stringstream matrix_market_mass_name;
        // matrix_market_mass_name << "M_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << "_" << IterationNumber << ".mm";
        // TSparseSpace::WriteMatrixMarketMatrix((char *)(matrix_market_mass_name.str()).c_str(), rM, false);

        // std::stringstream matrix_market_damping_name;
        // matrix_market_name << "S_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << "_" << IterationNumber << ".mm";
        // TSparseSpace::WriteMatrixMarketMatrix((char *)(matrix_market_damping_name.str()).c_str(), rS, false);

        // std::stringstream matrix_market_vectname;
        // matrix_market_vectname << "RHS_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << ".mm.rhs";
        // TSparseSpace::WriteMatrixMarketVector((char *)(matrix_market_vectname.str()).c_str(), *mpRHS);
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
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}

  private:
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

  protected:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{
    typename TLinearSolver::Pointer mpLinearSolver; /// The pointer to the linear solver considered
    typename OfflineStrategyType::Pointer mpOfflineStrategy;

    bool mExpandSolution;

    ComplexDenseMatrixPointerType mpA; /// The system matrix in reduced space
    ComplexDenseVectorPointerType mpDx; /// The result in reduced space
    ComplexDenseVectorPointerType mpRHS; /// The RHS vector in reduced space

    ComplexDenseVectorPointerType mpResult; /// The result expanded to original space
    std::complex<double> mScalarResult;

    unsigned int mSystemSize;
    unsigned int mSystemSizeR;

    bool mSolutionStepIsInitialized; /// Flag to set as initialized the solution step

    bool mInitializeWasPerformed; /// Flag to set as initialized the strategy

    ///@}
    ///@name Private Operators
    ///@{

    void AssignVariables(ComplexDenseVectorType& rDisplacement, int step=0)
    {
        // KRATOS_WATCH(rDisplacement)
        auto& r_model_part = BaseType::GetModelPart();
        for( auto& node : r_model_part.Nodes() )
        {
            ModelPart::NodeType::DofsContainerType& rNodeDofs = node.GetDofs();

            for( auto it_dof = std::begin(rNodeDofs); it_dof != std::end(rNodeDofs); it_dof++ )
            {
                if( !(*it_dof)->IsFixed() )
                {
                    (*it_dof)->GetSolutionStepValue(step) = std::real( rDisplacement((*it_dof)->EquationId()) );
                    (*it_dof)->GetSolutionStepReactionValue(step) = std::imag( rDisplacement((*it_dof)->EquationId()) );
                }
                else
                {
                    (*it_dof)->GetSolutionStepValue(step) = 0.0;
                }
            }
        }
    }

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

    MorOnlineSecondOrderStrategy(const MorOnlineSecondOrderStrategy &Other){};

    ///@}

}; /* Class MorOnlineSecondOrderStrategy */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos. */

#endif /* MOR_ONLINE_SECOND_ORDER_STRATEGY  defined */