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

#if !defined(MOR_ONLINE_STRATEGY)
#define MOR_ONLINE_STRATEGY

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "utilities/builtin_timer.h"

//default builder and solver
#include "custom_strategies/custom_builder_and_solvers/system_matrix_builder_and_solver.hpp"
#include "custom_strategies/custom_strategies/mor_offline_strategy.hpp"

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
 * @class MorOnlineStrategy
 * @ingroup KratosCore
 * @brief This is the MOR online strategy
 * @details The strategy calls an offline MOR strategy and solves the reduced system
 * @author Quirin Aumann
 */
template <class TSparseSpace,
          class TDenseSpace,  // = DenseSpace<double>,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
class MorOnlineStrategy
    : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
  public:
    ///@name Type Definitions
    ///@{
    // Counted pointer of ClassName
    KRATOS_CLASS_POINTER_DEFINITION(MorOnlineStrategy);

    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    // typedef MassAndStiffnessBuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver > TBuilderAndSolverType;

    // typedef LinearMorMatrixOutputStrategy< TSparseSpace, TDenseSpace, TLinearSolver > OfflineStrategyType;

    typedef MorOfflineStrategy< TSparseSpace, TDenseSpace, TLinearSolver > OfflineStrategyType;

    typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;

    typedef typename BaseType::TDataType TDataType;

    typedef TSparseSpace SparseSpaceType;

    typedef typename BaseType::TSchemeType TSchemeType;

    //typedef typename BaseType::DofSetType DofSetType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;

    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;

    ///@}
    ///@name Life Cycle

    ///@{

    /**
     * Default constructor
     * @param rModelPart The model part of the problem
     * @param pLinearSolver The linear solver
     * @param pOfflineStrategy The offline MOR strategy
     */
    MorOnlineStrategy(
        ModelPart& rModelPart,
        typename TLinearSolver::Pointer pLinearSolver,
        typename OfflineStrategyType::Pointer pOfflineStrategy)
        : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, false)
    {
        KRATOS_TRY;

        // Setting up the default scheme
        mpScheme = typename TSchemeType::Pointer(
            new TSchemeType());

        // Setting up the default builder and solver
        mpBuilderAndSolver = typename TBuilderAndSolverType::Pointer(
            new TBuilderAndSolverType(pLinearSolver));

        // Saving the linear solver
        mpLinearSolver = pLinearSolver;

        // Setting the offline strategy
        mpOfflineStrategy = pOfflineStrategy;

        // Set flags to correcty start the calculations
        mSolutionStepIsInitialized = false;
        mInitializeWasPerformed = false;
        mReformDofSetAtEachStep = false;

        // Tells to the builder and solver if the reactions have to be Calculated or not
        GetBuilderAndSolver()->SetCalculateReactionsFlag(false);

        // Tells to the Builder And Solver if the system matrix and vectors need to
        // be reshaped at each step or not
        GetBuilderAndSolver()->SetReshapeMatrixFlag(mReformDofSetAtEachStep);

        // Set EchoLevel to the default value (only time is displayed)
        SetEchoLevel(1);

        // Do not rebuild the stiffness matrix
        this->SetRebuildLevel(0);

        mpA = TSparseSpace::CreateEmptyMatrixPointer(); //stiffness
        mpM = TSparseSpace::CreateEmptyMatrixPointer(); //mass
        mpS = TSparseSpace::CreateEmptyMatrixPointer(); //Damping
        mpRHS = TSparseSpace::CreateEmptyVectorPointer();

        KRATOS_CATCH("");
    }

    /**
     * @brief Destructor.
     * @details In trilinos third party library, the linear solver's preconditioner should be freed before the system matrix. We control the deallocation order with Clear().
     */
    ~MorOnlineStrategy() override
    {
        Clear();
    }

    /**
     * @brief Set method for the time scheme
     * @param pScheme The pointer to the time scheme considered
     */
    void SetScheme(typename TSchemeType::Pointer pScheme)
    {
        mpScheme = pScheme;
    };

    /**
     * @brief Get method for the time scheme
     * @return mpScheme: The pointer to the time scheme considered
     */
    typename TSchemeType::Pointer GetScheme()
    {
        return mpScheme;
    };

    /**
     * @brief Set method for the builder and solver
     * @param pNewBuilderAndSolver The pointer to the builder and solver considered
     */
    void SetBuilderAndSolver(typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver)
    {
        mpBuilderAndSolver = pNewBuilderAndSolver;
    };

    /**
     * @brief Get method for the builder and solver
     * @return mpBuilderAndSolver: The pointer to the builder and solver considered
     */
    typename TBuilderAndSolverType::Pointer GetBuilderAndSolver()
    {
        return mpBuilderAndSolver;
    };

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

    /**
     * @brief This method sets the flag mReformDofSetAtEachStep
     * @param Flag The flag that tells if each time step the system is rebuilt
     */
    void SetReformDofSetAtEachStepFlag(bool Flag)
    {
        mReformDofSetAtEachStep = Flag;
        GetBuilderAndSolver()->SetReshapeMatrixFlag(mReformDofSetAtEachStep);
    }

    /**
     * @brief This method returns the flag mReformDofSetAtEachStep
     * @return The flag that tells if each time step the system is rebuilt
     */
    bool GetReformDofSetAtEachStepFlag()
    {
        return mReformDofSetAtEachStep;
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
        // GetBuilderAndSolver()->SetEchoLevel(Level);
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
            typename OfflineStrategyType::Pointer p_offline_strategy = GetOfflineStrategy();
            p_offline_strategy->Initialize();
            p_offline_strategy->Solve();

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

        // if the preconditioner is saved between solves, it
        // should be cleared here.
        // GetBuilderAndSolver()->GetLinearSystemSolver()->Clear();

        if (mpA != nullptr)
            SparseSpaceType::Clear(mpA);
        if (mpM != nullptr)
            SparseSpaceType::Clear(mpM);
        if (mpS != nullptr)
            SparseSpaceType::Clear(mpS);
        if (mpRHS != nullptr)
            SparseSpaceType::Clear(mpRHS);

        //setting to zero the internal flag to ensure that the dof sets are recalculated
        // GetBuilderAndSolver()->SetDofSetIsInitializedFlag(false);
        // GetBuilderAndSolver()->Clear();
        GetScheme()->Clear();

        mInitializeWasPerformed = false;

        KRATOS_CATCH("");
    }

    /**
     * @brief Solves the current step. This function returns true if a solution has been found, false otherwise.
     */
    bool SolveSolutionStep() override
    {
        // typename TSchemeType::Pointer p_scheme = GetScheme();
        // typename TBuilderAndSolverType::Pointer p_builder_and_solver = GetBuilderAndSolver();
        // TSystemMatrixType& rA  = *mpA;
        // TSystemMatrixType& rM = *mpM;
        // TSystemVectorType& rRHS  = *mpRHS;

        // TSystemVectorType tmp(rA.size1(), 0.0);

        // p_builder_and_solver->BuildStiffnessMatrix(p_scheme, BaseType::GetModelPart(), rA, tmp);

        // // Applying the Dirichlet Boundary conditions
        // p_builder_and_solver->ApplyDirichletConditions(p_scheme, BaseType::GetModelPart(), rA, tmp, rRHS);  

        // p_builder_and_solver->BuildMassMatrix(p_scheme, BaseType::GetModelPart(), rM, tmp);

        // p_builder_and_solver->ApplyDirichletConditionsForMassMatrix(p_scheme, BaseType::GetModelPart(), rM);

        // p_builder_and_solver->BuildRHS(p_scheme, BaseType::GetModelPart(), rRHS);

        // EchoInfo(0);

        typename OfflineStrategyType::Pointer p_offline_strategy = GetOfflineStrategy();
        auto& r_process_info = BaseType::GetModelPart().GetProcessInfo();
        double excitation_frequency = r_process_info[TIME];

        TSystemMatrixType r_Mr = p_offline_strategy->GetMr();
        TSystemMatrixType r_Ar = p_offline_strategy->GetAr();
        TSystemVectorType r_RHSr = p_offline_strategy->GetRHSr();
        TSystemMatrixType r_basis = p_offline_strategy->GetBasis();
        TSystemMatrixType r_Sr = p_offline_strategy->GetSr();

        const unsigned int system_size = r_basis.size1();
        const unsigned int system_size_r = r_basis.size2();
        // KRATOS_WATCH(system_size)
        // KRATOS_WATCH(system_size_r)

        TSystemVectorType displacement_reduced;
        displacement_reduced.resize( system_size_r, false );
        displacement_reduced = ZeroVector( system_size_r );
        auto kdyn = SparseSpaceType::CreateEmptyMatrixPointer();
        auto& r_kdyn = *kdyn;
        SparseSpaceType::Resize(r_kdyn, system_size_r, system_size_r);

        r_kdyn = r_Ar - ( std::pow( excitation_frequency, 2.0 ) * r_Mr );           // Without Damping
       // r_kdyn = r_Ar - ( std::pow( excitation_frequency, 2.0 ) * r_Mr ) - r_Sr;    // With Damping

        GetBuilderAndSolver()->GetLinearSystemSolver()->Solve( r_kdyn, displacement_reduced, r_RHSr );

        TSystemVectorType displacement;
        displacement.resize( system_size, false );
        displacement = ZeroVector( system_size );

        displacement = prod( r_basis, displacement_reduced );
        AssignVariables(displacement);

        return false;
    }

    /**
     * @brief Function to perform expensive checks.
     * @details It is designed to be called ONCE to verify that the input is correct.
     */
    int Check() override
    {
        KRATOS_TRY

        BaseType::Check();

        GetBuilderAndSolver()->Check(BaseType::GetModelPart());

        GetScheme()->Check(BaseType::GetModelPart());

        return 0;

        KRATOS_CATCH("")
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

    /**
     * @brief This method returns the LHS matrix
     * @return The LHS matrix
     */
    TSystemMatrixType &GetSystemMatrix()
    {
        TSystemMatrixType &mA = *mpA;

        return mA;
    }

    /**
     * @brief This method returns the RHS vector
     * @return The RHS vector
     */
    TSystemVectorType& GetSystemVector()
    {
        TSystemVectorType& mb = *mpRHS;

        return mb;
    }

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
    typename TSchemeType::Pointer mpScheme; /// The pointer to the time scheme employed
    typename TBuilderAndSolverType::Pointer mpBuilderAndSolver; /// The pointer to the builder and solver employe
    typename OfflineStrategyType::Pointer mpOfflineStrategy;

    TSystemVectorPointerType mpRHS; /// The RHS vector of the system of equations
    TSystemMatrixPointerType mpA; /// The Stiffness matrix of the system of equations
    TSystemMatrixPointerType mpM; /// The Mass matrix of the system of equations
    TSystemMatrixPointerType mpS; /// The Damping matrix of the system of equations
    
    // reduced matrices
    // TSystemVectorPointerType mpRHSr; //reduced RHS
    // TSystemMatrixPointerType mpAr;
    // TSystemMatrixPointerType mpMr;
    // TSystemMatrixPointerType mpBasis;

    /**
     * @brief Flag telling if it is needed to reform the DofSet at each
    solution step or if it is possible to form it just once
    * @details Default = false
        - true  : Reforme at each time step
        - false : Form just one (more efficient)
     */
    bool mReformDofSetAtEachStep;

    bool mSolutionStepIsInitialized; /// Flag to set as initialized the solution step

    bool mInitializeWasPerformed; /// Flag to set as initialized the strategy

    ///@}
    ///@name Private Operators
    ///@{

    /**
     * @brief This method returns the components of the system of equations depending of the echo level
     * @param IterationNumber The non linear iteration in the solution loop
     */
    virtual void EchoInfo(const unsigned int IterationNumber)
    {
        TSystemMatrixType& rA  = *mpA;
        TSystemMatrixType& rM = *mpM;
        TSystemVectorType& rRHS  = *mpRHS;
        TSystemMatrixType& rS  = *mpS;

        if (this->GetEchoLevel() == 2) //if it is needed to print the debug info
        {
            KRATOS_INFO("RHS") << "RHS  = " << rRHS << std::endl;
        }
        else if (this->GetEchoLevel() == 3) //if it is needed to print the debug info
        {
            KRATOS_INFO("LHS") << "SystemMatrix = " << rA << std::endl;
            KRATOS_INFO("Dx")  << "Mass Matrix = " << mpM << std::endl;
            KRATOS_INFO("Sx") << "Damping Matrix = " << rS << std::endl;
            KRATOS_INFO("RHS") << "RHS  = " << rRHS << std::endl;
        }
        std::stringstream matrix_market_name;
        matrix_market_name << "A_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << "_" << IterationNumber << ".mm";
        TSparseSpace::WriteMatrixMarketMatrix((char *)(matrix_market_name.str()).c_str(), rA, false);

        std::stringstream matrix_market_mass_name;
        matrix_market_mass_name << "M_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << "_" << IterationNumber << ".mm";
        TSparseSpace::WriteMatrixMarketMatrix((char *)(matrix_market_mass_name.str()).c_str(), rM, false); 

        std::stringstream matrix_market_damping_name;
        matrix_market_name << "S_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << "_" << IterationNumber << ".mm";
        TSparseSpace::WriteMatrixMarketMatrix((char *)(matrix_market_damping_name.str()).c_str(), rS, false);           

        std::stringstream matrix_market_vectname;
        matrix_market_vectname << "RHS_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << "_" << IterationNumber << ".mm.rhs";
        TSparseSpace::WriteMatrixMarketVector((char *)(matrix_market_vectname.str()).c_str(), rRHS);
    }

    /**
     * @brief This method prints information after reach the max number of iterations
     * @todo Replace by logger
     */

    virtual void MaxIterationsExceeded()
    {
        if (this->GetEchoLevel() != 0 && BaseType::GetModelPart().GetCommunicator().MyPID() == 0)
        {
            std::cout << "***************************************************" << std::endl;
            std::cout << "******* ATTENTION: max iterations exceeded ********" << std::endl;
            std::cout << "***************************************************" << std::endl;
        }
    }

    void AssignVariables(TSystemVectorType& rDisplacement, int step=0)
    {
        auto& r_model_part = BaseType::GetModelPart();
        for( auto& node : r_model_part.Nodes() )
        {
            ModelPart::NodeType::DofsContainerType& rNodeDofs = node.GetDofs();
            
            for( auto it_dof = std::begin(rNodeDofs); it_dof != std::end(rNodeDofs); it_dof++ )
            {
                if( !it_dof->IsFixed() )
                {
                    it_dof->GetSolutionStepValue(step) = rDisplacement(it_dof->EquationId());
                }
                else
                {
                    it_dof->GetSolutionStepValue(step) = 0.0;
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

    MorOnlineStrategy(const MorOnlineStrategy &Other){};

    ///@}

}; /* Class MorOnlineStrategy */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos. */

#endif /* MOR_ONLINE_STRATEGY  defined */