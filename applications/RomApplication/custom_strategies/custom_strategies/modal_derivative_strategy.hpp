//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Altug Emiroglu, https://github.com/emiroglu
//
//

#if !defined(KRATOS_MODAL_DERIVATIVE_STRATEGY )
#define  KRATOS_MODAL_DERIVATIVE_STRATEGY

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/builtin_timer.h"
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "rom_application_variables.h"

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

/** @brief Solving strategy base class
 * @details This is the base class from which we will derive all the strategies (line-search, NR, etc...)
 */

template<class TSparseSpace,
         class TDenseSpace,
         class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         >
class ModalDerivativeStrategy
    : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    ///@name Type Definitions
    ///@{

//     typedef std::set<Dof::Pointer,ComparePDof>                                    DofSetType;

    // Base type definition
    typedef SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver> BaseType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;

    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::TSchemeType TSchemeType;

    typedef typename BaseType::TSchemeType::Pointer TSchemePointerType;

    typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;

    typedef typename BaseType::TBuilderAndSolverType::Pointer TBuilderAndSolverPointerType;

    typedef typename BaseType::TDofType TDofType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::DofIteratorType DofIteratorType;

    typedef typename BaseType::DofConstantIteratorType DofConstantIteratorType;

    typedef typename BaseType::NodesArrayType NodesArrayType;

    typedef typename BaseType::ElementsArrayType ElementsArrayType;

    typedef typename BaseType::ConditionsArrayType ConditionsArrayType;

    /** Counted pointer of ClassName */
    KRATOS_CLASS_POINTER_DEFINITION(ModalDerivativeStrategy);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    ModalDerivativeStrategy(
        ModelPart& rModelPart,
        TSchemePointerType pScheme,
        TBuilderAndSolverPointerType pBuilderAndSolver,
        bool DerivativeType
        )
        : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart)
    {
        KRATOS_TRY

        // assign the scheme
        mpScheme = pScheme;

        // assign the builder & solver
        mpBuilderAndSolver = pBuilderAndSolver;

        // assign the derivative type. false : static, true : dynamic
        mDerivativeType = false;

        // ensure initialization of system matrices in InitializeSolutionStep()
        mpBuilderAndSolver->SetDofSetIsInitializedFlag(false);

        // default echo level (mute)
        this->SetEchoLevel(0);

        // default rebuild level (build at each solution step)
        this->SetRebuildLevel(1);
        
        // TSystemMatrixType* Auxmp = new TSystemMatrixType;
        // mpA = Kratos::shared_ptr<TSystemMatrixType>(AuxpA);
        TSystemVectorType* AuxpInitialVariables = new TSystemVectorType;
        mpInitialVariables = Kratos::shared_ptr<TSystemVectorType>(AuxpInitialVariables);
        // TSystemVectorType* Auxpb = new TSystemVectorType;
        // mpb = Kratos::shared_ptr<TSystemVectorType>(Auxpb);

        mSolutionStepIsInitialized = false;

        rModelPart.GetProcessInfo()[DERIVATIVE_INDEX] = rModelPart.GetProcessInfo()[EIGENVALUE_VECTOR].size();
        
        KRATOS_CATCH("")
    }

    // /// Deleted copy constructor.
    // ModalDerivativeStrategy(const ModalDerivativeStrategy& Other) = delete;

    /// Destructor.
    ~ModalDerivativeStrategy() override
    {
        // Clear() controls order of deallocation to avoid invalid memory access
        // in some special cases.
        // warning: BaseType::GetModelPart() may be invalid here.
        this->Clear();
    }
    
    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * Initialization to be performed once before using the strategy.
     */
    void Initialize() override
    {
        KRATOS_TRY

        ModelPart& r_model_part = BaseType::GetModelPart();
        const int rank = r_model_part.GetCommunicator().MyPID();

        KRATOS_INFO_IF("ModalDerivativeStrategy", BaseType::GetEchoLevel() > 2 && rank == 0)
            <<  "Entering Initialize" << std::endl;

        if (mInitializeWasPerformed == false)
        {
            TSchemePointerType& p_scheme = this->pGetScheme();

            if (p_scheme->SchemeIsInitialized() == false)
                p_scheme->Initialize(r_model_part);

            if (p_scheme->ElementsAreInitialized() == false)
                p_scheme->InitializeElements(r_model_part);

            if (p_scheme->ConditionsAreInitialized() == false)
                p_scheme->InitializeConditions(r_model_part);

            mInitializeWasPerformed = true;
        }

        KRATOS_INFO_IF("ModalDerivativeStrategy", BaseType::GetEchoLevel() > 2 && rank == 0)
            <<  "Exiting Initialize" << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * @brief Clears the internal storage
     */
    void Clear() override
    {
        KRATOS_TRY;

        // if the preconditioner is saved between solves, it
        // should be cleared here.
        this->pGetBuilderAndSolver()->GetLinearSystemSolver()->Clear();

        if (mpA != nullptr){
            TSparseSpace::Clear(mpA);
            mpA = nullptr;
        }
        if (mpDx != nullptr)
        {
            TSparseSpace::Clear(mpDx);
            mpDx = nullptr;
        }
        if (mpb != nullptr)
        {
            TSparseSpace::Clear(mpb);
            mpb = nullptr;
        }

        //setting to zero the internal flag to ensure that the dof sets are recalculated
        this->pGetBuilderAndSolver()->Clear();
        this->pGetScheme()->Clear();

        mInitializeWasPerformed = false;
        mSolutionStepIsInitialized = false;

        KRATOS_CATCH("");
    }

    /**
     * @brief Performs all the required operations that should be done (for each step) before solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     */
    void InitializeSolutionStep() override
    {
        KRATOS_TRY

        if (mSolutionStepIsInitialized == false)
        {
            //pointers needed in the solution
            TSchemePointerType& p_scheme = this->pGetScheme();
            TBuilderAndSolverPointerType p_builder_and_solver = this->pGetBuilderAndSolver();
            
            const int rank = BaseType::GetModelPart().GetCommunicator().MyPID();

            //set up the system, operation performed just once unless it is required
            //to reform the dof set at each iteration
            BuiltinTimer system_construction_time;
            // if (p_builder_and_solver->GetDofSetIsInitializedFlag() == false ||
            //         mReformDofSetAtEachStep == true)
            if (p_builder_and_solver->GetDofSetIsInitializedFlag() == false)
            {
                //setting up the list of the DOFs to be solved
                BuiltinTimer setup_dofs_time;
                p_builder_and_solver->SetUpDofSet(p_scheme, BaseType::GetModelPart());
                KRATOS_INFO_IF("Setup Dofs Time", BaseType::GetEchoLevel() > 0 && rank == 0)
                    << setup_dofs_time.ElapsedSeconds() << std::endl;

                //shaping correctly the system
                BuiltinTimer setup_system_time;
                p_builder_and_solver->SetUpSystem(BaseType::GetModelPart());
                KRATOS_INFO_IF("Setup System Time", BaseType::GetEchoLevel() > 0 && rank == 0)
                    << setup_system_time.ElapsedSeconds() << std::endl;

                //setting up the Vectors involved to the correct size
                BuiltinTimer system_matrix_resize_time;
                p_builder_and_solver->ResizeAndInitializeVectors(p_scheme, mpA, mpDx, mpb,
                                                                 BaseType::GetModelPart());

                //setting up initial SolutionStepValue
                TSparseSpace::Resize(*mpInitialVariables, mpA->size1());
                TSparseSpace::SetToZero(*mpInitialVariables);
                this->StoreVariables();
                KRATOS_INFO_IF("System Matrix Resize Time", BaseType::GetEchoLevel() > 0 && rank == 0)
                    << system_matrix_resize_time.ElapsedSeconds() << std::endl;

                
            }

            KRATOS_INFO_IF("System Construction Time", BaseType::GetEchoLevel() > 0 && rank == 0)
                << system_construction_time.ElapsedSeconds() << std::endl;

            TSystemMatrixType& rA  = *mpA;
            TSystemVectorType& rDx = *mpDx;
            TSystemVectorType& rb  = *mpb;

            //initial operations ... things that are constant over the Solution Step
            p_builder_and_solver->InitializeSolutionStep(BaseType::GetModelPart(), rA, rDx, rb);

            //initial operations ... things that are constant over the Solution Step
            p_scheme->InitializeSolutionStep(BaseType::GetModelPart(), rA, rDx, rb);

            mSolutionStepIsInitialized = true;
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief Performs all the required operations that should be done (for each step) after solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     */
    void FinalizeSolutionStep() override
    {
        KRATOS_TRY;

        const int rank = BaseType::GetModelPart().GetCommunicator().MyPID();
        KRATOS_INFO_IF("ModalDerivativeStrategy", BaseType::GetEchoLevel() > 2 && rank == 0)
            <<  "Entering FinalizeSolutionStep" << std::endl;

        TSchemePointerType& p_scheme = this->pGetScheme();
        TBuilderAndSolverPointerType p_builder_and_solver = this->pGetBuilderAndSolver();

        TSystemMatrixType &rA  = *mpA;
        TSystemVectorType &rDx = *mpDx;
        TSystemVectorType &rb  = *mpb;

        //Finalisation of the solution step,
        //operations to be done after achieving convergence, for example the
        //Final Residual Vector (mb) has to be saved in there
        //to avoid error accumulation

        p_builder_and_solver->FinalizeSolutionStep(BaseType::GetModelPart(), rA, rDx, rb);
        p_scheme->FinalizeSolutionStep(BaseType::GetModelPart(), rA, rDx, rb);
        
        //Cleaning memory after the solution
        p_scheme->Clean();

        this->Clear();

        KRATOS_INFO_IF("ModalDerivativeStrategy", BaseType::GetEchoLevel() > 2 && rank == 0)
            <<  "Exiting FinalizeSolutionStep" << std::endl;

        KRATOS_CATCH("");
    }

    /**
     * @brief Solves the current step. This function returns true if a solution has been found, false otherwise.
     */
    bool SolveSolutionStep() override
    {
        KRATOS_TRY

        // Implementation of this function considers only the static derivatives

        ModelPart& r_model_part = BaseType::GetModelPart();

        TSchemePointerType& p_scheme = this->pGetScheme();
        TSystemMatrixType& rA = *mpA;
        TSystemVectorType& rb = *mpb;
        TSystemVectorType& rDx = *mpDx;

        // Get eigenvalues vector
        LocalSystemVectorType& r_eigenvalues = r_model_part.GetProcessInfo()[EIGENVALUE_VECTOR];
        const std::size_t num_eigenvalues = r_eigenvalues.size();

        // if static derivatives then build the stiffness matrix only once
        if (!mDerivativeType)
        {
            // Build stiffness matrix contribution
            r_model_part.GetProcessInfo()[BUILD_LEVEL] = 1;
            this->pGetBuilderAndSolver()->BuildLHS(p_scheme,r_model_part,rA);
        } else
        {
            KRATOS_ERROR <<"Invalid DerivativeType" << std::endl;
        }
        
        unsigned int basis_j_start_index;
        for (unsigned int basis_i = 0; basis_i < num_eigenvalues; basis_i++)
        {
            // Set the EIGENVALUE_I counter for use in the scheme
            r_model_part.GetProcessInfo()[EIGENVALUE_I] = basis_i;

            // Shift the derivative start index due to symmetry of static derivatives
            basis_j_start_index = basis_i; 
            
            // Reset all the dofs to initial value
            this->ResetVariables();
            
            // if dynamic derivatives then build system matrix at each step
            if (mDerivativeType)
            {
                KRATOS_ERROR <<"Invalid DerivativeType" << std::endl;
                // // Reset LHS for each iteration
                // TSparseSpace::SetToZero(rA);

                // // Build stiffness matrix contribution
                // r_model_part.GetProcessInfo()[BUILD_LEVEL] = 1;
                // this->pGetBuilderAndSolver()->BuildLHS(p_scheme,r_model_part,rA);
                
                // // Build mass matrix contribution
                // r_model_part.GetProcessInfo()[BUILD_LEVEL] = 2;
                // this->pGetBuilderAndSolver()->BuildLHS(p_scheme,r_model_part,rA);

                // // Reset the derivative start index due to unsymmetry of dynamic derivatives
                // basis_j_start_index = 0;
            }

            for (unsigned int basis_j = basis_j_start_index; basis_j < num_eigenvalues; basis_j++)
            {
                // Set the EIGENVALUE_J counter for use in the scheme
                r_model_part.GetProcessInfo()[EIGENVALUE_J] = basis_j;

                // Reset RHS and solution vector at each step
                TSparseSpace::SetToZero(rb);
                TSparseSpace::SetToZero(rDx);

                // Build RHS for static derivative contribution
                r_model_part.GetProcessInfo()[BUILD_LEVEL] = 1;
                this->pGetBuilderAndSolver()->BuildRHS(p_scheme, r_model_part, rb);
                // Build RHS for dynamic derivative contribution
                if (mDerivativeType)
                {
                    KRATOS_ERROR <<"Invalid DerivativeType" << std::endl;
                    // r_model_part.GetProcessInfo()[BUILD_LEVEL] = 2;
                    // this->pGetBuilderAndSolver()->BuildRHS(p_scheme, r_model_part, rb);
                }

                if(r_model_part.MasterSlaveConstraints().size() != 0) {
                    this->pGetBuilderAndSolver()->ApplyRHSConstraints(p_scheme, r_model_part, rb);
                }

                this->pGetBuilderAndSolver()->ApplyDirichletConditions(p_scheme, r_model_part, rA, rDx, rb);

                // Apply dynamic derivative constraint
                // this->ApplyDynamicDerivativeConstraint();

                // Solve the system
                this->pGetBuilderAndSolver()->SystemSolve(rA, rDx, rb);
                // this->pGetBuilderAndSolver()->SystemSolveWithPhysics(rA, rDx, rb, r_model_part);
                
                // Assign solution to ROM_BASIS
                this->AssignVariables(rDx);
            }
        }

        return true;
        KRATOS_CATCH("")
    }

    /**
     * @brief This function stores the initial SolutionStepValue.
     * @details
     * { 
     * If an a priori analysis (e.g. a nonlinear analysis) is performed 
     * the last equilibrium state is stored as the equilibrium state at which the 
     * linearization of the nonlinear problem is performed and a starting point
     * to compute modal derivatives.
     * } 
     */
    void StoreVariables()
    {
        KRATOS_TRY

        DofsArrayType& r_dof_set = this->pGetBuilderAndSolver()->GetDofSet();
        for (auto dof : r_dof_set)
        {
            (*mpInitialVariables)(dof.EquationId()) = dof.GetSolutionStepValue();
        }
        KRATOS_CATCH("")
    }

    /**
     * @brief This function resets the SolutionStepValue to the stored initial SolutionStepValue.
     * @details
     * { 
     * Modal derivatives are computed using a semi-analytic scheme using FD. 
     * Thus, the SolutionStepValue has to be reset to the initial value at each loop
     * } 
     */    
    void ResetVariables()
    {
        KRATOS_TRY
        DofsArrayType& r_dof_set = this->pGetBuilderAndSolver()->GetDofSet();
        for (auto dof : r_dof_set)
        {
            dof.GetSolutionStepValue() = (*mpInitialVariables)(dof.EquationId());
        }
        KRATOS_CATCH("")
    }

    void AssignVariables(TSystemVectorType& rDx)
    {
        KRATOS_TRY
        
        ModelPart& r_model_part = BaseType::GetModelPart();

        const auto& r_dof_set = this->pGetBuilderAndSolver()->GetDofSet();

        std::size_t derivative_index = r_model_part.GetProcessInfo()[DERIVATIVE_INDEX];
        
        for (ModelPart::NodeIterator itNode = r_model_part.NodesBegin(); itNode!= r_model_part.NodesEnd(); itNode++) {
            ModelPart::NodeType::DofsContainerType& NodeDofs = itNode->GetDofs();
            const std::size_t NumNodeDofs = NodeDofs.size();
            Matrix& rRomBasis = itNode->GetValue(ROM_BASIS);

            // fill the ROM_BASIS
            for (std::size_t iDOF = 0; iDOF < NumNodeDofs; iDOF++)
            {
                const auto itDof = std::begin(NodeDofs) + iDOF;
                bool is_active = !(r_dof_set.find(**itDof) == r_dof_set.end());
                if ((*itDof)->IsFree() && is_active) {
                   rRomBasis(iDOF,derivative_index) = rDx((*itDof)->EquationId());
                }
                else {
                   rRomBasis(iDOF,derivative_index) = 0.0;
                }
            }            
        }

        r_model_part.GetProcessInfo()[DERIVATIVE_INDEX] += 1;

        KRATOS_CATCH("")

    }

    /**
     * @brief This sets the level of echo for the solving strategy
     * @param Level of echo for the solving strategy
     * @details 
     * {
     * 0 -> Mute... no echo at all
     * 1 -> Printing time and basic informations
     * 2 -> Printing linear solver data
     * 3 -> Print of debug informations: Echo of stiffness matrix, Dx, b...
     * }
     */
    void SetEchoLevel(const int Level) override
    {
        BaseType::SetEchoLevel(Level);
        this->pGetBuilderAndSolver()->SetEchoLevel(Level);
    }

    TBuilderAndSolverPointerType& pGetBuilderAndSolver()
    {
        return mpBuilderAndSolver;
    };

    TSchemePointerType& pGetScheme()
    {
        return mpScheme;
    };

    /**
     * @brief Function to perform expensive checks.
     * @details It is designed to be called ONCE to verify that the input is correct.
     */
    int Check() override
    {
        KRATOS_TRY

        ModelPart& r_model_part = BaseType::GetModelPart();
        const int rank = r_model_part.GetCommunicator().MyPID();

        KRATOS_INFO_IF("ModalDerivativeStrategy", BaseType::GetEchoLevel() > 2 && rank == 0)
            <<  "Entering Check" << std::endl;

        // check the model part
        BaseType::Check();

        // check the scheme
        this->pGetScheme()->Check(r_model_part);

        // check the builder and solver
        this->pGetBuilderAndSolver()->Check(r_model_part);

        KRATOS_INFO_IF("ModalDerivativeStrategy", BaseType::GetEchoLevel() > 2 && rank == 0)
            <<  "Exiting Check" << std::endl;

        return 0;

        KRATOS_CATCH("")
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ModalDerivativeStrategy";
    }

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

private:

    ///@}
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    TSchemePointerType mpScheme;

    TBuilderAndSolverPointerType mpBuilderAndSolver;

    bool mDerivativeType;

    TSystemMatrixPointerType mpA;

    TSystemVectorPointerType mpDx;
    
    TSystemVectorPointerType mpb;

    TSystemMatrixPointerType mpMassMatrix;

    TSystemMatrixPointerType mpStiffnessMatrix;

    TSystemVectorPointerType mpInitialVariables;

    bool mSolutionStepIsInitialized;

    bool mInitializeWasPerformed;

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

    /** Copy constructor.
     */
    ModalDerivativeStrategy(const ModalDerivativeStrategy& Other);


    ///@}

}; /* Class ModalDerivativeStrategy */

///@}

///@name Type Definitions
///@{


///@}

} /* namespace Kratos.*/

#endif /* KRATOS_MODAL_DERIVATIVE_STRATEGY  defined */
