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
#include "utilities/timer.h"
#include "utilities/openmp_utils.h"
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
        bool DerivativeTypeFlag,
        bool MassOrthonormalizeFlag
        )
        : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart),
        mDerivativeTypeFlag(DerivativeTypeFlag), 
        mMassOrthonormalizeFlag(MassOrthonormalizeFlag)
    {
        KRATOS_TRY

        // assign the scheme
        mpScheme = pScheme;

        // assign the builder & solver
        mpBuilderAndSolver = pBuilderAndSolver;

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
        if (mpStiffnessMatrix != nullptr){
            TSparseSpace::Clear(mpStiffnessMatrix);
            mpStiffnessMatrix = nullptr;
        }
        if (mpMassMatrix != nullptr){
            TSparseSpace::Clear(mpMassMatrix);
            mpMassMatrix = nullptr;
        }
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
        if (mpInitialVariables != nullptr)
        {
            TSparseSpace::Clear(mpInitialVariables);
            mpInitialVariables = nullptr;
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
                if (mDerivativeTypeFlag || mMassOrthonormalizeFlag)
                    p_builder_and_solver->ResizeAndInitializeVectors(p_scheme, mpMassMatrix, mpDx, mpb,
                                                                 BaseType::GetModelPart());
                if (mDerivativeTypeFlag)
                    p_builder_and_solver->ResizeAndInitializeVectors(p_scheme, mpStiffnessMatrix, mpDx, mpb,
                                                                 BaseType::GetModelPart());

                //setting up initial SolutionStepValue
                DofsArrayType& r_dof_set = p_builder_and_solver->GetDofSet();
                TSparseSpace::Resize(*mpInitialVariables, r_dof_set.size());
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

        ModelPart& r_model_part = BaseType::GetModelPart();
        TSchemePointerType& p_scheme = this->pGetScheme();
        TSystemMatrixType& rA = *mpA;
        TSystemMatrixType& rStiffnessMatrix = *mpStiffnessMatrix;
        TSystemMatrixType& rMassMatrix = *mpMassMatrix;
        TSystemVectorType& rb = *mpb;
        TSystemVectorType& rDx = *mpDx;
        TSystemVectorType basis;
        TSparseSpace::Resize(basis, rA.size1());
       
        // Get eigenvalues vector
        LocalSystemVectorType& r_eigenvalues = r_model_part.GetProcessInfo()[EIGENVALUE_VECTOR];
        const std::size_t num_eigenvalues = r_eigenvalues.size();

        // Build system matrices
        // Build mass matrix
        if (mMassOrthonormalizeFlag || mDerivativeTypeFlag)
        {
            r_model_part.GetProcessInfo()[BUILD_LEVEL] = 1;
            this->pGetBuilderAndSolver()->BuildLHS(p_scheme,r_model_part,rMassMatrix);
        }

        // Build stiffness matrix
        r_model_part.GetProcessInfo()[BUILD_LEVEL] = 2;
        if (mDerivativeTypeFlag)
        {   // Dynamic derivatives
            
            // Build stiffness matrix separately            
            this->pGetBuilderAndSolver()->BuildLHS(p_scheme,r_model_part,rStiffnessMatrix);            
        } 
        else
        {   // Static derivatives
            
            // Build the stiffness matrix into system matrix directly
            this->pGetBuilderAndSolver()->BuildLHS(p_scheme,r_model_part,rA);
        }
        
        // Derivative of basis_i
        std::size_t basis_j_start_index;
        for (std::size_t basis_i = 0; basis_i < num_eigenvalues; basis_i++)
        {
            // Set the BASIS_I counter for use in the scheme
            r_model_part.GetProcessInfo()[BASIS_I] = basis_i;

            if (mDerivativeTypeFlag)
            {   // Dynamic derivatives

                const double eigenvalue_i = r_eigenvalues[basis_i];

                // If dynamic derivatives then build system matrix for each eigenvalue
                rA = rStiffnessMatrix - (eigenvalue_i * rMassMatrix);
                
                // Dynamic derivatives are unsymmetric
                basis_j_start_index = 0;

                // Get basis_i
                this->GetBasis(basis_i, basis);
            } 
            else 
            {   // Static derivatives

                // Shift the derivative start index due to symmetry of static derivatives
                basis_j_start_index = basis_i;
            }                

            // Derivative wrt basis_j
            for (std::size_t basis_j = basis_j_start_index; basis_j < num_eigenvalues; basis_j++)
            {
                // Set the BASIS_J counter for using in the scheme
                r_model_part.GetProcessInfo()[BASIS_J] = basis_j;

                // Reset RHS and solution vector at each step
                TSparseSpace::SetToZero(rb);
                TSparseSpace::SetToZero(rDx);

                // Compute RHS and solve
                if (mDerivativeTypeFlag)
                {   // Dynamic derivatives

                    const double start_build_rhs = OpenMPUtils::GetCurrentTime();
                    Timer::Start("BuildRHS");
                    
                    // Build first stiffness contribution
                    this->pGetBuilderAndSolver()->BuildRHS(p_scheme, r_model_part, rb);

                    // Compute RHS for dynamic derivative
                    rb -= prec_inner_prod(basis, rb)* prec_prod(rMassMatrix, basis);

                    Timer::Stop("BuildRHS");
                    const double stop_build_rhs = OpenMPUtils::GetCurrentTime();
                    KRATOS_INFO_IF("ModalDerivativeStrategy", (this->GetEchoLevel() >=1 && r_model_part.GetCommunicator().MyPID() == 0)) << "Build time RHS: " << stop_build_rhs - start_build_rhs << std::endl;

                    // Builder and Solver routines
                    if(r_model_part.MasterSlaveConstraints().size() != 0) {
                        this->pGetBuilderAndSolver()->ApplyRHSConstraints(p_scheme, r_model_part, rb);
                    }

                    Timer::Start("ApplyRHSConstraints");

                    // Apply dynamic derivative constraint
                    std::size_t fixed_dof_index;
                    this->ApplyDynamicDerivativeConstraint(basis, fixed_dof_index);

                    Timer::Stop("ApplyRHSConstraints");

                    // Apply Dirichlet conditions
                    this->pGetBuilderAndSolver()->ApplyDirichletConditions(p_scheme, r_model_part, rA, rDx, rb);

                    const double start_solve = OpenMPUtils::GetCurrentTime();
                    Timer::Start("Solve");

                    // Compute null space solution
                    this->pGetBuilderAndSolver()->SystemSolve(rA, rDx, rb);
                    // this->pGetBuilderAndSolver()->SystemSolveWithPhysics(rA, rDx, rb, r_model_part);

                    // Compute particular solution
                    this->ComputeParticularSolution(rDx, basis, fixed_dof_index);

                    Timer::Stop("Solve");
                    const double stop_solve = OpenMPUtils::GetCurrentTime();

                    KRATOS_INFO_IF("ModalDerivativeStrategy", (this->GetEchoLevel() >=1 && r_model_part.GetCommunicator().MyPID() == 0)) << "System solve time: " << stop_solve - start_solve << std::endl;
                } 
                else 
                {   // Static derivatives

                    // Build only stiffness contribution
                    this->pGetBuilderAndSolver()->BuildRHSAndSolve(p_scheme, r_model_part, rA, rDx, rb);
                }
                
                // Mass orthonormalization
                if (mMassOrthonormalizeFlag)
                    this->MassOrthonormalize(rDx);

                // Assign solution to ROM_BASIS
                this->AssignVariables(rDx);

                // Update the derivative index
                r_model_part.GetProcessInfo()[DERIVATIVE_INDEX] += 1;
            }

            // Reset all the dofs to initial value
            this->ResetVariables();
        }

        // this->ComputeReducedMatrices();

        return true;
        KRATOS_CATCH("")
    }

    /**
     * @brief This function stores the initial SolutionStepValue.
     * @details
     * { 
     * If an a priori analysis (e.g. a nonlinear analysis) is performed 
     * the last equilibrium state is stored as the equilibrium state and 
     * a starting point to compute modal derivatives,
     * at which the linearization of the nonlinear problem is performed.
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

    /**
     * @brief This function applies mass orthonormalization of the given vector
     * @details
     * { 
     * Given vector is mass-orthonormalized w.r.t. the existing and/or previously computed 
     * ROM basis vectors using the modified Gram-Schmidt method
     * } 
     */  
    void MassOrthonormalize(TSystemVectorType& rDx)
    {
        ModelPart& r_model_part = BaseType::GetModelPart();
        TSystemMatrixType& rMassMatrix = *mpMassMatrix;
        TSystemVectorType basis;
        TSparseSpace::Resize(basis, mpA->size1());
        std::size_t current_basis_index = r_model_part.GetProcessInfo()[DERIVATIVE_INDEX];

        // Mass-Orthogonalization using modified Gram-Schmidt method
        for (std::size_t basis_index = 0; basis_index < current_basis_index; basis_index++){

            // Get previous basis
            this->GetBasis(basis_index, basis);

            // Apply Mass-Orthogonalization w.r.t. previous bases
            rDx -= (prec_inner_prod(prec_prod(rMassMatrix, basis), rDx) * basis);

        }

        // Mass-Normalization
        rDx /= std::sqrt(prec_inner_prod(prec_prod(rMassMatrix, rDx), rDx));
        
    }

    /**
     * @brief This function retrieves the basis of given index
     * @details
     * { 
     * } 
     */  
    void GetBasis(std::size_t basis_index, TSystemVectorType& rBasis)
    {
        ModelPart& r_model_part = BaseType::GetModelPart();

        for (ModelPart::NodeIterator itNode = r_model_part.NodesBegin(); itNode!= r_model_part.NodesEnd(); itNode++) {
            ModelPart::NodeType::DofsContainerType& NodeDofs = itNode->GetDofs();
            const std::size_t NumNodeDofs = NodeDofs.size();
            Matrix& rRomBasis = itNode->GetValue(ROM_BASIS);

            // fill the basis vector
            auto itDof = std::begin(NodeDofs);
            for (std::size_t iDOF = 0; iDOF < NumNodeDofs; iDOF++)
            {
                rBasis((*itDof)->EquationId()) = rRomBasis(iDOF,basis_index);
                itDof++;
            }
        }
    }

    /**
     * @brief This function applies the dynamic derivative constraint
     * @details
     * { 
     * The dynamic derivative LHS is singular since A_ij = K - lambda_i*M ,
     * thus the dynamic derivative constraint has to be applied with; Basis_i^T @ M @ dBasis_i_dalpha_j = 0
     * This function applies a Dirichlet constraint on the DOF with the maximum absolute value
     * } 
     */  
    void ApplyDynamicDerivativeConstraint(TSystemVectorType& rBasis, std::size_t& fixed_dof_index)
    {
        KRATOS_TRY

        DofsArrayType& r_dof_set = this->pGetBuilderAndSolver()->GetDofSet();

        // Find the DOF with the maximum absolute value in the considered basis
        double max_abs_value = 0.0;
        double temp_abs_value = 0.0;
        for (auto dof_i : r_dof_set)
        {
            temp_abs_value = std::abs(rBasis[dof_i.EquationId()]);
            if (temp_abs_value > max_abs_value)
            {
                max_abs_value = temp_abs_value;
                fixed_dof_index = dof_i.EquationId();
            }
        }

        // Fix the found DOF
        auto p_dof = r_dof_set.begin()+fixed_dof_index;
        p_dof->FixDof();

        // // Below is a self manifactured implementation which does not necessarily give optimum results for some cases
        // // It is still kept here for future correspondence
        // TSystemMatrixType& rA = *mpA;
        // TSystemMatrixType& rMassMatrix = *mpMassMatrix;

        // DofsArrayType& r_dof_set = this->pGetBuilderAndSolver()->GetDofSet();

        // TSystemVectorType constraint;
        // constraint = prec_prod(rMassMatrix, rBasis);

        // for (auto dof_i : r_dof_set)
        // {
        //     // If the DOF is free
        //     if (dof_i.IsFree())
        //     {
        //         // then add the constraint only to the first free DOF's row and column
        //         for (auto dof_j : r_dof_set)
        //         {
        //             // Set row
        //             rA(dof_i.EquationId(),dof_j.EquationId()) += constraint[dof_j.EquationId()];
        //             // Set column
        //             rA(dof_j.EquationId(),dof_i.EquationId()) += constraint[dof_j.EquationId()];
        //         }
        //         return;
        //     }
        // }
        
        KRATOS_CATCH("")
    }

    /**
     * @brief This function comnputes the particular solution from the null space solution
     * @details
     * { 
     * } 
     */ 
    void ComputeParticularSolution(TSystemVectorType& rDx, TSystemVectorType& rBasis, std::size_t& rFixedDofIndex)
    {
        KRATOS_TRY

        TSystemMatrixType& rMassMatrix = *mpMassMatrix;
        DofsArrayType& r_dof_set = this->pGetBuilderAndSolver()->GetDofSet();

        // Free the DOF that is previously fixed
        auto p_dof = r_dof_set.begin()+rFixedDofIndex;
        p_dof->FreeDof();

        // Component c for the null space solution
        double c = prec_inner_prod(rDx, prec_prod(rMassMatrix, rBasis));

        // Particular solution
        rDx += c*rBasis;

        KRATOS_CATCH("")
    }

    /**
     * @brief This function reduces the mass and stiffness matrices using the computed basis
     * @details
     * { 
     * This is implemented purely for debugging purposes
     * } 
     */  
    void ComputeReducedMatrices()
    {
        KRATOS_TRY

        ModelPart& r_model_part = BaseType::GetModelPart();
        std::size_t derivative_index = r_model_part.GetProcessInfo()[DERIVATIVE_INDEX];

        Matrix Mred;
        Mred.resize(derivative_index,derivative_index, false);
        Matrix Kred;
        Kred.resize(derivative_index,derivative_index, false);
        TSystemVectorType basis_i;
        TSparseSpace::Resize(basis_i, mpA->size1());
        TSystemVectorType basis_j;
        TSparseSpace::Resize(basis_j, mpA->size1());
        TSystemMatrixType& rA = *mpA;
        TSystemMatrixType& rStiffnessMatrix = *mpStiffnessMatrix;
        TSystemMatrixType& rMassMatrix = *mpMassMatrix;

        for (std::size_t i = 0; i < derivative_index; i++)
        {
            this->GetBasis(i, basis_i);
            for (std::size_t j = 0; j < derivative_index; j++)
            {
                this->GetBasis(j, basis_j);
                Mred(i,j) = prec_inner_prod(basis_i, prec_prod(rMassMatrix, basis_j));
                if (mDerivativeTypeFlag)
                    Kred(i,j) = prec_inner_prod(basis_i, prec_prod(rStiffnessMatrix, basis_j));
                else
                    Kred(i,j) = prec_inner_prod(basis_i, prec_prod(rA, basis_j));
            }
        }

        KRATOS_WATCH(Mred)
        KRATOS_WATCH(Kred)

        KRATOS_CATCH("")
    }

    /**
     * @brief This function assigns the solution vector to ROM_BASIS nodal variable
     * @details
     * { 
     * } 
     */  
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

    bool mDerivativeTypeFlag;

    bool mMassOrthonormalizeFlag;

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
