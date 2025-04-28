//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//                   Riccardo Rossi
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "utilities/builtin_timer.h"

#ifdef KRATOS_USE_FUTURE
#include "future/linear_solvers/linear_solver.h"
#include "future/solving_strategies/schemes/implicit_scheme.h"
#include "future/solving_strategies/strategies/strategy.h"
#endif

namespace Kratos::Future
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name Enum's
///@{

///@}
///@name Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class ImplicitStrategy
 * @ingroup KratosCore
 * @brief This is the base class for all implicit strategies
 * @author Ruben Zorrilla
 * @author Riccardo Rossi
 */
template <class TMatrixType, class TVectorType, class TSparseGraphType>
class ImplicitStrategy : public Strategy
{
public:
    ///@name Type Definitions */
    ///@{

    /** Counted pointer of ClassName */
    KRATOS_CLASS_POINTER_DEFINITION(ImplicitStrategy);

    // Scheme pointer type definition
    using SchemePointerType = typename ImplicitScheme<TMatrixType, TVectorType, TSparseGraphType>::Pointer;

    // Linear solver pointer definition
    using LinearSolverPointerType = typename Future::LinearSolver<TMatrixType, TVectorType>::Pointer;

    /// Index type definition
    using IndexType = typename TMatrixType::IndexType;

    /// Data type definition
    using DataType = typename TMatrixType::DataType;

    /// DOFs array type definition
    using DofsArrayType = typename ModelPart::DofsArrayType;

    /// Effective DOFs map type definition
    using EffectiveDofsMapType = std::unordered_map<typename Dof<DataType>::Pointer, IndexType>;

    /// Matrix type definition
    using SystemMatrixType = TMatrixType;

    /// Matrix pointer type definition
    using SystemMatrixPointerType = typename TMatrixType::Pointer;

    /// Vector type definition
    using SystemVectorType = TVectorType;

    /// Vector type definition
    using SystemVectorPointerType = typename SystemVectorType::Pointer;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     */
    explicit ImplicitStrategy() = default;

    /**
     * @brief Default constructor. (with parameters)
     * @param rModelPart The model part of the problem
     * @param ThisParameters The configuration parameters
     */
    explicit ImplicitStrategy(
        ModelPart& rModelPart,
        Parameters ThisParameters)
        : Strategy(rModelPart, ThisParameters)
    {
        // Validate and assign defaults
        ThisParameters.ValidateAndAssignDefaults(GetDefaultParameters());
        this->AssignSettings(ThisParameters);

        //TODO: In here we should leverage the Registry to construct the scheme and linear solver from the json input settings
    }

    /**
     * @brief Default constructor
     * @param rModelPart The model part of the problem
     * @param pScheme The integration scheme
     * @param pNewLinearSolver The linear solver employed
     * @param CalculateReactionFlag The flag for the reaction calculation
     * @param ReformDofSetAtEachStep The flag that allows to compute the modification of the DOF
     * @param CalculateNormDxFlag The flag sets if the norm of Dx is computed
     * @param MoveMeshFlag The flag that allows to move the mesh
     */
    explicit ImplicitStrategy(
        ModelPart& rModelPart,
        SchemePointerType pScheme,
        LinearSolverPointerType pLinearSolver,
        bool ComputeReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool CalculateNormDxFlag = false,
        bool MoveMeshFlag = false)
        : Strategy(rModelPart)
        , mpScheme(pScheme)
        , mpLinearSolver(pLinearSolver)
        , mReformDofsAtEachStep(ReformDofSetAtEachStep)
        , mComputeReactions(ComputeReactions)
    {
        KRATOS_TRY

        // Set EchoLevel to the default value (only time is displayed) //TODO: Get this from input settings
        this->SetEchoLevel(1);

        KRATOS_CATCH("")
    }

    /** Copy constructor.
     */
    ImplicitStrategy(const ImplicitStrategy &Other) = delete;

    /**
     * @brief Destructor.
     * @details In trilinos third party library, the linear solver's preconditioner should be freed before the system matrix. We control the deallocation order with Clear().
     */
    virtual ~ImplicitStrategy()
    {
        // Clear current class
        this->Clear();
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    typename Strategy::Pointer Create(
        ModelPart &rModelPart,
        Parameters ThisParameters) const override
    {
        return Kratos::make_shared<ImplicitStrategy<TMatrixType, TVectorType, TSparseGraphType>>(rModelPart, ThisParameters);
    }

    void Initialize() override
    {
        KRATOS_TRY

        // Initialize scheme (this has to be done once)
        pGetScheme()->Initialize();

        KRATOS_CATCH("")
    }

    void InitializeSolutionStep() override
    {
        KRATOS_TRY

        // Call the scheme InitializeSolutionStep
        pGetScheme()->InitializeSolutionStep(mDofSet, mEffectiveDofMap, mpA, mpdx, mpb, mConstraintsRelationMatrix, mConstraintsConstantVector, mReformDofsAtEachStep);

        KRATOS_CATCH("")
    }

    void Predict() override
    {
        KRATOS_TRY

        // Call the time scheme predict (note that this also updates the mesh if needed)
        pGetScheme()->Predict(*mpA, *mpdx, *mpb);

        //TODO: Constraints implementation
        // // Applying constraints if needed
        // const auto &r_comm = GetModelPart().GetCommunicator().GetDataCommunicator();
        // auto& r_constraints_array = BaseType::GetModelPart().MasterSlaveConstraints();
        // const int local_number_of_constraints = r_constraints_array.size();
        // const int global_number_of_constraints = r_comm.SumAll(local_number_of_constraints);
        // if(global_number_of_constraints != 0) {
        //     const auto& r_process_info = BaseType::GetModelPart().GetProcessInfo();

        //     block_for_each(r_constraints_array, [&r_process_info](MasterSlaveConstraint& rConstraint){
        //         rConstraint.ResetSlaveDofs(r_process_info);
        //     });
        //     block_for_each(r_constraints_array, [&r_process_info](MasterSlaveConstraint& rConstraint){
        //         rConstraint.Apply(r_process_info);
        //     });

        //     //the following is needed since we need to eventually compute time derivatives after applying
        //     //Master slave relations
        //     TSparseSpace::SetToZero(rDx);
        //     this->pGetScheme()->Update(BaseType::GetModelPart(), r_dof_set, rA, rDx, rb);
        // }

        KRATOS_CATCH("")
    }

    bool SolveSolutionStep() override
    {
        KRATOS_ERROR << "\'ImplicitStategy\' does not implement \'SolveSolutionStep\'. Please call derived class." << std::endl;
    }

    void FinalizeSolutionStep() override
    {
        KRATOS_TRY;

        // Finalisation of the solution step (operations to be done after achieving convergence)
        pGetScheme()->FinalizeSolutionStep(*mpA, *mpdx, *mpb);

        // Reset flags for next step
        mStiffnessMatrixIsBuilt = false; // By default we always rebuilt (if not implement a derived strategy)

        // Clear if needed (note that this deallocates the system arrays)
        if (mReformDofsAtEachStep) {
            Clear();
        }

        KRATOS_CATCH("");
    }

    void Clear() override
    {
        KRATOS_TRY;

        // Clear the linear solver
        if (mpLinearSolver != nullptr) {
            mpLinearSolver->Clear();
        }

        // Clearing the system of equations
        if (mpA != nullptr) {
            mpA->Clear();
        }
        if (mpdx != nullptr) {
            mpdx->Clear();
        }
        if (mpb != nullptr) {
            mpb->Clear();
        }
        if (mpEffectiveRhs != nullptr) {
            mpEffectiveRhs->Clear();
        }
        if (mpEffectiveLhs != nullptr) {
            mpEffectiveLhs->Clear();
        }
        if (mpEffectiveDx != nullptr) {
            mpEffectiveDx->Clear();
        }

        // Clear the constraint-related arrays
        mConstraintsRelationMatrix.Clear();
        mConstraintsConstantVector.Clear();

        // Clearing scheme
        // Note that this resets the DOF set
        if (pGetScheme() != nullptr) {
            pGetScheme()->Clear();
        }

        KRATOS_CATCH("");
    }

    int Check() override
    {
        KRATOS_TRY

        pGetScheme()->Check();

        //TODO: think if we do this in debug mode only
        //TODO: this should probably go to the Scheme check (based on previous comment)
        // If reactions are to be calculated, we check if all the dofs have reactions defined
        if (mComputeReactions) {
            for (auto& r_dof : mDofSet) {
                KRATOS_ERROR_IF_NOT(r_dof.HasReaction())
                    << "Reaction variable not set for the following: " << std::endl
                    << "- Node: "<< r_dof.Id() << std::endl
                    << "- DOF: " << r_dof << std::endl;
            }
        }

        return 0;

        KRATOS_CATCH("")
    }

    void CalculateOutputData() override
    {
        pGetScheme()->CalculateOutputData(*mpA, *mpdx, *mpb);
    }

    Parameters GetDefaultParameters() const override
    {
        // Current class default parameters
        Parameters default_parameters = Parameters(R"({
            "name" : "implicit_strategy",
            "echo_level" : 0,
            "move_mesh" : false,
            "compute_reactions" : false,
            "reform_dofs_at_each_step" : false,
            "linear_solver_settings" : {},
            "scheme_settings" : {}
        })");

        return default_parameters;
    }

    /**
     * @brief Returns the name of the class as used in the settings (snake_case format)
     * @return The name of the class
     */
    static std::string Name()
    {
        return "linear_strategy";
    }

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief It sets the level of echo for the solving strategy
     * @param Level The level to set
     * @details The different levels of echo are:
     * - 0: Mute... no echo at all
     * - 1: Printing time and basic information
     * - 2: Printing linear solver data
     * - 3: Print of debug information: Echo of stiffness matrix, Dx, b...
     */
     void SetEchoLevel(const int EchoLevel) //TODO: I'd also unify the echo level with the scheme
     {
         mEchoLevel = EchoLevel;
         pGetScheme()->SetEchoLevel(EchoLevel);
     }

    /**
     * @brief This method sets the flag mComputeReactions
     * @param ComputeReactions The flag that tells if the reactions are computed
     */
    void SetComputeReactions(const bool ComputeReactions)
    {
        mComputeReactions = ComputeReactions;
    }

    /**
     * @brief This method sets the flag mReformDofsAtEachStep
     * @param ReformDofsAtEachStep The flag that indicates if the DOF set is reformed at each step
     */
    void SetReformDofsAtEachStep(const bool ReformDofsAtEachStep)
    {
        mReformDofsAtEachStep = ReformDofsAtEachStep;
    }

    /**
     * @brief This method sets the flag mStiffnessMatrixIsBuilt
     * @param StiffnessMatrixIsBuilt The flag that indicates if the stiffness matrix is built
     */
    void SetStiffnessMatrixIsBuilt(const bool StiffnessMatrixIsBuilt)
    {
        mStiffnessMatrixIsBuilt = StiffnessMatrixIsBuilt;
    }

    /**
     * @brief Set method for the time scheme
     * @param pScheme The pointer to the time scheme considered
     */
    void SetScheme(SchemePointerType pScheme)
    {
        mpScheme = pScheme;
    }

    /**
     * @brief Get method for the time scheme
     * @return mpScheme: The pointer to the time scheme considered
     */
    SchemePointerType pGetScheme()
    {
        return mpScheme;
    }

    /**
     * @brief Set method for the linear solver
     * @param pLinearSolver The pointer to the linear solver
     */
    void SetLinearSolver(LinearSolverPointerType pLinearSolver)
    {
        mpLinearSolver = pLinearSolver;
    }

    /**
     * @brief Get method for the linear solver
     * @return mpLinearSolver The pointer to the linear solver
     */
    LinearSolverPointerType pGetLinearSolver()
    {
        return mpLinearSolver;
    }

    /**
     * @brief This method returns the LHS matrix
     * @return The LHS matrix
     */
    SystemMatrixType& GetSystemMatrix()
    {
        return *mpA;
    }

    /**
     * @brief This method returns the LHS matrix
     * @return The LHS matrix
     */
    typename SystemMatrixType::Pointer pGetSystemMatrix()
    {
        return mpA;
    }

    /**
     * @brief This method returns the effective LHS matrix
     * @return The effective LHS matrix
     */
    SystemMatrixPointerType pGetEffectiveSystemMatrix()
    {
        return mpEffectiveLhs;
    }

    /**
     * @brief This method returns the RHS vector
     * @return The RHS vector
     */
    SystemVectorType& GetSystemVector()
    {
        return *mpb;
    }

    /**
     * @brief This method returns the RHS vector
     * @return The RHS vector
     */
    typename SystemVectorType::Pointer pGetSystemVector()
    {
        return mpb;
    }

    /**
     * @brief This method returns the effective RHS vector
     * @return The effective RHS vector
     */
    typename TVectorType::Pointer pGetEffectiveSystemVector()
    {
        return mpEffectiveRhs;
    }

    /**
     * @brief This method returns the constraints constant vector (i.e., b)
     * @return Pointer to the constaints constant vector
     */
    SystemVectorType& GetConstraintsConstantVector()
    {
        return mConstraintsConstantVector;
    }

    /**
     * @brief This method returns the constraints relation matrix (i.e., T)
     * @return Pointer to the constraints relation matrix
     */
    SystemMatrixType& GetConstraintsRelationMatrix()
    {
        return mConstraintsRelationMatrix;
    }

    /**
     * @brief This method returns the solution vector
     * @return The Dx vector
     */
    SystemVectorType& GetSolutionVector()
    {
        return *mpdx;
    }

    /**
     * @brief This method returns the solution vector
     * @return The Dx vector
     */
    typename SystemVectorType::Pointer pGetSolutionVector()
    {
        return mpdx;
    }

    /**
     * @brief This method returns the effective solution vector
     * @return The effective Dx vector
     */
    typename SystemVectorType::Pointer pGetEffectiveSolutionVector()
    {
        return mpEffectiveDx;
    }

    /**
     * @brief It allows to get the list of Dofs from the element
     */
    DofsArrayType& GetDofSet()
    {
        return mDofSet;
    }

    /**
     * @brief It allows to get the list of Dofs from the element
     */
    const DofsArrayType& GetDofSet() const
    {
        return mDofSet;
    }

    /**
     * @brief It allows to get the map of effective DOFs
     */
    EffectiveDofsMapType& GetEffectiveDofMap()
    {
        return mEffectiveDofMap;
    }

    /**
     * @brief It allows to get the map of effective DOFs
     */
    const EffectiveDofsMapType &GetEffectiveDofMap() const
    {
        return mEffectiveDofMap;
    }

    ///@}
    ///@name Inquiry
    ///@{

    /**
     * @brief This returns the level of echo for the solving strategy
     * @details
     * {
     * 0 -> Mute... no echo at all
     * 1 -> Printing time and basic information
     * 2 -> Printing linear solver data
     * 3 -> Print of debug information: Echo of stiffness matrix, Dx, b...
     * }
     * @return Level of echo for the solving strategy
     */
    int GetEchoLevel() const
    {
        return mEchoLevel;
    }

    /**
     * @brief This method returns the flag mComputeReactions
     * @return The flag that tells if the reactions are computed
     */
    bool GetComputeReactions() const
    {
        return mComputeReactions;
    }

    /**
     * @brief This method returns the flag mReformDofsAtEachStep
     * @return The flag that indicates if the DOF set is reformed at each step
     */
    bool GetReformDofsAtEachStep() const
    {
        return mReformDofsAtEachStep;
    }

    /**
     * @brief This method returns the flag mStiffnessMatrixIsBuilt
     * @return The flag that indicates if the stiffness matrix is built
     */
    bool GetStiffnessMatrixIsBuilt() const
    {
        return mStiffnessMatrixIsBuilt;
    }

    /**
     * @brief This method returns the residual norm
     * @return The residual norm
     */
    double GetResidualNorm()
    {
        if (mpb->size() != 0) {
            return mpb->Norm();
        } else {
            return 0.0;
        }
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ImplicitStrategy";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

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

    LinearSolverPointerType mpLinearSystemSolver = nullptr; /// Pointer to the linear solver

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief This method assigns settings to member variables
     * @param ThisParameters Parameters that are assigned to the member variables
     */
    void AssignSettings(const Parameters ThisParameters)
    {
        // BaseType::AssignSettings(ThisParameters); // TODO: Once we have a base class
        mEchoLevel = ThisParameters["echo_level"].GetInt();
        mComputeReactions = ThisParameters["compute_reactions"].GetBool();
        mReformDofsAtEachStep = ThisParameters["reform_dofs_at_each_step"].GetBool();

        // Saving the scheme
        if (ThisParameters["scheme_settings"].Has("name")) {
            KRATOS_ERROR << "IMPLEMENTATION PENDING IN CONSTRUCTOR WITH PARAMETERS" << std::endl;
        }
    }

    /**
     * @brief This method returns the components of the system of equations depending of the echo level
     */
    virtual void EchoInfo()
    {
        const auto& r_A = *mpA;
        const auto& r_b = *mpb;
        const auto& r_dx = *mpdx;

        if (GetEchoLevel() == 3) { //if it is needed to print the debug info
            KRATOS_INFO("LHS") << "LHS = " << r_A << std::endl;
            KRATOS_INFO("RHS") << "RHS = " << r_b << std::endl;
            KRATOS_INFO("Dx")  << "Solution obtained = " << r_dx << std::endl;
        }

        if (this->GetEchoLevel() == 4) { //print to matrix market file
            //FIXME: Impllement the MatrixMarket output

            // std::stringstream matrix_market_name;
            // matrix_market_name << "A_" << BaseType::GetModelPart().GetProcessInfo()[TIME] <<  ".mm";
            // TSparseSpace::WriteMatrixMarketMatrix((char*) (matrix_market_name.str()).c_str(), rA, false);

            // std::stringstream matrix_market_vectname;
            // matrix_market_vectname << "b_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << ".mm.rhs";
            // TSparseSpace::WriteMatrixMarketVector((char*) (matrix_market_vectname.str()).c_str(), rb);
        }
    }

    ///@}
    ///@name Protected  Access
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

    ModelPart* mpModelPart = nullptr;

    SchemePointerType mpScheme = nullptr; /// The pointer to the scheme

    LinearSolverPointerType mpLinearSolver = nullptr; /// The pointer to the linear solver

    DofsArrayType mDofSet; /// The set containing the DOFs of the system

    EffectiveDofsMapType mEffectiveDofMap; /// The map containing the effective DOFs of the system

    //TODO: Should these be unique_ptr?
    SystemVectorPointerType mpdx = nullptr; /// The incremement in the solution //TODO: use naming convention

    //TODO: Should these be unique_ptr?
    SystemVectorPointerType mpb = nullptr; /// The RHS vector of the system of equations //TODO: use naming convention (mpRHS)

    //TODO: Should these be unique_ptr?
    SystemMatrixPointerType mpA = nullptr; /// The LHS matrix of the system of equations //TODO: use naming convention (mpLHS)

    //TODO: Should these be unique_ptr? --> This we cannot for the case without constraints (effective and original point to the same matrix)
    SystemMatrixPointerType mpEffectiveLhs = nullptr; /// The LHS matrix of the system of equations

    //TODO: Should these be unique_ptr? --> This we cannot for the case without constraints (effective and original point to the same vector)
    SystemVectorPointerType mpEffectiveRhs = nullptr; /// The RHS vector of the system of equations

    //TODO: Should these be unique_ptr? --> This we cannot for the case without constraints (effective and original point to the same vector)
    SystemVectorPointerType mpEffectiveDx = nullptr; /// The effective solution increment for the constrained system

    SystemMatrixType mConstraintsRelationMatrix; // Constraints relation matrix

    SystemVectorType mConstraintsConstantVector; // Constraints constant vector

    int mEchoLevel;

    bool mReformDofsAtEachStep = false;

    bool mStiffnessMatrixIsBuilt = false;

    bool mComputeReactions; // Flag to activate the reactions calculation from the residual vector

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

    ///@}
}; /* Class ImplicitStrategy */

///@}

///@name Type Definitions */
///@{


///@}

} // namespace Kratos::Future

