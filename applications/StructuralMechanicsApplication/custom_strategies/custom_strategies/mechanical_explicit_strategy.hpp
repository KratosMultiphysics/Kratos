//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license:
//kratos/license.txt
//
//  Main authors:    Klaus B Sautter (based on the work of JMCarbonel)
//
//

#if !defined(KRATOS_MECHANICAL_EXPLICIT_STRATEGY)
#define KRATOS_MECHANICAL_EXPLICIT_STRATEGY

/* System includes */

/* Project includes */
#include "solving_strategies/strategies/solving_strategy.h"
#include "structural_mechanics_application_variables.h"
#include "utilities/variable_utils.h"

namespace Kratos {
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
 * @class MechanicalExplicitStrategy
 * @ingroup StructuralMechanicsApplciation
 * @brief This strategy is used for the explicit time integration
 * @author Klauss B Sautter (based on the work of JMCarbonel)
 */
template <class TSparseSpace,
          class TDenseSpace,  // = DenseSpace<double>,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
class MechanicalExplicitStrategy
    : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> {
public:
    ///@name Type Definitions
    ///@{

    // Base class definition
    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    /// Some definitions from the base class
    typedef typename BaseType::TSchemeType TSchemeType;
    typedef typename BaseType::DofsArrayType DofsArrayType;
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;
    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;
    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;
    typedef typename BaseType::NodesArrayType NodesArrayType;
    typedef typename BaseType::ElementsArrayType ElementsArrayType;
    typedef typename BaseType::ConditionsArrayType ConditionsArrayType;
    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    /// Counted pointer of MechanicalExplicitStrategy
    KRATOS_CLASS_POINTER_DEFINITION(MechanicalExplicitStrategy);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Default constructor
     * @param rModelPart The model part of the problem
     * @param pScheme The integration scheme
     * @param CalculateReactions The flag for the reaction calculation
     * @param ReformDofSetAtEachStep The flag that allows to compute the modification of the DOF
     * @param MoveMeshFlag The flag that allows to move the mesh
     */
    MechanicalExplicitStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = true)
        : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, MoveMeshFlag) {
        KRATOS_TRY

        // set flags to default values
        this->mCalculateReactionsFlag = CalculateReactions;
        this->mReformDofSetAtEachStep = ReformDofSetAtEachStep;

        // saving the scheme
        this->mpScheme = pScheme;

        // set EchoLevel to the default value (only time is displayed)
        BaseType::SetEchoLevel(1);

        // set RebuildLevel to the default value
        BaseType::SetRebuildLevel(0);

        KRATOS_CATCH("")
    }

    /** Destructor.
    */
    virtual ~MechanicalExplicitStrategy()
    {
        Clear();
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

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

    // Set and Get Flags

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
     * @brief This method sets the flag mCalculateReactionsFlag
     * @param CalculateReactionsFlag The flag that tells if the reactions are computed
     */
    void SetCalculateReactionsFlag(bool CalculateReactionsFlag)
    {
        mCalculateReactionsFlag = CalculateReactionsFlag;
    }

    /**
     * @brief This method returns the flag mCalculateReactionsFlag
     * @return The flag that tells if the reactions are computed
     */
    bool GetCalculateReactionsFlag()
    {
        return mCalculateReactionsFlag;
    }

    /**
     * @brief This method sets the flag mReformDofSetAtEachStep
     * @param Flag The flag that tells if each time step the system is rebuilt
     */
    void SetReformDofSetAtEachStepFlag(bool Flag)
    {
        mReformDofSetAtEachStep = Flag;
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
     * @brief This function is designed to move the mesh
     * @note Be careful it just consider displacements, derive this method to adapt to your own strategies (ALE, FSI, etc...)
     */
    void MoveMesh() override
    {
        KRATOS_TRY

        KRATOS_ERROR_IF(BaseType::GetModelPart().NodesBegin()->SolutionStepsDataHas(DISPLACEMENT_X) == false) << "It is impossible to move the mesh since the DISPLACEMENT var is not in the Model Part. Either use SetMoveMeshFlag(False) or add DISPLACEMENT to the list of variables" << std::endl;

        NodesArrayType& r_nodes_array = BaseType::GetModelPart().Nodes();
        const auto it_node_begin = r_nodes_array.begin();

        #pragma omp parallel for
        for(int i = 0; i < static_cast<int>(r_nodes_array.size()); ++i) {
            auto it_node = it_node_begin + i;

            noalias(it_node->Coordinates()) = it_node->GetInitialPosition().Coordinates();
            noalias(it_node->Coordinates()) += it_node->FastGetSolutionStepValue(DISPLACEMENT);
        }

        KRATOS_INFO_IF("SolvingStrategy", this->GetEchoLevel() > 1 && BaseType::GetModelPart().GetCommunicator().MyPID() == 0) << " MESH MOVED " << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * @brief Initialization of member variables and prior operations
     */
    void Initialize() override
    {
        KRATOS_TRY

        if (!this->mInitializeWasPerformed){
            // Pointers needed in the solution
            typename TSchemeType::Pointer pScheme = GetScheme();
            ModelPart& r_model_part = BaseType::GetModelPart();

            TSystemMatrixType matrix_a_dummy = TSystemMatrixType();

            // Initialize The Scheme - OPERATIONS TO BE DONE ONCE
            if (!pScheme->SchemeIsInitialized())pScheme->Initialize(r_model_part);

            // Initialize The Elements - OPERATIONS TO BE DONE ONCE
            if (!pScheme->ElementsAreInitialized())pScheme->InitializeElements(r_model_part);

            // Initialize The Conditions- OPERATIONS TO BE DONE ONCE
            if (!pScheme->ConditionsAreInitialized())pScheme->InitializeConditions(r_model_part);

            // Set Nodal Mass to zero
            NodesArrayType& r_nodes = r_model_part.Nodes();
            VariableUtils().SetNonHistoricalVariable(NODAL_MASS, 0.0, r_nodes);

            // Iterate over the elements
            ElementsArrayType& r_elements = r_model_part.Elements();
            const auto it_elem_begin = r_elements.begin();
            ProcessInfo& r_current_process_info = r_model_part.GetProcessInfo();

            Vector dummy_vector;
            // If we consider the rotation DoF
            const bool has_dof_for_rot_z = r_model_part.Nodes().begin()->HasDofFor(ROTATION_Z);
            if (has_dof_for_rot_z) {
                const array_1d<double, 3> zero_array = ZeroVector(3);
                VariableUtils().SetNonHistoricalVariable(NODAL_INERTIA, zero_array, r_nodes);

                #pragma omp parallel for firstprivate(dummy_vector)
                for (int i = 0; i < static_cast<int>(r_elements.size()); ++i) {
                    // Getting nodal mass and inertia from element
                    // this function needs to be implemented in the respective
                    // element to provide inertias and nodal masses
                    auto it_elem = it_elem_begin + i;
                    it_elem->AddExplicitContribution(dummy_vector, RESIDUAL_VECTOR, NODAL_INERTIA, r_current_process_info);
                }
            } else { // Only NODAL_MASS is needed
                #pragma omp parallel for firstprivate(dummy_vector)
                for (int i = 0; i < static_cast<int>(r_elements.size()); ++i) {
                    // Getting nodal mass and inertia from element
                    // this function needs to be implemented in the respective
                    // element to provide nodal masses
                    auto it_elem = it_elem_begin + i;
                    it_elem->AddExplicitContribution(dummy_vector, RESIDUAL_VECTOR, NODAL_MASS, r_current_process_info);
                }
            }

            this->mInitializeWasPerformed = true;
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief Performs all the required operations that should be done (for each step) before solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     */
    void InitializeSolutionStep() override {
        KRATOS_TRY

        typename TSchemeType::Pointer pScheme = GetScheme();
        ModelPart& r_model_part = BaseType::GetModelPart();

        TSystemMatrixType matrix_a_dummy = TSystemMatrixType();
        TSystemVectorType mDx = TSystemVectorType();
        TSystemVectorType mb = TSystemVectorType();

        // initial operations ... things that are constant over the Solution Step
        pScheme->InitializeSolutionStep(r_model_part, matrix_a_dummy, mDx, mb);

        if (BaseType::mRebuildLevel > 0) {
            ProcessInfo& r_current_process_info = r_model_part.GetProcessInfo();
            ElementsArrayType& r_elements = r_model_part.Elements();
            const auto it_elem_begin = r_elements.begin();

            Vector dummy_vector;
            // If we consider the rotation DoF
            const bool has_dof_for_rot_z = r_model_part.Nodes().begin()->HasDofFor(ROTATION_Z);
            if (has_dof_for_rot_z) {
                #pragma omp parallel for firstprivate(dummy_vector)
                for (int i = 0; i < static_cast<int>(r_elements.size()); ++i) {
                    // Getting nodal mass and inertia from element
                    // this function needs to be implemented in the respective
                    // element to provide inertias and nodal masses
                    auto it_elem = it_elem_begin + i;
                    it_elem->AddExplicitContribution(dummy_vector, RESIDUAL_VECTOR, NODAL_INERTIA, r_current_process_info);
                }
            } else { // Only NODAL_MASS is needed
                #pragma omp parallel for firstprivate(dummy_vector)
                for (int i = 0; i < static_cast<int>(r_elements.size()); ++i) {
                    // Getting nodal mass and inertia from element
                    // this function needs to be implemented in the respective
                    // element to provide nodal masses
                    auto it_elem = it_elem_begin + i;
                    it_elem->AddExplicitContribution(dummy_vector, RESIDUAL_VECTOR, NODAL_MASS, r_current_process_info);
                }
            }
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief This method add the contributions of the residual
     * @param pScheme The integration scheme considered
     * @param rModelPart The model of the problem to solve
     */
    void CalculateAndAddRHS(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart
        )
    {
        KRATOS_TRY

        ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
        ConditionsArrayType& r_conditions = rModelPart.Conditions();
        ElementsArrayType& r_elements = rModelPart.Elements();

        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(r_conditions.size()); ++i) {
            auto it_cond = r_conditions.begin() + i;
            LocalSystemVectorType RHS_Condition_Contribution = LocalSystemVectorType(0);
            Element::EquationIdVectorType equation_id_vector_dummy; // Dummy

            pScheme->Condition_Calculate_RHS_Contribution(
                (*it_cond.base()), RHS_Condition_Contribution,
                equation_id_vector_dummy, r_current_process_info);
        }

        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(r_elements.size()); ++i) {
            auto it_elem = r_elements.begin() + i;
            LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);
            Element::EquationIdVectorType equation_id_vector_dummy; // Dummy

            pScheme->Calculate_RHS_Contribution(
                (*it_elem.base()), RHS_Contribution, equation_id_vector_dummy,
                r_current_process_info);
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief Solves the current step. This function returns true if a solution has been found, false otherwise.
     */
    bool SolveSolutionStep() override
    {
        typename TSchemeType::Pointer pScheme = GetScheme();
        ModelPart& r_model_part = BaseType::GetModelPart();
        DofsArrayType dof_set_dummy;
        TSystemMatrixType mA = TSystemMatrixType();
        TSystemVectorType mDx = TSystemVectorType();
        TSystemVectorType mb = TSystemVectorType();

        pScheme->InitializeNonLinIteration(BaseType::GetModelPart(), mA, mDx, mb);

        this->CalculateAndAddRHS(pScheme, r_model_part);

        pScheme->Update(r_model_part, dof_set_dummy, mA, mDx,
                        mb); // Explicitly integrates the equation of motion.
        return true;
    }

    /**
     * @brief Performs all the required operations that should be done (for each step) after solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     */
    void FinalizeSolutionStep() override
    {
        typename TSchemeType::Pointer pScheme = GetScheme();
        ModelPart& r_model_part = BaseType::GetModelPart();
        TSystemMatrixType mA = TSystemMatrixType();
        TSystemVectorType mDx = TSystemVectorType();
        TSystemVectorType mb = TSystemVectorType();
        // Finalisation of the solution step,
        // operations to be done after achieving convergence, for example the
        // Final Residual Vector (mb) has to be saved in there
        // to avoid error accumulation
        pScheme->FinalizeSolutionStep(r_model_part, mA, mDx, mb);

        // Move the mesh if needed
        if (BaseType::MoveMeshFlag())
            MoveMesh();

        // Cleaning memory after the solution
        pScheme->Clean();
    }

    /**
     * @brief Clears the internal storage
     */
    void Clear() override
    {
        KRATOS_TRY

        KRATOS_INFO("MechanicalExplicitStrategy") << "Clear function used" << std::endl;

        GetScheme()->Clear();
        mInitializeWasPerformed = false;

        KRATOS_CATCH("")
    }

    /**
     * @brief This function is designed to be called once to perform all the checks needed
     * on the input provided.
     * @details Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param rModelPart The model of the problem to solve
     * @return Zero means  all ok
     */
    int Check() override
    {
        KRATOS_TRY

        BaseType::Check();

        GetScheme()->Check(BaseType::GetModelPart());

        return 0;

        KRATOS_CATCH("")
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{
private:
    ///@name Static Member Variables
    ///@{
    ///@}
    ///@name Member Variables
    ///@{

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

    ///@}
    ///@name Serialization
    ///@{

    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    typename TSchemeType::Pointer mpScheme;

    /**
    Flag telling if it is needed to reform the DofSet at each
    solution step or if it is possible to form it just once
    - true  => reforme at each time step
    - false => form just one (more efficient)

    Default = false
    */
    bool mReformDofSetAtEachStep;

    /**
    Flag telling if it is needed or not to compute the reactions

    default = true
    */
    bool mCalculateReactionsFlag;

    bool mInitializeWasPerformed = false;

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

    /** Copy constructor.
    */
    MechanicalExplicitStrategy(const MechanicalExplicitStrategy& Other){};

    ///@}

}; /* Class MechanicalExplicitStrategy */

///@}

} /* namespace Kratos.*/

#endif /* KRATOS_EXPLICIT_STRATEGY  defined */
