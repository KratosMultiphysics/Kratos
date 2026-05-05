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
//
//

#pragma once

// System includes
#include <set>
#include <unordered_set>

// External includes

// Project includes
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"
#include "utilities/constraint_utilities.h"
#include "utilities/reduction_utilities.h"
#include "utilities/atomic_utilities.h"
#include "includes/kratos_parameters.h"
#include "factories/factory.h"

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
 * @class ExplicitBuilder
 * @ingroup KratosCore
 * @brief Current class provides an implementation for the base explicit builder and solving operations.
 * @details The RHS is constituted by the unbalanced loads (residual)
 * Degrees of freedom are reordered putting the restrained degrees of freedom at
 * the end of the system ordered in reverse order with respect to the DofSet.
 * @author Ruben Zorrilla
 */
template<class TSparseSpace, class TDenseSpace >
class ExplicitBuilder
{
public:
    ///@name Type Definitions
    ///@{

    /// Definition of the size type
    using SizeType = std::size_t;

    /// Definition of the index type
    using IndexType = std::size_t;

    /// Definition of the data type
    using TDataType = typename TSparseSpace::DataType;

    ///Definition of the sparse matrix
    using TSystemMatrixType = typename TSparseSpace::MatrixType;

    /// Definition of the vector size
    using TSystemVectorType = typename TSparseSpace::VectorType;

    /// Definition of the pointer to the sparse matrix
    using TSystemMatrixPointerType = typename TSparseSpace::MatrixPointerType;

    /// Definition of the pointer to the vector
    using TSystemVectorPointerType = typename TSparseSpace::VectorPointerType;

    /// The local matrix definition
    using LocalSystemMatrixType = typename TDenseSpace::MatrixType;

    /// The local vector definition
    using LocalSystemVectorType = typename TDenseSpace::VectorType;

    /// Definition of the DoF class
    using DofType = ModelPart::DofType;

    /// Definition of the DoF array type
    using DofsArrayType = ModelPart::DofsArrayType;

    /// Definition of the DoF vector type
    using DofsVectorType = ModelPart::DofsVectorType;

    /// The definition of the DoF set type
    using DofSetType = typename std::unordered_set<DofType::Pointer, DofPointerHasher>;

    /// The containers of the entities
    using NodesArrayType = ModelPart::NodesContainerType;
    using ElementsArrayType = ModelPart::ElementsContainerType;
    using ConditionsArrayType = ModelPart::ConditionsContainerType;

    /// The definition of the element container type
    using ElementsContainerType = PointerVectorSet<Element, IndexedObject>;

    /// The definition of the current class
    using ClassType = ExplicitBuilder<TSparseSpace, TDenseSpace>;

    /// Pointer definition of ExplicitBuilder
    KRATOS_CLASS_POINTER_DEFINITION(ExplicitBuilder);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Construct a new Explicit Builder object
     * Default constructor with Parameters
     * @param ThisParameters The configuration parameters
     */
    explicit ExplicitBuilder(Parameters ThisParameters)
    {
        // Validate and assign defaults
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);
    }

    /**
     * @brief Construct a new Explicit Builder object
     * Default empty constructor
     */
    explicit ExplicitBuilder() = default;

    /**
     * @brief Destroy the Explicit Builder object
     * Default destructor
     */
    virtual ~ExplicitBuilder() = default;


    /**
     * @brief Create method
     * @param ThisParameters The configuration parameters
     */
    virtual typename ClassType::Pointer Create(Parameters ThisParameters) const
    {
        return Kratos::make_shared<ClassType>(ThisParameters);
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method returns the flag mCalculateReactionsFlag
     * @return The flag that tells if the reactions are computed
     */
    bool GetCalculateReactionsFlag() const
    {
        return mCalculateReactionsFlag;
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
     * @brief This method returns the flag mDofSetIsInitialized
     * @return The flag that tells if the dof set is initialized
     */
    bool GetDofSetIsInitializedFlag() const
    {
        return mDofSetIsInitialized;
    }

    /**
     * @brief This method sets the flag mDofSetIsInitialized
     * @param DofSetIsInitialized The flag that tells if the dof set is initialized
     */
    void SetDofSetIsInitializedFlag(bool DofSetIsInitialized)
    {
        mDofSetIsInitialized = DofSetIsInitialized;
    }

    /**
     * @brief This method returns the flag mReshapeMatrixFlag
     * @return The flag that tells if we need to reset the DOF set
     */
    bool GetResetDofSetFlag() const
    {
        return mResetDofSetFlag;
    }

    /**
     * @brief This method sets the flag mResetDofSetFlag
     * @param mResetDofSetFlag The flag that tells if we need to reset the DOF set
     */
    void SetResetDofSetFlag(bool ResetDofSetFlag)
    {
        mResetDofSetFlag = ResetDofSetFlag;
    }

    /**
     * @brief This method returns the flag GetResetLumpedMassVectorFlag
     * @return The flag that tells if we need to reset the lumped mass vector
     */
    bool GetResetLumpedMassVectorFlag() const
    {
        return mResetLumpedMassVectorFlag;
    }

    /**
     * @brief This method sets the flag mResetLumpedMassVectorFlag
     * @param ResetLumpedMassVectorFlag The flag that tells if we need to reset the lumped mass vector
     */
    void SetResetLumpedMassVectorFlag(bool ResetLumpedMassVectorFlag)
    {
        mResetLumpedMassVectorFlag = ResetLumpedMassVectorFlag;
    }

    /**
     * @brief This method returns the value mEquationSystemSize
     * @return Size of the system of equations
     */
    unsigned int GetEquationSystemSize() const
    {
        return mEquationSystemSize;
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
     * @brief Get the lumped mass matrix vector pointer
     * It allows to get the lumped mass matrix vector pointer
     * @return TSystemVectorPointerType& The lumped mass matrix vector pointer
     */
    TSystemVectorPointerType& pGetLumpedMassMatrixVector()
    {
        return mpLumpedMassVector;
    }

    /**
     * @brief Get the lumped mass matrix vector
     * It allows to get the lumped mass matrix vector
     * @return TSystemVectorType& The lumped mass matrix vector
     */
    TSystemVectorType& GetLumpedMassMatrixVector()
    {
        KRATOS_ERROR_IF_NOT(mpLumpedMassVector) << "Lumped mass matrix vector is not initialized!" << std::endl;
        return (*mpLumpedMassVector);
    }

    /**
     * @brief Function to perform the build of the RHS.
     * The vector could be sized as the total number of dofs or as the number of unrestrained ones
     * @param rModelPart The model part to compute
     */
    virtual void BuildRHS(ModelPart& rModelPart)
    {
        KRATOS_TRY

        // Build the Right Hand Side without Dirichlet conditions
        // We skip setting the Dirichlet nodes residual to zero for the sake of efficiency
        // Note that this is not required as the Dirichlet conditions are set in the strategy
        BuildRHSNoDirichlet(rModelPart);

        KRATOS_CATCH("")
    }

    /**
     * @brief Function to perform the build of the RHS.
     * The vector could be sized as the total number of dofs or as the number of unrestrained ones
     * @param rModelPart The model part to compute
     */
    virtual void BuildRHSNoDirichlet(ModelPart& rModelPart)
    {
        KRATOS_TRY

        // Initialize the reaction (residual)
        InitializeDofSetReactions();

        // Gets the array of elements, conditions and constraints from the modeler
        const auto &r_elements_array = rModelPart.Elements();
        const auto &r_conditions_array = rModelPart.Conditions();
        const std::size_t n_elems = r_elements_array.size();
        const std::size_t n_conds = r_conditions_array.size();

        const auto& r_process_info = rModelPart.GetProcessInfo();

        // Assemble all elements
        IndexPartition<std::size_t>(n_elems).for_each([&](std::size_t i_elem) {
            auto it_elem = r_elements_array.begin() + i_elem;
            // If the element is active
            if (it_elem->IsActive()) {
                // Calculate elemental explicit residual contribution
                // The explicit builder and solver assumes that the residual contribution is assembled in the REACTION variables
                it_elem->AddExplicitContribution(r_process_info);
            }
        });

        // Assemble all conditions
        IndexPartition<std::size_t>(n_conds).for_each([&](std::size_t i_cond) {
            auto it_cond = r_conditions_array.begin() + i_cond;
            // If the condition is active
            if (it_cond->IsActive()) {
                // Calculate condition explicit residual contribution
                // The explicit builder and solver assumes that the residual contribution is assembled in the REACTION variables
                it_cond->AddExplicitContribution(r_process_info);
            }
        });

        KRATOS_CATCH("")
    }

    /**
     * @brief Applies the constraints
     * @param rModelPart The model part to compute
     * @param rA The LHS matrix of the system of equations
     * @param rb The RHS vector of the system of equations
     */
    virtual void ApplyConstraints(ModelPart& rModelPart)
    {
        // First we reset the slave dofs
        ConstraintUtilities::ResetSlaveDofs(rModelPart);

        // Now we apply the constraints
        ConstraintUtilities::ApplyConstraints(rModelPart);
    }

    /**
     * @brief It applied those operations that are expected to be executed once
     * @param rModelPart
     */
    virtual void Initialize(ModelPart& rModelPart)
    {
        if (!mDofSetIsInitialized) {
            // Initialize the DOF set and the equation ids
            this->SetUpDofSet(rModelPart);
            this->SetUpDofSetEquationIds();
            // Set up the lumped mass vector
            this->SetUpLumpedMassVector(rModelPart);
        } else if (!mLumpedMassVectorIsInitialized) {
            KRATOS_WARNING("ExplicitBuilder") << "Calling Initialize() with already initialized DOF set. Initializing lumped mass vector." << std::endl;;
            // Only set up the lumped mass vector
            this->SetUpLumpedMassVector(rModelPart);
        } else {
            KRATOS_WARNING("ExplicitBuilder") << "Calling Initialize() with already initialized DOF set and lumped mass vector." << std::endl;;
        }
    }

    /**
     * @brief It applies certain operations at the system of equations at the beginning of the solution step
     * @param rModelPart The model part to compute
     */
    virtual void InitializeSolutionStep(ModelPart& rModelPart)
    {
        // Check the operations that are required to be done
        if (mResetDofSetFlag) {
            // If required (e.g. topology changes) reset the DOF set
            // Note that we also set lumped mass vector in this case
            this->SetUpDofSet(rModelPart);
            this->SetUpDofSetEquationIds();
            this->SetUpLumpedMassVector(rModelPart);
        } else if (mResetLumpedMassVectorFlag) {
            // Only reset the lumped mass vector
            this->SetUpLumpedMassVector(rModelPart);
        }

        // Initialize the reactions (residual)
        this->InitializeDofSetReactions();
    }

    /**
     * @brief It applies certain operations at the system of equations at the end of the solution step
     * @param rModelPart The model part to compute
     */
    virtual void FinalizeSolutionStep(ModelPart& rModelPart)
    {
        // If required, calculate the reactions
        if (mCalculateReactionsFlag) {
            this->CalculateReactions();
        }
    }

    /**
     * @brief This function is intended to be called at the end of the solution step to clean up memory
     * storage not needed
     */
    virtual void Clear()
    {
        this->mDofSet = DofsArrayType();
        this->mpLumpedMassVector.reset();

        KRATOS_INFO_IF("ExplicitBuilder", this->GetEchoLevel() > 0) << "Clear Function called" << std::endl;
    }

    /**
     * @brief This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param rModelPart The model part to compute
     * @return 0 all ok
     */
    virtual int Check(const ModelPart& rModelPart) const
    {
        KRATOS_TRY

        return 0;

        KRATOS_CATCH("");
    }

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    virtual Parameters GetDefaultParameters() const
    {
        const Parameters default_parameters = Parameters(R"(
        {
            "name" : "explicit_builder"
        })");
        return default_parameters;
    }

    /**
     * @brief Returns the name of the class as used in the settings (snake_case format)
     * @return The name of the class
     */
    static std::string Name()
    {
        return "explicit_builder";
    }

    /**
     * @brief It sets the level of echo for the solving strategy
     * @param Level The level to set
     * @details The different levels of echo are:
     * - 0: Mute... no echo at all
     * - 1: Printing time and basic information
     * - 2: Printing linear solver data
     * - 3: Print of debug information: Echo of stiffness matrix, Dx, b...
     * - 4: Print of stiffness matrix, b to Matrix Market
     */
    void SetEchoLevel(int Level)
    {
        mEchoLevel = Level;
    }

    /**
     * @brief It returns the echo level
     * @return The echo level of the builder and solver
     */
    int GetEchoLevel() const
    {
        return mEchoLevel;
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

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "ExplicitBuilder";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
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

    DofsArrayType mDofSet; /// The set containing the DoF of the system

    TSystemVectorPointerType mpLumpedMassVector; // The lumped mass vector associated to the DOF set

    bool mResetDofSetFlag = false;  /// If the DOF set requires to be set at each time step

    bool mResetLumpedMassVectorFlag = false;  /// If the lumped mass vector requires to be set at each time step

    bool mDofSetIsInitialized = false; /// Flag taking care if the dof set was initialized ot not

    bool mLumpedMassVectorIsInitialized = false; /// Flag taking care if the lumped mass vector was initialized or not

    bool mCalculateReactionsFlag = false; /// Flag taking in account if it is needed or not to calculate the reactions

    unsigned int mEquationSystemSize; /// Number of degrees of freedom of the problem to be solve

    int mEchoLevel = 0;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief Builds the list of the DofSets involved in the problem by "asking" to each element and condition its Dofs.
     * @details The list of dofs is stores inside the ExplicitBuilder as it is closely connected to the way the matrix and RHS are built
     * @param rModelPart The model part to compute
     */
    virtual void SetUpDofSet(const ModelPart& rModelPart)
    {
        KRATOS_TRY;

        KRATOS_INFO_IF("ExplicitBuilder", this->GetEchoLevel() > 1) << "Setting up the dofs" << std::endl;

        // Gets the array of elements, conditions and constraints from the modeler
        const auto &r_elements_array = rModelPart.Elements();
        const auto &r_conditions_array = rModelPart.Conditions();
        const auto &r_constraints_array = rModelPart.MasterSlaveConstraints();
        const int n_elems = static_cast<int>(r_elements_array.size());
        const int n_conds = static_cast<int>(r_conditions_array.size());
        const int n_constraints = static_cast<int>(r_constraints_array.size());

        // Global dof set
        DofSetType dof_global_set;
        dof_global_set.reserve(n_elems*20);

        // Auxiliary DOFs list
        DofsVectorType dof_list;
        DofsVectorType second_dof_list; // The second_dof_list is only used on constraints to include master/slave relations

#pragma omp parallel firstprivate(dof_list, second_dof_list)
        {
            const auto& r_process_info = rModelPart.GetProcessInfo();

            // We cleate the temporal set and we reserve some space on them
            DofSetType dofs_tmp_set;
            dofs_tmp_set.reserve(20000);

            // Get the DOFs list from each element and insert it in the temporary set
#pragma omp for schedule(guided, 512) nowait
            for (int i_elem = 0; i_elem < n_elems; ++i_elem) {
                const auto it_elem = r_elements_array.begin() + i_elem;
                it_elem->GetDofList(dof_list, r_process_info);
                dofs_tmp_set.insert(dof_list.begin(), dof_list.end());
            }

            // Get the DOFs list from each condition and insert it in the temporary set
#pragma omp for schedule(guided, 512) nowait
            for (int i_cond = 0; i_cond < n_conds; ++i_cond) {
                const auto it_cond = r_conditions_array.begin() + i_cond;
                it_cond->GetDofList(dof_list, r_process_info);
                dofs_tmp_set.insert(dof_list.begin(), dof_list.end());
            }

            // Get the DOFs list from each constraint and insert it in the temporary set
#pragma omp for  schedule(guided, 512) nowait
            for (int i_const = 0; i_const < n_constraints; ++i_const) {
                auto it_const = r_constraints_array.begin() + i_const;
                it_const->GetDofList(dof_list, second_dof_list, r_process_info);
                dofs_tmp_set.insert(dof_list.begin(), dof_list.end());
                dofs_tmp_set.insert(second_dof_list.begin(), second_dof_list.end());
            }

            // We merge all the sets in one thread
#pragma omp critical
            {
                dof_global_set.insert(dofs_tmp_set.begin(), dofs_tmp_set.end());
            }
        }

        KRATOS_INFO_IF("ExplicitBuilder", ( this->GetEchoLevel() > 2)) << "Initializing ordered array filling\n" << std::endl;

        // Ordering the global DOF set
        mDofSet = DofsArrayType();
        DofsArrayType temp_dof_set;
        temp_dof_set.reserve(dof_global_set.size());
        for (auto it_dof = dof_global_set.begin(); it_dof != dof_global_set.end(); ++it_dof) {
            temp_dof_set.push_back(*it_dof);
        }
        temp_dof_set.Sort();
        mDofSet = temp_dof_set;
        mEquationSystemSize = mDofSet.size();

        // DoFs set checks
        // Throws an exception if there are no Degrees Of Freedom involved in the analysis
        KRATOS_ERROR_IF(mDofSet.size() == 0) << "No degrees of freedom!" << std::endl;

        // Check if each DOF has a reaction. Note that the explicit residual is stored in these
        for (auto it_dof = mDofSet.begin(); it_dof != mDofSet.end(); ++it_dof) {
            KRATOS_ERROR_IF_NOT(it_dof->HasReaction()) << "Reaction variable not set for the following : " << std::endl
                << "Node : " << it_dof->Id() << std::endl
                << "Dof : " << (*it_dof) << std::endl << "Not possible to calculate reactions." << std::endl;
        }

        mDofSetIsInitialized = true;

        KRATOS_INFO_IF("ExplicitBuilder", ( this->GetEchoLevel() > 2)) << "Number of degrees of freedom:" << mDofSet.size() << std::endl;
        KRATOS_INFO_IF("ExplicitBuilder", ( this->GetEchoLevel() > 2 && rModelPart.GetCommunicator().MyPID() == 0)) << "Finished setting up the dofs" << std::endl;
        KRATOS_INFO_IF("ExplicitBuilder", ( this->GetEchoLevel() > 2)) << "End of setup dof set\n" << std::endl;

        KRATOS_CATCH("");
    }

    /**
     * @brief Set the Up Dof Set Equation Ids object
     * Set up the DOF set equation ids
     */
    virtual void SetUpDofSetEquationIds()
    {
        // Firstly check that the DOF set is initialized
        KRATOS_ERROR_IF_NOT(mDofSetIsInitialized) << "Trying to set the equation ids. before initializing the DOF set. Please call the SetUpDofSet() before." << std::endl;
        KRATOS_ERROR_IF(mEquationSystemSize == 0) << "Trying to set the equation ids. in an empty DOF set (equation system size is 0)." << std::endl;

        // Loop the DOF set to assign the equation ids
        IndexPartition<unsigned int>(mEquationSystemSize).for_each(
            [&](unsigned int i_dof){
                auto it_dof = mDofSet.begin() + i_dof;
                it_dof->SetEquationId(i_dof);
            }
        );
    }

    /**
     * @brief Set the Up Lumped Mass Vector object
     * This method sets up the lumped mass vector used in the explicit update.
     * Note that it requires that the equation ids. are already set and the
     * implementation of the mass contributions to be done in the element level
     * in the CalculateLumpedMassVector method.
     * @param rModelPart The model part to compute
     */
    virtual void SetUpLumpedMassVector(const ModelPart &rModelPart)
    {
        KRATOS_TRY;

        KRATOS_INFO_IF("ExplicitBuilder", this->GetEchoLevel() > 1) << "Setting up the lumped mass matrix vector" << std::endl;

        // Initialize the lumped mass matrix vector
        // Note that the lumped mass matrix vector size matches the dof set one
        mpLumpedMassVector = TSystemVectorPointerType(new TSystemVectorType(GetDofSet().size()));
        TDenseSpace::SetToZero(*mpLumpedMassVector);

        // Loop the elements to get the lumped mass vector
        const auto &r_elements_array = rModelPart.Elements();
        const auto &r_process_info = rModelPart.GetProcessInfo();
        const std::size_t n_elems = r_elements_array.size();

        // Auxiliary definitions
        struct TLS {
            LocalSystemVectorType elem_mass_vector;
            std::vector<std::size_t> elem_equation_id;
        };

        // Iterate over elements
        IndexPartition<std::size_t>(n_elems).for_each(TLS(), [&](std::size_t i_elem, TLS& rTLS){
            const auto it_elem = r_elements_array.begin() + i_elem;

            // Calculate the elemental lumped mass vector
            it_elem->CalculateLumpedMassVector(rTLS.elem_mass_vector, r_process_info);
            it_elem->EquationIdVector(rTLS.elem_equation_id, r_process_info);

            // Update value of lumped mass vector
            for (IndexType i = 0; i < rTLS.elem_equation_id.size(); ++i) {
                AtomicAdd((*mpLumpedMassVector)[rTLS.elem_equation_id[i]], rTLS.elem_mass_vector(i));
            }
        });

        // Set the lumped mass vector flag as true
        mLumpedMassVectorIsInitialized = true;

        KRATOS_CATCH("");
    }

    /**
     * @brief Initialize the DOF set reactions
     * For an already initialized dof set (mDofSet), this method sets to
     * zero the corresponding reaction variable values. Note that in the
     * explicit build the reactions are used as residual container.
     */
    virtual void InitializeDofSetReactions()
    {
        // Firstly check that the DOF set is initialized
        KRATOS_ERROR_IF_NOT(mDofSetIsInitialized) << "Trying to initialize the explicit residual but the DOFs set is not initialized yet." << std::endl;
        KRATOS_ERROR_IF(mEquationSystemSize == 0) << "Trying to set the equation ids. in an empty DOF set (equation system size is 0)." << std::endl;

        // Loop the reactions to initialize them to zero
        block_for_each(
            mDofSet,
            [](DofType& rDof){
                rDof.GetSolutionStepReactionValue() = 0.0;
            }
        );
    }

    /**
     * @brief It computes the reactions of the system
     * @param rModelPart The model part to compute
     */
    virtual void CalculateReactions()
    {
        if (mCalculateReactionsFlag) {
            // Firstly check that the DOF set is initialized
            KRATOS_ERROR_IF_NOT(mDofSetIsInitialized) << "Trying to initialize the explicit residual but the DOFs set is not initialized yet." << std::endl;
            KRATOS_ERROR_IF(mEquationSystemSize == 0) << "Trying to set the equation ids. in an empty DOF set (equation system size is 0)." << std::endl;

            // Calculate the reactions as minus the current value
            // Note that we take advantage of the fact that Kratos always works with a residual based formulation
            // This means that the reactions are minus the residual. As we use the reaction as residual container
            // during the explicit resolution of the problem, the calculate reactions is as easy as switching the sign
            block_for_each(
                mDofSet,
                [](DofType& rDof){
                    auto& r_reaction_value = rDof.GetSolutionStepReactionValue();
                    r_reaction_value *= -1.0;
                }
            );
        }
    }

    /**
     * @brief This method validate and assign default parameters
     * @param rParameters Parameters to be validated
     * @param DefaultParameters The default parameters
     * @return Returns validated Parameters
     */
    virtual Parameters ValidateAndAssignParameters(
        Parameters ThisParameters,
        const Parameters DefaultParameters
        ) const
    {
        ThisParameters.ValidateAndAssignDefaults(DefaultParameters);
        return ThisParameters;
    }

    /**
     * @brief This method assigns settings to member variables
     * @param ThisParameters Parameters that are assigned to the member variables
     */
    virtual void AssignSettings(const Parameters ThisParameters)
    {
    }

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

    static std::vector<Internals::RegisteredPrototypeBase<ClassType>> msPrototypes;

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
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}
}; /* Class ExplicitBuilder */
///@}
///@name Type Definitions
///@{

///@}
} /* namespace Kratos.*/
