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

#pragma once

// System includes

// External includes

// Project includes
#include "containers/csr_matrix.h"
#include "containers/system_vector.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "utilities/builtin_timer.h"
#include "utilities/dof_utilities/assembly_helper.h" //TODO: we should move this to somewhere else
#include "utilities/dof_utilities/block_build_dof_array_utility.h"
#include "utilities/dof_utilities/elimination_build_dof_array_utility.h"
#include "utilities/entities_utilities.h"
#include "utilities/openmp_utils.h" //TODO: SOME FILES INCLUDING scheme.h RELY ON THIS. LEAVING AS FUTURE TODO.
#include "utilities/parallel_utilities.h"
#include "utilities/timer.h"

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
 * @class NewScheme
 * @ingroup KratosCore
 * @brief This class provides the implementation of the basic tasks that are needed by the solution strategy.
 * @details It is intended to be the place for tailoring the solution strategies to problem specific tasks.
 * @author Ruben Zorrilla
 */
//TODO: Think about the template parameters
template<class TDataType = double, class TIndexType=std::size_t, class TMatrixType=CsrMatrix<TDataType, TIndexType>, class TVectorType=SystemVector<TDataType, TIndexType>>
class NewScheme
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of NewScheme
    KRATOS_CLASS_POINTER_DEFINITION(NewScheme);

    /// The definition of the current class
    using ClassType = NewScheme;

    /// Index type definition
    using IndexType = TIndexType;

    /// Size type definition
    using SizeType = std::size_t;

    /// Data type definition
    using DataType = TDataType;

    /// Matrix type definition
    using SystemMatrixType = TMatrixType;

    /// Matrix pointer type definition
    using SystemMatrixPointerType = typename SystemMatrixType::Pointer;

    /// Vector type definition
    using SystemVectorType = TVectorType;

    /// Vector type definition
    using SystemVectorPointerType = typename SystemVectorType::Pointer;

    /// Dense space definition
    using DenseSpaceType = UblasSpace<double, Matrix, Vector>;

    /// Local system matrix type definition
    using LocalSystemMatrixType = DenseSpaceType::MatrixType;

    /// Local system vector type definition
    using LocalSystemVectorType = DenseSpaceType::VectorType;

    /// DoF type definition
    using DofType = Dof<double>;

    /// DoF array type definition
    using DofsArrayType = ModelPart::DofsArrayType;

    /// Elements containers definition
    using ElementsArrayType = ModelPart::ElementsContainerType;

    /// Conditions containers definition
    using ConditionsArrayType = ModelPart::ConditionsContainerType;

    /// Build type definition
    enum class BuildType
    {
        Block,
        Elimination,
        Custom
    };

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default Constructor
     * @details Initializes the flags
     */
    explicit NewScheme()
    {
        mSchemeIsInitialized = false;
        mElementsAreInitialized = false;
        mConditionsAreInitialized = false;
    }

    /**
     * @brief Constructor with Parameters
     */
    explicit NewScheme(Parameters ThisParameters)
    {
        // Validate default parameters
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);

        mSchemeIsInitialized = false;
        mElementsAreInitialized = false;
        mConditionsAreInitialized = false;
    }

    /** Copy Constructor.
     */
    explicit NewScheme(NewScheme& rOther)
      :mSchemeIsInitialized(rOther.mSchemeIsInitialized)
      ,mElementsAreInitialized(rOther.mElementsAreInitialized)
      ,mConditionsAreInitialized(rOther.mConditionsAreInitialized)
    {
    }

    /** Destructor.
     */
    virtual ~NewScheme()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Create method
     * @param ThisParameters The configuration parameters
     */
    virtual typename ClassType::Pointer Create(Parameters ThisParameters) const
    {
        return Kratos::make_shared<ClassType>(ThisParameters);
    }

    /**
     * @brief Clone method
     * @return The pointer of the cloned NewScheme
     */
    virtual Pointer Clone()
    {
        return Kratos::make_shared<NewScheme>(*this) ;
    }

    /**
     * @brief This is the place to initialize the NewScheme.
     * @details This is intended to be called just once when the strategy is initialized
     * @param rModelPart The model part of the problem to solve
     */
    virtual void Initialize(ModelPart& rModelPart)
    {
        KRATOS_TRY

        mSchemeIsInitialized = true;

        KRATOS_CATCH("")
    }

    /**
     * @brief This method returns if the scheme is initialized
     * @return True if initialized, false otherwise
     */
    bool SchemeIsInitialized()
    {
        return mSchemeIsInitialized;
    }

    /**
     * @brief This method sets if the elements have been initialized or not (true by default)
     * @param ElementsAreInitializedFlag If the flag must be set to true or false
     */
    void SetSchemeIsInitialized(bool SchemeIsInitializedFlag = true)
    {
        mSchemeIsInitialized = SchemeIsInitializedFlag;
    }

    /**
     * @brief This method returns if the elements are initialized
     * @return True if initialized, false otherwise
     */
    bool ElementsAreInitialized()
    {
        return mElementsAreInitialized;
    }

    /**
     * @brief This method sets if the elements have been initialized or not (true by default)
     * @param ElementsAreInitializedFlag If the flag must be set to true or false
     */
    void SetElementsAreInitialized(bool ElementsAreInitializedFlag = true)
    {
        mElementsAreInitialized = ElementsAreInitializedFlag;
    }

    /**
     * @brief This method returns if the conditions are initialized
     * @return True if initialized, false otherwise
     */
    bool ConditionsAreInitialized()
    {
        return mConditionsAreInitialized;
    }

    /**
     * @brief This method sets if the conditions have been initialized or not (true by default)
     * @param ConditionsAreInitializedFlag If the flag must be set to true or false
     */
    void SetConditionsAreInitialized(bool ConditionsAreInitializedFlag = true)
    {
        mConditionsAreInitialized = ConditionsAreInitializedFlag;
    }

    /**
     * @brief This is the place to initialize the elements.
     * @details This is intended to be called just once when the strategy is initialized
     * @param rModelPart The model part of the problem to solve
     */
    virtual void InitializeElements( ModelPart& rModelPart)
    {
        KRATOS_TRY

        EntitiesUtilities::InitializeEntities<Element>(rModelPart);

        SetElementsAreInitialized();

        KRATOS_CATCH("")
    }

    /**
     * @brief This is the place to initialize the conditions.
     * @details This is intended to be called just once when the strategy is initialized
     * @param rModelPart The model part of the problem to solve
     */
    virtual void InitializeConditions(ModelPart& rModelPart)
    {
        KRATOS_TRY

        KRATOS_ERROR_IF_NOT(mElementsAreInitialized) << "Before initializing Conditions, initialize Elements FIRST" << std::endl;

        EntitiesUtilities::InitializeEntities<Condition>(rModelPart);

        SetConditionsAreInitialized();

        KRATOS_CATCH("")
    }

    /**
     * @brief Function called once at the beginning of each solution step.
     * @details The basic operations to be carried in there are the following:
     * - managing variables to be kept constant over the time step (for example time-Scheme constants depending on the actual time step)
     * @param rModelPart The model part of the problem to solve
     * @param A LHS matrix
     * @param Dx Incremental update of primary variables
     * @param b RHS Vector
     */
    virtual void InitializeSolutionStep(
        ModelPart& rModelPart,
        SystemMatrixType& A,
        SystemVectorType& Dx,
        SystemVectorType& b)
    {
        KRATOS_TRY

        // Initializes solution step for all of the elements, conditions and constraints
        EntitiesUtilities::InitializeSolutionStepAllEntities(rModelPart);

        KRATOS_CATCH("")
    }

    /**
     * @brief Function called once at the end of a solution step, after convergence is reached if an iterative process is needed
     * @param rModelPart The model part of the problem to solve
     * @param A LHS matrix
     * @param Dx Incremental update of primary variables
     * @param b RHS Vector
     */
    virtual void FinalizeSolutionStep(
        ModelPart& rModelPart,
        SystemMatrixType& A,
        SystemVectorType& Dx,
        SystemVectorType& b)
    {
        KRATOS_TRY

        // Finalizes solution step for all of the elements, conditions and constraints
        EntitiesUtilities::FinalizeSolutionStepAllEntities(rModelPart);

        KRATOS_CATCH("")
    }

    /**
     * @brief unction to be called when it is needed to initialize an iteration. It is designed to be called at the beginning of each non linear iteration
     * @note Take care: the elemental function with the same name is NOT called here.
     * @warning Must be defined in derived classes
     * @details The function is called in the builder for memory efficiency
     * @param rModelPart The model part of the problem to solve
     * @param A LHS matrix
     * @param Dx Incremental update of primary variables
     * @param b RHS Vector
     */
    virtual void InitializeNonLinIteration(
        ModelPart& rModelPart,
        SystemMatrixType& A,
        SystemVectorType& Dx,
        SystemVectorType& b)
    {
        KRATOS_TRY

        // Finalizes non-linear iteration for all of the elements, conditions and constraints
        EntitiesUtilities::InitializeNonLinearIterationAllEntities(rModelPart);

        KRATOS_CATCH("")
    }

    /**
     * @brief Function to be called when it is needed to finalize an iteration. It is designed to be called at the end of each non linear iteration
     * @param rModelPart The model part of the problem to solve
     * @param A LHS matrix
     * @param Dx Incremental update of primary variables
     * @param b RHS Vector
     */
    virtual void FinalizeNonLinIteration(
        ModelPart& rModelPart,
        SystemMatrixType& A,
        SystemVectorType& Dx,
        SystemVectorType& b)
    {
        KRATOS_TRY

        // Finalizes non-linear iteration for all of the elements, conditions and constraints
        EntitiesUtilities::FinalizeNonLinearIterationAllEntities(rModelPart);

        KRATOS_CATCH("")
    }

    virtual void SetUpDofArray(const ModelPart& rModelPart)
    {
        //TODO: Discuss on how to provide the custom DOF set function
        // Call external utility to perform the build
        if (mBuildType == BuildType::Block) {
            BlockBuildDofArrayUtility::SetUpDofArray(rModelPart, mDofSet, mEchoLevel, mCalculateReactionsFlag);
        } else if (mBuildType == BuildType::Elimination) {
            EliminationBuildDofArrayUtility::SetUpDofArray(rModelPart, mDofSet, mEchoLevel, mCalculateReactionsFlag);
        } else {
            KRATOS_ERROR << "Build type not supported." << std::endl;
        }

        // Save the size of the system based on the number of DOFs
        mEquationSystemSize = mDofSet.size();

        // Set the flag as already initialized
        mDofSetIsInitialized = true;

        KRATOS_INFO_IF("NewScheme", mEchoLevel >= 2) << "Finished DOFs array set up." << std::endl;
    }

    //TODO: think about renaming this to SetUpEquationIds (maybe we can just keep it because we are already used to it)
    void SetUpSystem(const ModelPart& rModelPart)
    {
        // Set up the DOFs equation global ids
        IndexPartition<IndexType>(mEquationSystemSize).for_each([&, this](IndexType Index) {
            auto dof_iterator = this->mDofSet.begin() + Index;
            dof_iterator->SetEquationId(Index);
        });

        KRATOS_INFO_IF("NewScheme", mEchoLevel >= 2) << "Finished system set up." << std::endl;
    }

    //TODO: I think this belongs to the strategy
    // // TODO: think about renaming this to ResizeAndInitializeSystem (maybe we can just keep it because we are already used to it)
    // void ResizeAndInitializeVectors(
    //     const ModelPart& rModelPart,
    //     SystemMatrixPointerType& rpLHS,
    //     SystemVectorPointerType& rpDx,
    //     SystemVectorPointerType& rpRHS)
    // {
    //     KRATOS_TRY

    //     if (rpLHS == nullptr) {
    //         auto p_aux = Kratos::make_shared<SystemMatrixPointerType>();
    //     }

    //     // if (pA == NULL) // if the pointer is not initialized initialize it to an empty matrix
    //     // {
    //     //     TSystemMatrixPointerType pNewA = TSystemMatrixPointerType(new TSystemMatrixType(0, 0));
    //     //     pA.swap(pNewA);
    //     // }
    //     // if (pDx == NULL) // if the pointer is not initialized initialize it to an empty matrix
    //     // {
    //     //     TSystemVectorPointerType pNewDx = TSystemVectorPointerType(new TSystemVectorType(0));
    //     //     pDx.swap(pNewDx);
    //     // }
    //     // if (pb == NULL) // if the pointer is not initialized initialize it to an empty matrix
    //     // {
    //     //     TSystemVectorPointerType pNewb = TSystemVectorPointerType(new TSystemVectorType(0));
    //     //     pb.swap(pNewb);
    //     // }

    //     // TSystemMatrixType &A = *pA;
    //     // TSystemVectorType &Dx = *pDx;
    //     // TSystemVectorType &b = *pb;

    //     // // resizing the system vectors and matrix
    //     // if (A.size1() == 0 || BaseType::GetReshapeMatrixFlag() == true) // if the matrix is not initialized
    //     // {
    //     //     A.resize(BaseType::mEquationSystemSize, BaseType::mEquationSystemSize, false);
    //     //     ConstructMatrixStructure(pScheme, A, rModelPart);
    //     // }
    //     // else
    //     // {
    //     //     if (A.size1() != BaseType::mEquationSystemSize || A.size2() != BaseType::mEquationSystemSize)
    //     //     {
    //     //         KRATOS_ERROR << "The equation system size has changed during the simulation. This is not permitted." << std::endl;
    //     //         A.resize(BaseType::mEquationSystemSize, BaseType::mEquationSystemSize, true);
    //     //         ConstructMatrixStructure(pScheme, A, rModelPart);
    //     //     }
    //     // }
    //     // if (Dx.size() != BaseType::mEquationSystemSize)
    //     //     Dx.resize(BaseType::mEquationSystemSize, false);
    //     // TSparseSpace::SetToZero(Dx);
    //     // if (b.size() != BaseType::mEquationSystemSize)
    //     // {
    //     //     b.resize(BaseType::mEquationSystemSize, false);
    //     // }
    //     // TSparseSpace::SetToZero(b);

    //     // ConstructMasterSlaveConstraintsStructure(rModelPart);

    //     KRATOS_INFO_IF("NewScheme", mEchoLevel >= 2) << "Finished system initialization." << std::endl;

    //     KRATOS_CATCH("")
    // }

    virtual void Build(
        const ModelPart& rModelPart,
        SystemMatrixType& rLHS,
        SystemVectorType& rRHS)
    {
        Timer::Start("Build");

        const auto timer = BuiltinTimer();

        struct TLS
        {
            LocalSystemMatrixType LocLhs; // Local LHS contribution
            LocalSystemVectorType LocRhs; // Local RHS constribution
            Element::EquationIdVectorType LocEqIds; // Vector containing the localization in the system of the different terms
        };

        const auto elem_func = [](ModelPart::ElementConstantIterator ItElem, const ProcessInfo& rProcessInfo, TLS& rTLS){
            // Calculate local LHS and RHS contributions
            ItElem->CalculateLocalSystem(rTLS.LocLhs, rTLS.LocRhs, rProcessInfo);
        };

        const auto cond_func = [](ModelPart::ConditionConstantIterator ItCond, const ProcessInfo& rProcessInfo, TLS& rTLS){
            // Calculate local LHS and RHS contributions
            ItCond->CalculateLocalSystem(rTLS.LocLhs, rTLS.LocRhs, rProcessInfo);
        };

        TLS aux_tls;
        AssemblyHelper<TLS> assembly_helper(rModelPart);
        assembly_helper.SetElementAssemblyFunction(elem_func);
        assembly_helper.SetConditionAssemblyFunction(cond_func);
        assembly_helper.AssembleLocalSystem(rLHS, rRHS, aux_tls);
        // assembly_helper.AssembleLocalSystemElements<>(rModelPart, elem_func, rLHS, rRHS, aux_tls);
        // assembly_helper.AssembleLocalSystemConditions<>(rModelPart, cond_func, rLHS, rRHS, aux_tls);

        KRATOS_INFO_IF("NewScheme", mEchoLevel >= 1) << "Build time: " << timer << std::endl;
        KRATOS_INFO_IF("NewScheme", mEchoLevel >= 2) << "Finished parallel building" << std::endl;

        Timer::Stop("Build");


        // AssemblyHelper<TLS1,TLS2>::Assemble(ElmeConainter,
        //                                     [](ModelPart::ElementConstantIterator ItElem, const ProcessInfo& rProcessInfo, TLS& rTLS){
        //                                         // Calculate local LHS and RHS contributions
        //                                         ItElem->CalculateLocalSystem(rTLS.LocLhs, rTLS.LocRhs, rProcessInfo),
        //                                     [](ModelPart::ElementConstantIterator ItElem, const ProcessInfo& rProcessInfo, TLS& rTLS){
        //                                         // Calculate local LHS and RHS contributions
        //                                         ItElem->CalculateLocalSystem(rTLS.LocLhs, rTLS.LocRhs, rProcessInfo),
        //                                     CondFunc)

        // AssemblyHelper<TLS1,TLS2>::Assemble(ElmeConainter, ElemFunc, CondContainer, CondFunc)

    }

    /**
     * @brief Performing the prediction of the solution.
     * @warning Must be defined in derived classes
     * @param rModelPart The model part of the problem to solve
     * @param A LHS matrix
     * @param Dx Incremental update of primary variables
     * @param b RHS Vector
     */
    virtual void Predict(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        SystemMatrixType& A,
        SystemVectorType& Dx,
        SystemVectorType& b)
    {
        KRATOS_TRY
        KRATOS_CATCH("")
    }

    /**
     * @brief Performing the update of the solution.
     * @warning Must be defined in derived classes
     * @param rModelPart The model part of the problem to solve
     * @param rDofSet Set of all primary variables
     * @param A LHS matrix
     * @param Dx Incremental update of primary variables
     * @param b RHS Vector
     */
    virtual void Update(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        SystemMatrixType& A,
        SystemVectorType& Dx,
        SystemVectorType& b)
    {
        KRATOS_TRY
        KRATOS_CATCH("")
    }

    /**
     * @brief Functions to be called to prepare the data needed for the output of results.
     * @warning Must be defined in derived classes
     * @param rModelPart The model part of the problem to solve
     * @param rDofSet Set of all primary variables
     * @param A LHS matrix
     * @param Dx Incremental update of primary variables
     * @param b RHS Vector
     */
    virtual void CalculateOutputData(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        SystemMatrixType& A,
        SystemVectorType& Dx,
        SystemVectorType& b)
    {
        KRATOS_TRY
        KRATOS_CATCH("")
    }

    /**
     * @brief Functions that cleans the results data.
     * @warning Must be implemented in the derived classes
     */
    virtual void CleanOutputData()
    {
        KRATOS_TRY
        KRATOS_CATCH("")
    }

    /**
     * @brief This function is intended to be called at the end of the solution step to clean up memory storage not needed after the end of the solution step
     * @warning Must be implemented in the derived classes
     */
    virtual void Clean()
    {
        KRATOS_TRY
        KRATOS_CATCH("")
    }

    /**
     * @brief Liberate internal storage.
     * @warning Must be implemented in the derived classes
     */
    virtual void Clear()
    {
        KRATOS_TRY
        KRATOS_CATCH("")
    }

    /**
     * @brief This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @details Checks can be "expensive" as the function is designed
     * @param rModelPart The model part of the problem to solve
     * @return 0 all OK, 1 otherwise
     */
    virtual int Check(const ModelPart& rModelPart) const
    {
        KRATOS_TRY
        return 0;
        KRATOS_CATCH("");
    }

    virtual int Check(ModelPart& rModelPart)
    {
        // calling the const version for backward compatibility
        const NewScheme& r_const_this = *this;
        const ModelPart& r_const_model_part = rModelPart;
        return r_const_this.Check(r_const_model_part);
    }

    /**
     * @brief This function is designed to be called in the builder and solver to introduce the selected time integration scheme.
     * @details It "asks" the matrix needed to the element and performs the operations needed to introduce the selected time integration scheme. This function calculates at the same time the contribution to the LHS and to the RHS of the system
     * @param rElement The element to compute
     * @param LHS_Contribution The LHS matrix contribution
     * @param RHS_Contribution The RHS vector contribution
     * @param rEquationIdVector The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    virtual void CalculateSystemContributions(
        Element& rElement,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& rEquationIdVector,
        const ProcessInfo& rCurrentProcessInfo)
    {
        rElement.CalculateLocalSystem(LHS_Contribution, RHS_Contribution, rCurrentProcessInfo);
    }

    /**
     * @brief Functions totally analogous to the precedent but applied to the "condition" objects
     * @param rCondition The condition to compute
     * @param LHS_Contribution The LHS matrix contribution
     * @param RHS_Contribution The RHS vector contribution
     * @param rEquationIdVector The ID's of the condition degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    virtual void CalculateSystemContributions(
        Condition& rCondition,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& rEquationIdVector,
        const ProcessInfo& rCurrentProcessInfo)
    {
        rCondition.CalculateLocalSystem(LHS_Contribution, RHS_Contribution, rCurrentProcessInfo);
    }

    /**
     * @brief This function is designed to calculate just the RHS contribution
     * @param rElement The element to compute
     * @param RHS_Contribution The RHS vector contribution
     * @param rEquationIdVector The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    virtual void CalculateRHSContribution(
        Element& rElement,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& rEquationIdVector,
        const ProcessInfo& rCurrentProcessInfo)
    {
        rElement.CalculateRightHandSide(RHS_Contribution, rCurrentProcessInfo);
    }

    /**
     * @brief Functions totally analogous to the precedent but applied to the "condition" objects
     * @param rCondition The condition to compute
     * @param RHS_Contribution The RHS vector contribution
     * @param rEquationIdVector The ID's of the condition degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    virtual void CalculateRHSContribution(
        Condition& rCondition,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& rEquationIdVector,
        const ProcessInfo& rCurrentProcessInfo)
    {
        rCondition.CalculateRightHandSide(RHS_Contribution, rCurrentProcessInfo);
    }

    /**
     * @brief This function is designed to calculate just the LHS contribution
     * @param rElement The element to compute
     * @param LHS_Contribution The RHS vector contribution
     * @param rEquationIdVector The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    virtual void CalculateLHSContribution(
        Element& rElement,
        LocalSystemMatrixType& LHS_Contribution,
        Element::EquationIdVectorType& rEquationIdVector,
        const ProcessInfo& rCurrentProcessInfo)
    {
        rElement.CalculateLeftHandSide(LHS_Contribution, rCurrentProcessInfo);
    }

    /**
     * @brief Functions totally analogous to the precedent but applied to the "condition" objects
     * @param rCondition The condition to compute
     * @param LHS_Contribution The RHS vector contribution
     * @param rEquationIdVector The ID's of the condition degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    virtual void CalculateLHSContribution(
        Condition& rCondition,
        LocalSystemMatrixType& LHS_Contribution,
        Element::EquationIdVectorType& rEquationIdVector,
        const ProcessInfo& rCurrentProcessInfo)
    {
        rCondition.CalculateLeftHandSide(LHS_Contribution, rCurrentProcessInfo);
    }

    /**
     * @brief This method gets the eqaution id corresponding to the current element
     * @param rElement The element to compute
     * @param rEquationId The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    virtual void EquationId(
        const Element& rElement,
        Element::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo)
    {
        rElement.EquationIdVector(rEquationId, rCurrentProcessInfo);
    }

    /**
     * @brief Functions totally analogous to the precedent but applied to the "condition" objects
     * @param rCondition The condition to compute
     * @param rEquationId The ID's of the condition degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    virtual void EquationId(
        const Condition& rCondition,
        Element::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo)
    {
        rCondition.EquationIdVector(rEquationId, rCurrentProcessInfo);
    }

    /**
     * @brief Function that returns the list of Degrees of freedom to be assembled in the system for a Given element
     * @param pCurrentElement The element to compute
     * @param rDofList The list containing the element degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    virtual void GetDofList(
        const Element& rElement,
        Element::DofsVectorType& rDofList,
        const ProcessInfo& rCurrentProcessInfo)
    {
        rElement.GetDofList(rDofList, rCurrentProcessInfo);
    }

    /**
     * @brief Function that returns the list of Degrees of freedom to be assembled in the system for a Given condition
     * @param rCondition The condition to compute
     * @param rDofList The list containing the condition degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    virtual void GetDofList(
        const Condition& rCondition,
        Element::DofsVectorType& rDofList,
        const ProcessInfo& rCurrentProcessInfo)
    {
        rCondition.GetDofList(rDofList, rCurrentProcessInfo);
    }

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    virtual Parameters GetDefaultParameters() const
    {
        const Parameters default_parameters = Parameters(R"({
            "name" : "new_scheme",
            "build_type" : "block",
            "echo_level" : 0
        })");
        return default_parameters;
    }

    /**
     * @brief Returns the name of the class as used in the settings (snake_case format)
     * @return The name of the class
     */
    static std::string Name()
    {
        return "new_scheme";
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
        return "Scheme";
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

    bool mSchemeIsInitialized;      /// Flag to be used in controlling if the Scheme has been initialized or not
    bool mElementsAreInitialized;   /// Flag taking in account if the elements were initialized correctly or not
    bool mConditionsAreInitialized; /// Flag taking in account if the conditions were initialized correctly or not

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief This method validate and assign default parameters
     * @param rParameters Parameters to be validated
     * @param DefaultParameters The default parameters
     * @return Returns validated Parameters
     */
    virtual Parameters ValidateAndAssignParameters(
        Parameters ThisParameters,
        const Parameters DefaultParameters) const
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
        // Set build type
        const std::string build_type = ThisParameters["build_type"].GetString();
        if (build_type == "block") {
            mBuildType = BuildType::Block;
        } else if (build_type == "elimination") {
            mBuildType = BuildType::Elimination;
        } else if (build_type == "custom") {
            mBuildType = BuildType::Custom;
        } else {
            KRATOS_ERROR << "Provided 'build_type' is '" << build_type << "'. Available options are:\n"
            << "\t- 'block'\n"
            << "\t- 'elimination'\n"
            << "\t- 'custom'\n" << std::endl;
        }

        // Set verbosity level
        mEchoLevel = ThisParameters["echo_level"].GetInt();
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

    ///@}
    ///@name Member Variables
    ///@{

    BuildType mBuildType;

    DofsArrayType mDofSet;

    SizeType mEquationSystemSize; /// Number of degrees of freedom of the problem to be solve

    bool mReshapeMatrixFlag = false; /// If the matrix is reshaped each step

    bool mDofSetIsInitialized = false; /// Flag taking care if the dof set was initialized ot not

    bool mCalculateReactionsFlag = false; /// Flag taking in account if it is needed or not to calculate the reactions

    int mEchoLevel = 0;

    // TSystemVectorPointerType mpReactionsVector; //FIXME:

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

}; // Class Scheme

} // namespace Kratos.
