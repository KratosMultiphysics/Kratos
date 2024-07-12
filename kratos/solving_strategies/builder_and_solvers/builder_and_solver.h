//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//
//

#pragma once

/* System includes */
#include <set>

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "includes/kratos_parameters.h"

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
 * @class BuilderAndSolver
 * @ingroup KratosCore
 * @brief Current class provides an implementation for the base builder and solving operations.
 * @details The RHS is constituted by the unbalanced loads (residual)
 * Degrees of freedom are reordered putting the restrained degrees of freedom at
 * the end of the system ordered in reverse order with respect to the DofSet.
 * @author Riccardo Rossi
 */
template<class TSparseSpace,
         class TDenseSpace, // = DenseSpace<double>,
         class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         >
class BuilderAndSolver
{
public:
    ///@name Type Definitions
    ///@{

    /// The definition of the current class
    typedef BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> ClassType;

    /// Definition of the size type
    typedef std::size_t SizeType;

    /// Definition of the index type
    typedef std::size_t IndexType;

    /// Definition of the data type
    typedef typename TSparseSpace::DataType TDataType;

    ///Definition of the sparse matrix
    typedef typename TSparseSpace::MatrixType TSystemMatrixType;

    /// Definition of the vector size
    typedef typename TSparseSpace::VectorType TSystemVectorType;

    /// Definition of the pointer to the sparse matrix
    typedef typename TSparseSpace::MatrixPointerType TSystemMatrixPointerType;

    /// Definition of the pointer to the vector
    typedef typename TSparseSpace::VectorPointerType TSystemVectorPointerType;

    /// The local matrix definition
    typedef typename TDenseSpace::MatrixType LocalSystemMatrixType;

    /// The local vector definition
    typedef typename TDenseSpace::VectorType LocalSystemVectorType;

    /// Definition of the scheme type
    typedef Scheme<TSparseSpace, TDenseSpace> TSchemeType;

    /// Definition of the DoF class
    typedef ModelPart::DofType TDofType;

    /// Definition of the DoF array type
    typedef ModelPart::DofsArrayType DofsArrayType;

    /// The containers of the entities
    typedef ModelPart::NodesContainerType NodesArrayType;
    typedef ModelPart::ElementsContainerType ElementsArrayType;
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;

    /// The definition of the element container type
    typedef PointerVectorSet<Element, IndexedObject> ElementsContainerType;

    /// Pointer definition of BuilderAndSolver
    KRATOS_CLASS_POINTER_DEFINITION(BuilderAndSolver);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     */
    explicit BuilderAndSolver()
    {
    }

    /**
     * @brief Default constructor with Parameters
     * @param pNewLinearSystemSolver The linear solver for the system of equations
     * @param ThisParameters The configuration parameters
     */
    explicit BuilderAndSolver(
        typename TLinearSolver::Pointer pNewLinearSystemSolver,
        Parameters ThisParameters
        )
    {
        // Validate and assign defaults
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);

        // We set the other member variables
        mpLinearSystemSolver = pNewLinearSystemSolver;
    }

    /**
     * @brief Default constructor.
     * @param pNewLinearSystemSolver The linear solver for the system of equations
     */
    explicit BuilderAndSolver(typename TLinearSolver::Pointer pNewLinearSystemSolver)
    {
        mpLinearSystemSolver = pNewLinearSystemSolver;
    }

    /** Destructor.
     */
    virtual ~BuilderAndSolver()
    {
    }

    /**
     * @brief Create method
     * @param pNewLinearSystemSolver The linear solver for the system of equations
     * @param ThisParameters The configuration parameters
     */
    virtual ClassType::Pointer Create(
        typename TLinearSolver::Pointer pNewLinearSystemSolver,
        Parameters ThisParameters
        ) const
    {
        return Kratos::make_shared<ClassType>(pNewLinearSystemSolver,ThisParameters);
    }

    ///@}
    ///@name Operators
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
    void SetCalculateReactionsFlag(bool flag)
    {
        mCalculateReactionsFlag = flag;
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
     * @return The flag that tells if we need to reshape the LHS matrix
     */
    bool GetReshapeMatrixFlag() const
    {
        return mReshapeMatrixFlag;
    }

    /**
     * @brief This method sets the flag mReshapeMatrixFlag
     * @param ReshapeMatrixFlag The flag that tells if we need to reshape the LHS matrix
     */
    void SetReshapeMatrixFlag(bool ReshapeMatrixFlag)
    {
        mReshapeMatrixFlag = ReshapeMatrixFlag;
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
     * @brief This method return the linear solver used
     * @return mpLinearSystemSolver The linear solver used
     */
    typename TLinearSolver::Pointer GetLinearSystemSolver() const
    {
        return mpLinearSystemSolver;
    }

    /**
     * @brief This method sets the linear solver to be used
     * @param pLinearSystemSolver The linear solver to be used
     */
    void SetLinearSystemSolver(typename TLinearSolver::Pointer pLinearSystemSolver)
    {
       mpLinearSystemSolver = pLinearSystemSolver;
    }

    /**
     * @brief Function to perform the building of the LHS, depending on the implementation chosen the size of the matrix could be equal to the total number of Dofs or to the number unrestrained dofs
     * @param pScheme The pointer to the integration scheme
     * @param rModelPart The model part to compute
     * @param rA The LHS matrix of the system of equations
     */
    virtual void BuildLHS(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& rA
        )
    {
    }

    /**
     * @brief Function to perform the build of the RHS. The vector could be sized as the total number of dofs or as the number of unrestrained ones
     * @param pScheme The pointer to the integration scheme
     * @param rModelPart The model part to compute
     * @param rb The RHS vector of the system of equations
     */
    virtual void BuildRHS(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemVectorType& rb
        )
    {
    }

    virtual SizeType GetNumberOfROMModes() const noexcept
    {
    }

    virtual void ProjectRHS_ROM(
        ModelPart &rModelPart,
        TSystemVectorType &rb,
        TSystemVectorType &rb_rom)
    {
    }


    virtual void BuildRHSNoDirichlet(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemVectorType& rb)
    {

    }

    /**
     * @brief equivalent (but generally faster) then performing BuildLHS and BuildRHS
     * @param pScheme The pointer to the integration scheme
     * @param rModelPart The model part to compute
     * @param rA The LHS matrix of the system of equations
     * @param rb The RHS vector of the system of equations
     */
    virtual void Build(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rb
        )
    {
    }

    /**
     * @brief It builds a rectangular matrix of size n*N where "n" is the number of unrestrained degrees of freedom and "N" is the total number of degrees of freedom involved.
     * @details This matrix is obtained by building the total matrix without the lines corresponding to the fixed degrees of freedom (but keeping the columns!!)
     * @param pScheme The pointer to the integration scheme
     * @param rModelPart The model part to compute
     * @param rA The LHS matrix of the system of equations
     */
    virtual void BuildLHS_CompleteOnFreeRows(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& rA
        )
    {
    }

    /**
     * @brief It builds a  matrix of size N*N where "N" is the total number of degrees of freedom involved.
     * @details This matrix is obtained by building the total matrix including the lines and columns corresponding to the fixed degrees of freedom
     * @param pScheme The pointer to the integration scheme
     * @param rModelPart The model part to compute
     * @param rA The LHS matrix of the system of equations
     */
    virtual void BuildLHS_Complete(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& rA
        )
    {
    }

    /**
     * @brief This is a call to the linear system solver
     * @param rA The LHS matrix of the system of equations
     * @param rDx The vector of unknowns
     * @param rb The RHS vector of the system of equations
     */
    virtual void SystemSolve(
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
    )
    {
    }

    /**
     * @brief Function to perform the building and solving phase at the same time.
     * @details It is ideally the fastest and safer function to use when it is possible to solve just after building
     * @param pScheme The pointer to the integration scheme
     * @param rModelPart The model part to compute
     * @param rA The LHS matrix of the system of equations
     * @param rDx The vector of unknowns
     * @param rb The RHS vector of the system of equations
     */
    virtual void BuildAndSolve(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb)
    {
    }

    /**
     * @brief Function to perform the building and solving phase at the same time Linearizing with the database at the old iteration
     * @details It is ideally the fastest and safer function to use when it is possible to solve just after building
     * @param pScheme The pointer to the integration scheme
     * @param rModelPart The model part to compute
     * @param rA The LHS matrix of the system of equations
     * @param rDx The vector of unknowns
     * @param rb The RHS vector of the system of equations
     * @param MoveMesh tells if the update of the scheme needs  to be performed when calling the Update of the scheme
     */
    virtual void BuildAndSolveLinearizedOnPreviousIteration(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb,
        const bool MoveMesh
        )
    {
        KRATOS_ERROR << "No special implementation available for "
            << "BuildAndSolveLinearizedOnPreviousIteration "
            << " please use UseOldStiffnessInFirstIterationFlag=false in the settings of the strategy "
            << std::endl;
    }

    /**
     * @brief Corresponds to the previews, but the System's matrix is considered already built and only the RHS is built again
     * @param pScheme The pointer to the integration scheme
     * @param rModelPart The model part to compute
     * @param rA The LHS matrix of the system of equations
     * @param rDx The vector of unknowns
     * @param rb The RHS vector of the system of equations
     */
    virtual void BuildRHSAndSolve(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
        )
    {
    }

    /**
     * @brief It applies the dirichlet conditions. This operation may be very heavy or completely unexpensive depending on the implementation chosen and on how the System Matrix is built.
     * @details For explanation of how it works for a particular implementation the user should refer to the particular Builder And Solver chosen
     * @param pScheme The pointer to the integration scheme
     * @param rModelPart The model part to compute
     * @param rA The LHS matrix of the system of equations
     * @param rDx The vector of unknowns
     * @param rb The RHS vector of the system of equations
     */
    virtual void ApplyDirichletConditions(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
        )
    {
    }

    /**
     * @brief The same of the precedent but affecting only the LHS
     * @param pScheme The pointer to the integration scheme
     * @param rModelPart The model part to compute
     * @param rA The LHS matrix of the system of equations
     * @param rDx The vector of unknowns
     */
    virtual void ApplyDirichletConditions_LHS(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx
        )
    {
    }

    /**
     * @brief The same of the precedent but affecting only the RHS
     * @param pScheme The pointer to the integration scheme
     * @param rModelPart The model part to compute
     * @param rA The LHS matrix of the system of equations
     * @param rb The RHS vector of the system of equations
     */
    virtual void ApplyDirichletConditions_RHS(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
        )
    {
    }

    /**
     * @brief Applies the constraints with master-slave relation matrix (RHS only)
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param rb The RHS vector
     */
    virtual void ApplyRHSConstraints(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemVectorType& rb
        )
    {
    }

    /**
     * @brief Applies the constraints with master-slave relation matrix
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param rA The LHS matrix
     * @param rb The RHS vector
     */
    virtual void ApplyConstraints(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rb
        )
    {
    }

    /**
     * @brief Builds the list of the DofSets involved in the problem by "asking" to each element and condition its Dofs.
     * @details The list of dofs is stores inside the BuilderAndSolver as it is closely connected to the way the matrix and RHS are built
     * @param pScheme The pointer to the integration scheme
     * @param rModelPart The model part to compute
     */
    virtual void SetUpDofSet(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart
        )
    {
    }

    /**
     * @brief It allows to get the list of Dofs from the element
     */
    virtual DofsArrayType& GetDofSet()
    {
        return mDofSet;
    }

    /**
     * @brief It allows to get the list of Dofs from the element
     */
    virtual const DofsArrayType& GetDofSet() const
    {
        return mDofSet;
    }

    /**
     * @brief It organises the dofset in order to speed up the building phase
     * @param rModelPart The model part to compute
     */
    virtual void SetUpSystem(ModelPart& rModelPart)
    {
    }

    /**
     * @brief This method initializes and resizes the system of equations
     * @param pScheme The pointer to the integration scheme
     * @param rModelPart The model part to compute
     * @param pA The pointer to LHS matrix of the system of equations
     * @param pDx The pointer to  vector of unknowns
     * @param pb The pointer to  RHS vector of the system of equations
     */
    virtual void ResizeAndInitializeVectors(
        typename TSchemeType::Pointer pScheme,
        TSystemMatrixPointerType& pA,
        TSystemVectorPointerType& pDx,
        TSystemVectorPointerType& pb,
        ModelPart& rModelPart
        )
    {
    }

    /**
     * @brief It applies certain operations at the system of equations at the beginning of the solution step
     * @param rModelPart The model part to compute
     * @param rA The LHS matrix of the system of equations
     * @param rDx The vector of unknowns
     * @param rb The RHS vector of the system of equations
     */
    virtual void InitializeSolutionStep(
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
        )
    {
    }

    /**
     * @brief It applies certain operations at the system of equations at the end of the solution step
     * @param rModelPart The model part to compute
     * @param rA The LHS matrix of the system of equations
     * @param rDx The vector of unknowns
     * @param rb The RHS vector of the system of equations
     */
    virtual void FinalizeSolutionStep(
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
        )
    {
    }

    /**
     * @brief It computes the reactions of the system
     * @param pScheme The pointer to the integration scheme
     * @param rModelPart The model part to compute
     * @param rA The LHS matrix of the system of equations
     * @param rDx The vector of unknowns
     * @param rb The RHS vector of the system of equations
     */
    virtual void CalculateReactions(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
        )
    {
    }

    /**
     * @brief This function is intended to be called at the end of the solution step to clean up memory
    storage not needed
     */
    virtual void Clear()
    {
        this->mDofSet = DofsArrayType();
        this->mpReactionsVector.reset();
        if (this->mpLinearSystemSolver != nullptr) this->mpLinearSystemSolver->Clear();

        KRATOS_INFO_IF("BuilderAndSolver", this->GetEchoLevel() > 0) << "Clear Function called" << std::endl;
    }

    /**
     * @brief This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param rModelPart The model part to compute
     * @return 0 all ok
     */
    virtual int Check(ModelPart& rModelPart)
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
            "name"       : "builder_and_solver",
            "echo_level" : 1
        })" );
        return default_parameters;
    }

    /**
     * @brief Returns the name of the class as used in the settings (snake_case format)
     * @return The name of the class
     */
    static std::string Name()
    {
        return "builder_and_solver";
    }

    ///@}
    ///@name Operations
    ///@{

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

    /**
     * @brief This method returns constraint relation (T) matrix
     * @return The constraint relation (T) matrix
     */
    virtual typename TSparseSpace::MatrixType& GetConstraintRelationMatrix()
    {
        KRATOS_ERROR << "GetConstraintRelationMatrix is not implemented in base BuilderAndSolver" << std::endl;
    }

    /**
     * @brief This method returns constraint constant vector
     * @return The constraint constant vector
     */
    virtual typename TSparseSpace::VectorType& GetConstraintConstantVector()
    {
        KRATOS_ERROR << "GetConstraintConstantVector is not implemented in base BuilderAndSolver" << std::endl;
    }

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "BuilderAndSolver";
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

    typename TLinearSolver::Pointer mpLinearSystemSolver = nullptr; /// Pointer to the linear solver

    DofsArrayType mDofSet; /// The set containing the DoF of the system

    bool mReshapeMatrixFlag = false;  /// If the matrix is reshaped each step

    bool mDofSetIsInitialized = false; /// Flag taking care if the dof set was initialized ot not

    bool mCalculateReactionsFlag = false; /// Flag taking in account if it is needed or not to calculate the reactions

    unsigned int mEquationSystemSize; /// Number of degrees of freedom of the problem to be solve

    int mEchoLevel = 0;

    TSystemVectorPointerType mpReactionsVector;

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

}; /* Class BuilderAndSolver */

///@}

///@name Type Definitions
///@{


///@}

} /* namespace Kratos.*/
