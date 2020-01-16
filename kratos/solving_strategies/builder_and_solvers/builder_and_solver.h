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

#if !defined(KRATOS_BUILDER_AND_SOLVER )
#define  KRATOS_BUILDER_AND_SOLVER

/* System includes */
#include <set>

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"

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

    /// The definition of the DoF objects
    typedef typename DofsArrayType::iterator DofIterator;
    typedef typename DofsArrayType::const_iterator DofConstantIterator;

    /// The containers of the entities
    typedef ModelPart::NodesContainerType NodesArrayType;
    typedef ModelPart::ElementsContainerType ElementsArrayType;
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;

    /// The definition of the element container type
    typedef PointerVectorSet<Element, IndexedObject> ElementsContainerType;

    /// Pointer definition of BuilderAndSolver
    KRATOS_CLASS_POINTER_DEFINITION(BuilderAndSolver);

    /**
     * @brief This struct is used in the component wise calculation only is defined here and is used to declare a member variable in the component wise builder and solver
     * @details Private pointers can only be accessed by means of set and get functions this allows to set and not copy the Element_Variables and Condition_Variables which will be asked and set by another strategy object
     */
    struct GlobalSystemComponents
    {
    private:

      // Elements
      std::vector<TSystemMatrixType> *mpLHS_Element_Components;
      const std::vector< Variable< LocalSystemMatrixType > > *mpLHS_Element_Variables;

      std::vector<TSystemVectorType> *mpRHS_Element_Components;
      const std::vector< Variable< LocalSystemVectorType > > *mpRHS_Element_Variables;

      // Conditions
      std::vector<TSystemMatrixType> *mpLHS_Condition_Components;
      const std::vector< Variable< LocalSystemMatrixType > > *mpLHS_Condition_Variables;

      std::vector<TSystemVectorType> *mpRHS_Condition_Components;
      const std::vector< Variable< LocalSystemVectorType > > *mpRHS_Condition_Variables;

    public:

      void Initialize()
      {
          mpLHS_Element_Components = NULL;
          mpLHS_Element_Variables  = NULL;

          mpRHS_Element_Components = NULL;
          mpRHS_Element_Variables  = NULL;

          mpLHS_Condition_Components = NULL;
          mpLHS_Condition_Variables  = NULL;

          mpRHS_Condition_Components = NULL;
          mpRHS_Condition_Variables  = NULL;
      }

      // Setting pointer variables

      // Elements
      void SetLHS_Element_Components ( std::vector<TSystemMatrixType>& rLHS_Element_Components ) { mpLHS_Element_Components = &rLHS_Element_Components; };
      void SetLHS_Element_Variables     ( const std::vector< Variable< LocalSystemMatrixType > >& rLHS_Element_Variables ) { mpLHS_Element_Variables = &rLHS_Element_Variables; };
      void SetRHS_Element_Components ( std::vector<TSystemVectorType>& rRHS_Element_Components ) { mpRHS_Element_Components = &rRHS_Element_Components; };
      void SetRHS_Element_Variables     ( const std::vector< Variable< LocalSystemVectorType > >& rRHS_Element_Variables ) { mpRHS_Element_Variables = &rRHS_Element_Variables; };

      bool Are_LHS_Element_Components_Set() { if( mpLHS_Element_Variables == NULL ) return false; else return true; };
      bool Are_RHS_Element_Components_Set() { if( mpRHS_Element_Variables == NULL ) return false; else return true; };

      // Conditions
      void SetLHS_Condition_Components ( std::vector<TSystemMatrixType>& rLHS_Condition_Components ) { mpLHS_Condition_Components = &rLHS_Condition_Components; };
      void SetLHS_Condition_Variables     ( const std::vector< Variable< LocalSystemMatrixType > >& rLHS_Condition_Variables ) { mpLHS_Condition_Variables = &rLHS_Condition_Variables; };
      void SetRHS_Condition_Components ( std::vector<TSystemVectorType>& rRHS_Condition_Components ) { mpRHS_Condition_Components = &rRHS_Condition_Components; };
      void SetRHS_Condition_Variables     ( const std::vector< Variable< LocalSystemVectorType > >& rRHS_Condition_Variables ) { mpRHS_Condition_Variables = &rRHS_Condition_Variables; };

      bool Are_LHS_Condition_Components_Set() { if( mpLHS_Condition_Variables == NULL ) return false; else return true; };
      bool Are_RHS_Condition_Components_Set() { if( mpRHS_Condition_Variables == NULL ) return false; else return true; };

      //getting pointer variables

      // Elements
      std::vector<TSystemMatrixType>& GetLHS_Element_Components() { return *mpLHS_Element_Components; };
      const std::vector< Variable< LocalSystemMatrixType > >& GetLHS_Element_Variables() { return *mpLHS_Element_Variables; };
      std::vector<TSystemVectorType>& GetRHS_Element_Components() { return *mpRHS_Element_Components; };
      const std::vector< Variable< LocalSystemVectorType > >& GetRHS_Element_Variables() { return *mpRHS_Element_Variables; };

      // Conditions
      std::vector<TSystemMatrixType>& GetLHS_Condition_Components() { return *mpLHS_Condition_Components; };
      const std::vector< Variable< LocalSystemMatrixType > >& GetLHS_Condition_Variables() { return *mpLHS_Condition_Variables; };
      std::vector<TSystemVectorType>& GetRHS_Condition_Components() { return *mpRHS_Condition_Components; };
      const std::vector< Variable< LocalSystemVectorType > >& GetRHS_Condition_Variables() { return *mpRHS_Condition_Variables; };

    };

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor with Parameters
     * @param pNewLinearSystemSolver The linear solver for the system of equations
     * @param ThisParameters The configuration parameters
     */
    explicit BuilderAndSolver(typename TLinearSolver::Pointer pNewLinearSystemSolver, Parameters ThisParameters)
    {
        // Validate default parameters
        Parameters default_parameters = Parameters(R"(
        {
        })" );

        ThisParameters.ValidateAndAssignDefaults(default_parameters);

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


    ///@}
    ///@name Operators
    ///@{

    /**
     * Component wise components Get method
     */
    virtual GlobalSystemComponents& GetGlobalSystemComponents()
    {
        KRATOS_ERROR <<  "Asking for Global Components to the BUIDER and SOlVER base class which is not component wise and not contains this member variable" << std::endl;
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
    bool GetDofSetIsInitializedFlag()
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
    bool GetReshapeMatrixFlag()
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
    unsigned int GetEquationSystemSize()
    {
        return mEquationSystemSize;
    }

    /**
     * @brief This method return the linear solver used
     * @return mpLinearSystemSolver The linear solver used
     */
    typename TLinearSolver::Pointer GetLinearSystemSolver()
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
     * @brief Function to perform the building of the LHS, depending on the implementation choosen the size of the matrix could be equal to the total number of Dofs or to the number unrestrained dofs
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
     * @param rDx The vector of unkowns
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
     * @param rDx The vector of unkowns
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
     * @brief Corresponds to the previews, but the System's matrix is considered already built and only the RHS is built again
     * @param pScheme The pointer to the integration scheme
     * @param rModelPart The model part to compute
     * @param rA The LHS matrix of the system of equations
     * @param rDx The vector of unkowns
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
     * @brief It applies the dirichlet conditions. This operation may be very heavy or completely unexpensive depending on the implementation choosen and on how the System Matrix is built.
     * @details For explanation of how it works for a particular implementation the user should refer to the particular Builder And Solver choosen
     * @param pScheme The pointer to the integration scheme
     * @param rModelPart The model part to compute
     * @param rA The LHS matrix of the system of equations
     * @param rDx The vector of unkowns
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
     * @param rDx The vector of unkowns
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
     * @brief Applies the constraints
     * @param pScheme The pointer to the integration scheme
     * @param rModelPart The model part to compute
     * @param rb The RHS vector of the system of equations
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
     * @details The list of dofs is stores insde the BuilderAndSolver as it is closely connected to the way the matrix and RHS are built
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
     * @param pDx The pointer to  vector of unkowns
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
     * @brief It applies certain operations at the system of equations at the begining of the solution step
     * @param rModelPart The model part to compute
     * @param rA The LHS matrix of the system of equations
     * @param rDx The vector of unkowns
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
     * @param rDx The vector of unkowns
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
     * @param rDx The vector of unkowns
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
     * @brief It sets the level of echo for the solving strategy
     * @param Level The level to set
     * @details The different levels of echo are:
     * - 0: Mute... no echo at all
     * - 1: Printing time and basic informations
     * - 2: Printing linear solver data
     * - 3: Print of debug informations: Echo of stiffness matrix, Dx, b...
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
    int GetEchoLevel()
    {
        return mEchoLevel;
    }

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

    typename TLinearSolver::Pointer mpLinearSystemSolver; /// Pointer to the linear solver

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

#endif /* KRATOS_BUILDER_AND_SOLVER  defined */

