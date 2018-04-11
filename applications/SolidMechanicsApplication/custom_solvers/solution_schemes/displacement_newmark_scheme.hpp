//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:            November 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_DISPLACEMENT_NEWMARK_SCHEME_H_INCLUDED)
#define  KRATOS_DISPLACEMENT_NEWMARK_SCHEME_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_solvers/solution_schemes/solution_scheme.hpp"
#include "custom_solvers/time_integration_methods/newmark_method.hpp"

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

  /** @brief Newmark integration scheme (for dynamic problems)
   */
  template<class TSparseSpace,  class TDenseSpace >
  class DisplacementNewmarkScheme: public SolutionScheme<TSparseSpace,TDenseSpace>
  {
  protected:

    struct  GeneralMatrices
    {
      std::vector< Matrix > M;     // First derivative matrix  (usually mass matrix)
      std::vector< Matrix > D;     // Second derivative matrix (usually damping matrix)
    };

    struct GeneralVectors
    {
      std::vector< Vector > v;    // Velocity
      std::vector< Vector > a;    // Acceleration
      std::vector< Vector > ap;   // Previous acceleration
    };

    
  public:
    
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION( DisplacementNewmarkScheme );
    
    typedef SolutionScheme<TSparseSpace,TDenseSpace>                             BaseType;
    typedef typename BaseType::SolutionSchemePointerType                  BasePointerType;
     typedef typename BaseType::LocalFlagType                               LocalFlagType;

    typedef typename BaseType::NodeType                                          NodeType;
    typedef typename BaseType::DofsArrayType                                DofsArrayType;
    typedef typename BaseType::SystemMatrixType                          SystemMatrixType;
    typedef typename BaseType::SystemVectorType                          SystemVectorType;
    typedef typename BaseType::LocalSystemVectorType                LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType                LocalSystemMatrixType; 

    typedef typename BaseType::NodesContainerType                      NodesContainerType;
    typedef typename BaseType::ElementsContainerType                ElementsContainerType;
    typedef typename BaseType::ConditionsContainerType            ConditionsContainerType;
    typedef typename BaseType::IntegrationPointerType              IntegrationPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default Constructor.
    DisplacementNewmarkScheme()
      :BaseType()
    {
    }

    /// Constructor.
    DisplacementNewmarkScheme(Flags& rOptions)
      :BaseType(rOptions)
    {
    }
    
    /// Copy Constructor.
    DisplacementNewmarkScheme(DisplacementNewmarkScheme& rOther)
      :BaseType(rOther)
      ,mMatrix(rOther.mMatrix)
      ,mVector(rOther.mVector)
    {
    }

    /// Clone.
    BasePointerType Clone() override
    {
      return BasePointerType( new DisplacementNewmarkScheme(*this) );
    }

    /// Destructor.
    ~DisplacementNewmarkScheme() override {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
    this is the place to initialize the Scheme.
    This is intended to be called just once when the strategy is initialized
     */
    virtual void Initialize(ModelPart& rModelPart) override
    {
        KRATOS_TRY

	BaseType::Initialize(rModelPart);
	  
	ProcessInfo& rCurrentProcessInfo= rModelPart.GetProcessInfo();

	// Set integration method
	this->SetIntegrationMethod(rCurrentProcessInfo);
            
	// Allocate auxiliary memory
	const unsigned int NumThreads = OpenMPUtils::GetNumThreads();

	mMatrix.M.resize(NumThreads);
	mMatrix.D.resize(NumThreads);
	
	mVector.v.resize(NumThreads);
	mVector.a.resize(NumThreads);

	KRATOS_CATCH("")
    }

    
    /**
     * Performing the update of the solution
     * Incremental update within newton iteration. It updates the state variables at the end of the time step: u_{n+1}^{k+1}= u_{n+1}^{k}+ \Delta u
     * @param rModelPart: The model of the problem to solve
     * @param rDofSet: Set of all primary variables
     * @param rDx: incremental update of primary variables
      */

    void Update(ModelPart& rModelPart,
		DofsArrayType& rDofSet,
		SystemVectorType& rDx) override
    {
      KRATOS_TRY;

      const unsigned int NumThreads = OpenMPUtils::GetNumThreads();

      // Update of displacement (by DOF)
      OpenMPUtils::PartitionVector DofPartition;
      OpenMPUtils::DivideInPartitions(rDofSet.size(), NumThreads, DofPartition);

      const int ndof = static_cast<int>(rDofSet.size());
      typename DofsArrayType::iterator DofBegin = rDofSet.begin();

#pragma omp parallel for firstprivate(DofBegin)
      for(int i = 0;  i < ndof; i++)
        {
	  typename DofsArrayType::iterator itDof = DofBegin + i;

	  if (itDof->IsFree() )
            {
	      itDof->GetSolutionStepValue() += TSparseSpace::GetValue(rDx,itDof->EquationId());
            }
        }

      // Updating time derivatives (nodally for efficiency)
      OpenMPUtils::PartitionVector NodePartition;
      OpenMPUtils::DivideInPartitions(rModelPart.Nodes().size(), NumThreads, NodePartition);

      const int nnodes = static_cast<int>(rModelPart.Nodes().size());
      typename NodesContainerType::iterator NodeBegin = rModelPart.Nodes().begin();

#pragma omp parallel for firstprivate(NodeBegin)
      for(int i = 0;  i < nnodes; i++)
        {
	  typename NodesContainerType::iterator itNode = NodeBegin + i;

	  this->IntegrationMethodUpdate(*itNode);
        }

      this->MoveMesh(rModelPart);
      
      KRATOS_CATCH("")
    }

    /**
     * Performing the prediction of the solution
     * It predicts the solution for the current step: x = xold + vold * Dt
     * @param rModelPart: The model of the problem to solve
     * @param rDofSet set of all primary variables
     * @param rDx: Incremental update of primary variables
     */

    void Predict(ModelPart& rModelPart,
		 DofsArrayType& rDofSet,
		 SystemVectorType& rDx) override
    {
      KRATOS_TRY;

      // std::cout << " Prediction " << std::endl;

      // Updating time derivatives (nodally for efficiency)
      const unsigned int NumThreads = OpenMPUtils::GetNumThreads();
      OpenMPUtils::PartitionVector NodePartition;
      OpenMPUtils::DivideInPartitions(rModelPart.Nodes().size(), NumThreads, NodePartition);

      const int nnodes = static_cast<int>( rModelPart.Nodes().size() );
      typename NodesContainerType::iterator NodeBegin = rModelPart.Nodes().begin();

#pragma omp parallel for firstprivate(NodeBegin)
      for(int i = 0;  i< nnodes; i++)
        {
	  typename NodesContainerType::iterator itNode = NodeBegin + i;

	  this->IntegrationMethodPredict(*itNode);
        }

      this->MoveMesh(rModelPart);    

      KRATOS_CATCH("")
    }

    /**
     * This is the place to initialize the elements.
     * This is intended to be called just once when the strategy is initialized
     * @param rModelPart: The model of the problem to solve
     */
    void InitializeElements(ModelPart& rModelPart) override
    {
      KRATOS_TRY;

      const unsigned int NumThreads = OpenMPUtils::GetNumThreads();
      OpenMPUtils::PartitionVector ElementPartition;
      OpenMPUtils::DivideInPartitions(rModelPart.Elements().size(), NumThreads, ElementPartition);

      const int nelem = static_cast<int>(rModelPart.Elements().size());
      typename ElementsContainerType::iterator ElemBegin = rModelPart.Elements().begin();

#pragma omp parallel for
      for(int i = 0;  i < nelem; i++)
        {
	  typename ElementsContainerType::iterator itElem = ElemBegin + i;

	  itElem->Initialize(); //function to initialize the element
        }

      this->Set(LocalFlagType::ELEMENTS_INITIALIZED,true);

      KRATOS_CATCH("")
    }

    /**
     * This is the place to initialize the conditions. This is intended to be called just once when the strategy is initialized
     * @param rModelPart: The model of the problem to solve
     */

    void InitializeConditions(ModelPart& rModelPart) override
    {
      KRATOS_TRY;

      if(this->IsNot(LocalFlagType::ELEMENTS_INITIALIZED))
        {
	  KRATOS_ERROR << "Before initilizing Conditions, initialize Elements FIRST";
        }

      if(this->IsNot(LocalFlagType::CONDITIONS_INITIALIZED))
      {
      
        const unsigned int NumThreads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector ConditionPartition;
        OpenMPUtils::DivideInPartitions(rModelPart.Conditions().size(), NumThreads, ConditionPartition);

        const int ncond = static_cast<int>(rModelPart.Conditions().size());
        typename ConditionsContainerType::iterator CondBegin = rModelPart.Conditions().begin();

#pragma omp parallel for
        for(int i = 0;  i < ncond; i++)
        {
	  typename ConditionsContainerType::iterator itCond = CondBegin + i;

	  itCond->Initialize(); //function to initialize the condition
        }
      }
     
      this->Set(LocalFlagType::CONDITIONS_INITIALIZED, true);

      KRATOS_CATCH("")
    }

    /**
     * It initializes time step solution. Only for reasons if the time step solution is restarted
     * @param rModelPart: The model of the problem to solve
     * @param A: LHS matrix
     * @param Dx: Incremental update of primary variables
     * @param b: RHS Vector
     *
     */
    
    void InitializeSolutionStep(ModelPart& rModelPart) override
    {
      KRATOS_TRY;
      
      BaseType::InitializeSolutionStep(rModelPart);

      KRATOS_CATCH("")
    }

    /**
     * Function called once at the end of a solution step, after convergence is reached if
     * an iterative process is needed
     * @param rModelPart: The model of the problem to solve
     */

    void FinalizeSolutionStep(ModelPart& rModelPart) override
    {
      KRATOS_TRY;

      BaseType::FinalizeSolutionStep(rModelPart);

      KRATOS_CATCH("")
    }

    /**
     * It initializes a non-linear iteration
     * @param rModelPart: The model of the problem to solve
      */

    void InitializeNonLinearIteration(ModelPart& rModelPart) override
    {
      KRATOS_TRY;

      BaseType::InitializeNonLinearIteration(rModelPart);

      KRATOS_CATCH("")
    }

    /**
     * It finalizes a non-linear iteration
     * @param rModelPart: The model of the problem to solve
      */

    void FinalizeNonLinearIteration(ModelPart& rModelPart) override
    {
      KRATOS_TRY;

      BaseType::FinalizeNonLinearIteration(rModelPart);

      KRATOS_CATCH("")
    }
    
    /**
     * It initializes a non-linear iteration (for an individual condition)
     * @param rCurrentConditiont: The condition to compute
     * @param rCurrentProcessInfo: The current process info instance
     */

    void InitializeNonLinearIteration(Condition::Pointer rCurrentCondition,
				      ProcessInfo& rCurrentProcessInfo) override
    {
      (rCurrentCondition) -> InitializeNonLinearIteration(rCurrentProcessInfo);
    }

    /**
     * It initializes a non-linear iteration (for an individual element)
     * @param rCurrentElement: The element to compute
     * @param rCurrentProcessInfo: The current process info instance
     */

    void InitializeNonLinearIteration(Element::Pointer rCurrentElement,
				      ProcessInfo& rCurrentProcessInfo) override
    {
      (rCurrentElement) -> InitializeNonLinearIteration(rCurrentProcessInfo);
    }

    /**
     * This function is designed to be called in the builder and solver to introduce
     * @param rCurrentElement: The element to compute
     * @param LHS_Contribution: The LHS matrix contribution
     * @param RHS_Contribution: The RHS vector contribution
     * @param EquationId: The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo: The current process info instance
     */

    void CalculateSystemContributions(Element::Pointer rCurrentElement,
				      LocalSystemMatrixType& LHS_Contribution,
				      LocalSystemVectorType& RHS_Contribution,
				      Element::EquationIdVectorType& EquationId,
				      ProcessInfo& rCurrentProcessInfo) override
    {
      KRATOS_TRY;

      int thread = OpenMPUtils::ThisThread();

      //(rCurrentElement) -> InitializeNonLinearIteration(rCurrentProcessInfo);

      (rCurrentElement) -> CalculateLocalSystem(LHS_Contribution,RHS_Contribution, rCurrentProcessInfo);

      (rCurrentElement) -> EquationIdVector(EquationId, rCurrentProcessInfo);

      (rCurrentElement) -> CalculateMassMatrix(mMatrix.M[thread], rCurrentProcessInfo);

      (rCurrentElement) -> CalculateDampingMatrix(mMatrix.D[thread], rCurrentProcessInfo);

      AddDynamicsToLHS (LHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], rCurrentProcessInfo);

      AddDynamicsToRHS (rCurrentElement, RHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], rCurrentProcessInfo);

      //AssembleTimeSpaceLHS(rCurrentElement, LHS_Contribution, DampMatrix, MassMatrix, rCurrentProcessInfo);

      KRATOS_CATCH("")
    }

    /**
     * This function is designed to calculate just the RHS contribution
     * @param rCurrentElemen: The element to compute
     * @param RHS_Contribution: The RHS vector contribution
     * @param EquationId: The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo: The current process info instance
     */

    void Calculate_RHS_Contribution(Element::Pointer rCurrentElement,
				    LocalSystemVectorType& RHS_Contribution,
				    Element::EquationIdVectorType& EquationId,
				    ProcessInfo& rCurrentProcessInfo) override
    {

      KRATOS_TRY;

      int thread = OpenMPUtils::ThisThread();

      // Initializing the non linear iteration for the current element
      // (rCurrentElement) -> InitializeNonLinearIteration(rCurrentProcessInfo);

      // Basic operations for the element considered
      (rCurrentElement) -> CalculateRightHandSide(RHS_Contribution, rCurrentProcessInfo);

      (rCurrentElement) -> CalculateMassMatrix(mMatrix.M[thread], rCurrentProcessInfo);

      (rCurrentElement) -> CalculateDampingMatrix(mMatrix.D[thread], rCurrentProcessInfo);

      (rCurrentElement) -> EquationIdVector(EquationId, rCurrentProcessInfo);

      AddDynamicsToRHS (rCurrentElement, RHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], rCurrentProcessInfo);

      KRATOS_CATCH("")
    }

    /**
     * Functions totally analogous to the precedent but applied to the "condition" objects
     * @param rCurrentCondition: The condition to compute
     * @param LHS_Contribution: The LHS matrix contribution
     * @param RHS_Contribution: The RHS vector contribution
     * @param EquationId: The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo: The current process info instance
     */

    void Condition_CalculateSystemContributions(Condition::Pointer rCurrentCondition,
						LocalSystemMatrixType& LHS_Contribution,
						LocalSystemVectorType& RHS_Contribution,
						Element::EquationIdVectorType& EquationId,
						ProcessInfo& rCurrentProcessInfo) override
    {
      KRATOS_TRY;

      int thread = OpenMPUtils::ThisThread();

      // Initializing the non linear iteration for the current condition
      //(rCurrentCondition) -> InitializeNonLinearIteration(rCurrentProcessInfo);

      // Basic operations for the condition considered
      (rCurrentCondition) -> CalculateLocalSystem(LHS_Contribution,RHS_Contribution, rCurrentProcessInfo);

      (rCurrentCondition) -> EquationIdVector(EquationId, rCurrentProcessInfo);

      (rCurrentCondition) -> CalculateMassMatrix(mMatrix.M[thread], rCurrentProcessInfo);

      (rCurrentCondition) -> CalculateDampingMatrix(mMatrix.D[thread], rCurrentProcessInfo);

      AddDynamicsToLHS  (LHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], rCurrentProcessInfo);

      AddDynamicsToRHS  (rCurrentCondition, RHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], rCurrentProcessInfo);

      // AssembleTimeSpaceLHS_Condition(rCurrentCondition, LHS_Contribution,DampMatrix, MassMatrix, rCurrentProcessInfo);

      KRATOS_CATCH("")
    }

    /**
     * Functions that calculates the RHS of a "condition" object
     * @param rCurrentCondition: The condition to compute
     * @param RHS_Contribution: The RHS vector contribution
     * @param EquationId: The ID's of the condition degrees of freedom
     * @param rCurrentProcessInfo: The current process info instance
     */

    void Condition_Calculate_RHS_Contribution(Condition::Pointer rCurrentCondition,
					      LocalSystemVectorType& RHS_Contribution,
					      Element::EquationIdVectorType& EquationId,
					      ProcessInfo& rCurrentProcessInfo) override
    {
      KRATOS_TRY;

      int thread = OpenMPUtils::ThisThread();

      // Initializing the non linear iteration for the current condition
      //(rCurrentCondition) -> InitializeNonLinearIteration(rCurrentProcessInfo);

      // Basic operations for the condition considered
      (rCurrentCondition) -> CalculateRightHandSide(RHS_Contribution, rCurrentProcessInfo);

      (rCurrentCondition) -> EquationIdVector(EquationId, rCurrentProcessInfo);

      (rCurrentCondition) -> CalculateMassMatrix(mMatrix.M[thread], rCurrentProcessInfo);

      (rCurrentCondition) -> CalculateDampingMatrix(mMatrix.D[thread], rCurrentProcessInfo);

      // Adding the dynamic contributions (static is already included)
      AddDynamicsToRHS  (rCurrentCondition, RHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], rCurrentProcessInfo);

      KRATOS_CATCH("")
    }

    /**
     * Function that returns the list of Degrees of freedom to be assembled in the system for a Given Element
     * @param rCurrentElement: The element to compute
     * @param rElementalDofsList: The element dofs list 
     * @param rCurrentProcessInfo: The current process info instance
     */

    void GetElementalDofList(Element::Pointer rCurrentElement,
			     Element::DofsVectorType& rElementalDofList,
			     ProcessInfo& rCurrentProcessInfo) override
    {
      rCurrentElement->GetDofList(rElementalDofList, rCurrentProcessInfo);
    }

    /**
     * Function that returns the list of Degrees of freedom to be assembled in the system for a Given Element
     * @param rCurrentCondition: The condition to compute
     * @param rConditionDofsList: The condition dofs list 
     * @param rCurrentProcessInfo: The current process info instance
     */

    void GetConditionDofList(Condition::Pointer rCurrentCondition,
			     Element::DofsVectorType& rConditionDofList,
			     ProcessInfo& rCurrentProcessInfo) override
    {
      rCurrentCondition->GetDofList(rConditionDofList, rCurrentProcessInfo);
    }

    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param rModelPart: The model of the problem to solve
     * @return Zero means  all ok
     */

    virtual int Check(ModelPart& rModelPart) override
    {
      KRATOS_TRY;

      // Perform base base checks
      int ErrorCode = 0;
      ErrorCode  = BaseType::Check(rModelPart);

      // Check that all required variables have been registered
      KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT);
      KRATOS_CHECK_VARIABLE_KEY(VELOCITY);
      KRATOS_CHECK_VARIABLE_KEY(ACCELERATION);

      // Check that variables are correctly allocated
      for(ModelPart::NodesContainerType::iterator it=rModelPart.NodesBegin(); it!=rModelPart.NodesEnd(); ++it)
        {
	  // Nodal data
	  KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT,(*it));
	  KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY,(*it));
	  KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ACCELERATION,(*it));

	  // Nodal dofs
	  KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X,(*it));
	  KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y,(*it));
	  if( rModelPart.GetProcessInfo()[SPACE_DIMENSION] == 3 )
	    KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Z,(*it));
        }

      // Check for minimum value of the buffer index
      if (rModelPart.GetBufferSize() < 2)
        {
	  KRATOS_ERROR << "insufficient buffer size. Buffer size should be greater than 2. Current size is" << rModelPart.GetBufferSize() << std::endl;
        }

      if ( this->mpIntegrationMethod == NULL ) {
         ProcessInfo & rCurrentProcessInfo = rModelPart.GetProcessInfo();
         this->SetIntegrationMethod( rCurrentProcessInfo);
      }
      if ( this->mpIntegrationMethod == NULL ) {
	      KRATOS_ERROR << "scheme do not have a Time Integration Method " << this->mpIntegrationMethod << std::endl;
      }

      return ErrorCode;
      
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
    
    /// Turn back information as a string.
    virtual std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "Displacement NewmarkScheme";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Displacement NewmarkScheme";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "Displacement NewmarkScheme Data";     
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

    GeneralMatrices     mMatrix;

    GeneralVectors      mVector;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    virtual void SetIntegrationMethod(ProcessInfo& rCurrentProcessInfo)
    {      
      this->mpIntegrationMethod = IntegrationPointerType( new NewmarkMethod<Variable<array_1d<double, 3> >, array_1d<double,3> > );

      // Set scheme variables
      this->mpIntegrationMethod->SetVariables(DISPLACEMENT,VELOCITY,ACCELERATION);

      // Set scheme parameters
      this->mpIntegrationMethod->SetParameters(rCurrentProcessInfo);

      // Modify ProcessInfo scheme parameters
      this->mpIntegrationMethod->SetProcessInfoParameters(rCurrentProcessInfo);
    }

    virtual void IntegrationMethodUpdate(NodeType& rNode)
    {
      this->mpIntegrationMethod->Update(rNode);
    }

    virtual void IntegrationMethodPredict(NodeType& rNode)
    {
      this->mpIntegrationMethod->Predict(rNode);
    }


    /**
     * It adds the dynamic LHS contribution of the elements: M*c0 + D*c1 + K
     * @param rLHS_Contribution: The dynamic contribution for the LHS
     * @param rD: The damping matrix
     * @param rM: The mass matrix
     * @param rCurrentProcessInfo: The current process info instance
     */

    virtual void AddDynamicsToLHS(LocalSystemMatrixType& rLHS_Contribution,
				  LocalSystemMatrixType& rD,
				  LocalSystemMatrixType& rM,
				  ProcessInfo& rCurrentProcessInfo)
    {

      double parameter = 0;
      // Adding mass contribution to the dynamic stiffness
      if (rM.size1() != 0) // if M matrix declared
        {
	  parameter = this->mpIntegrationMethod->GetSecondDerivativeParameter(parameter);
	  noalias(rLHS_Contribution) += rM * parameter;
        }

      // Adding  damping contribution
      if (rD.size1() != 0) // if D matrix declared
        {
	  parameter = this->mpIntegrationMethod->GetFirstDerivativeParameter(parameter);
	  noalias(rLHS_Contribution) += rD * parameter;
        }
    }

    /**
     * It adds the dynamic RHS contribution of the elements: b - M*a - D*v
     * @param rCurrentElement: The element to compute
     * @param rRHS_Contribution: The dynamic contribution for the RHS
     * @param rD: The damping matrix
     * @param rM: The mass matrix
     * @param rCurrentProcessInfo: The current process info instance
     */

    virtual void AddDynamicsToRHS(Element::Pointer rCurrentElement,
				  LocalSystemVectorType& rRHS_Contribution,
				  LocalSystemMatrixType& rD,
				  LocalSystemMatrixType& rM,
				  ProcessInfo& rCurrentProcessInfo)
    {
      int thread = OpenMPUtils::ThisThread();

      // Adding inertia contribution
      if (rM.size1() != 0)
        {
	  rCurrentElement->GetSecondDerivativesVector(mVector.a[thread], 0);

	  noalias(rRHS_Contribution) -= prod(rM, mVector.a[thread]);
        }

      // Adding damping contribution
      if (rD.size1() != 0)
        {
	  rCurrentElement->GetFirstDerivativesVector(mVector.v[thread], 0);

	  noalias(rRHS_Contribution) -= prod(rD, mVector.v[thread]);
        }
    }

    /**
     * It adds the dynamic RHS contribution of the condition: b - M*a - D*v
     * @param rCurrentCondition: The condition to compute
     * @param rRHS_Contribution: The dynamic contribution for the RHS
     * @param rD: The damping matrix
     * @param rM: The mass matrix
     * @param rCurrentProcessInfo: The current process info instance
     */

    virtual void AddDynamicsToRHS(Condition::Pointer rCurrentCondition,
				  LocalSystemVectorType& rRHS_Contribution,
				  LocalSystemMatrixType& rD,
				  LocalSystemMatrixType& rM,
				  ProcessInfo& rCurrentProcessInfo)
    {
      int thread = OpenMPUtils::ThisThread();

      // Adding inertia contribution
      if (rM.size1() != 0)
        {
	  rCurrentCondition->GetSecondDerivativesVector(mVector.a[thread], 0);

	  noalias(rRHS_Contribution) -= prod(rM, mVector.a[thread]);
        }

      // Adding damping contribution
      // Damping contribution
      if (rD.size1() != 0)
        {
	  rCurrentCondition->GetFirstDerivativesVector(mVector.v[thread], 0);

	  noalias(rRHS_Contribution) -= prod(rD, mVector.v[thread]);
        }

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
    ///@name Serialization
    ///@{
  
    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{
  
    ///@}
  }; // Class DisplacementNewmarkScheme
  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{

  
  ///@}

  ///@} addtogroup block
  
}  // namespace Kratos.

#endif // KRATOS_DISPLACEMENT_NEWMARK_SCHEME_H_INCLUDED defined
