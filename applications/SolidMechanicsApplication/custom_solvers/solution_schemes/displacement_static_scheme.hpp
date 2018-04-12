//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:            November 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_DISPLACEMENT_STATIC_SCHEME_H_INCLUDED)
#define  KRATOS_DISPLACEMENT_STATIC_SCHEME_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_solvers/solution_schemes/solution_scheme.hpp"

#include "custom_solvers/time_integration_methods/static_method.hpp"

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

  /** @brief Static integration scheme (for static problems)
   */
  template<class TSparseSpace,  class TDenseSpace >
  class DisplacementStaticScheme: public SolutionScheme<TSparseSpace,TDenseSpace>
  {   
  public:
    
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION( DisplacementStaticScheme );

    typedef SolutionScheme<TSparseSpace,TDenseSpace>                             BaseType;
    typedef typename BaseType::SolutionSchemePointerType                  BasePointerType;

    typedef typename BaseType::NodeType                                          NodeType;
    typedef typename BaseType::DofsArrayType                                DofsArrayType;
    typedef typename Element::DofsVectorType                               DofsVectorType;
    typedef typename BaseType::SystemMatrixType                          SystemMatrixType;
    typedef typename BaseType::SystemVectorType                          SystemVectorType;
    typedef typename BaseType::LocalSystemVectorType                LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType                LocalSystemMatrixType;

    typedef ModelPart::NodesContainerType                              NodesContainerType;
    typedef ModelPart::ElementsContainerType                        ElementsContainerType;  
    typedef ModelPart::ConditionsContainerType                    ConditionsContainerType;

    typedef typename BaseType::IntegrationType                            IntegrationType;
    typedef typename BaseType::IntegrationPointerType              IntegrationPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default Constructor.
    DisplacementStaticScheme()
      :BaseType()
    {
    }

    /// Default Constructor.
    DisplacementStaticScheme(Flags& rOptions)
      :BaseType(rOptions)
    {
    }
    
    /// Copy Constructor.
    DisplacementStaticScheme(DisplacementStaticScheme& rOther)
      :BaseType(rOther)
    {
    }

    /// Clone.
    BasePointerType Clone() override
    {
      return BasePointerType( new DisplacementStaticScheme(*this) );
    }

    /// Destructor.
    ~DisplacementStaticScheme() override {}

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
      NodesContainerType::iterator NodeBegin = rModelPart.Nodes().begin();

#pragma omp parallel for firstprivate(NodeBegin)
      for(int i = 0;  i < nnodes; i++)
        {
	  NodesContainerType::iterator itNode = NodeBegin + i;

	  this->IntegrationMethodUpdate(*itNode);
        }

      this->MoveMesh(rModelPart);
      
      KRATOS_CATCH( "" );
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

      // Updating time derivatives (nodally for efficiency)
      const unsigned int NumThreads = OpenMPUtils::GetNumThreads();
      OpenMPUtils::PartitionVector NodePartition;
      OpenMPUtils::DivideInPartitions(rModelPart.Nodes().size(), NumThreads, NodePartition);

      const int nnodes = static_cast<int>( rModelPart.Nodes().size() );
      NodesContainerType::iterator NodeBegin = rModelPart.Nodes().begin();

#pragma omp parallel for firstprivate(NodeBegin)
      for(int i = 0;  i< nnodes; i++)
        {
	  NodesContainerType::iterator itNode = NodeBegin + i;

	  this->IntegrationMethodPredict(*itNode);
        }

      this->MoveMesh(rModelPart);
      
      KRATOS_CATCH( "" );
    }

    /**
     * This is the place to initialize the elements.
     * This is intended to be called just once when the strategy is initialized
     * @param rModelPart: The model of the problem to solve
     */
    void InitializeElements(ModelPart& rModelPart) override
    {
      KRATOS_TRY;

      BaseType::InitializeElements(rModelPart);

      KRATOS_CATCH( "" );
    }

    /**
     * This is the place to initialize the conditions. This is intended to be called just once when the strategy is initialized
     * @param rModelPart: The model of the problem to solve
     */

    void InitializeConditions(ModelPart& rModelPart) override
    {
      KRATOS_TRY;

      BaseType::InitializeConditions(rModelPart);
      
      KRATOS_CATCH( "" );
    }

    /**
     * It initializes time step solution. Only for reasons if the time step solution is restarted
     * @param rModelPart: The model of the problem to solve
      *
     */
    
    void InitializeSolutionStep(ModelPart& rModelPart) override
    {
      KRATOS_TRY;

      BaseType::InitializeSolutionStep(rModelPart);

      KRATOS_CATCH( "" );
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

      KRATOS_CATCH( "" );
    }

    /**
     * It initializes a non-linear iteration
     * @param rModelPart: The model of the problem to solve
     */

    void InitializeNonLinearIteration(ModelPart& rModelPart) override
    {
      KRATOS_TRY;

      BaseType::InitializeNonLinearIteration(rModelPart);

      KRATOS_CATCH( "" );
    }


    /**
     * It finalizes a non-linear iteration
     * @param rModelPart: The model of the problem to solve
     */

    void FinalizeNonLinearIteration(ModelPart& rModelPart) override
    {
      KRATOS_TRY;

      BaseType::FinalizeNonLinearIteration(rModelPart);

      KRATOS_CATCH( "" );
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
     * @param rLHS_Contribution: The LHS matrix contribution
     * @param rRHS_Contribution: The RHS vector contribution
     * @param rEquationId: The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo: The current process info instance
     */

    void CalculateSystemContributions(Element::Pointer rCurrentElement,
				      LocalSystemMatrixType& rLHS_Contribution,
				      LocalSystemVectorType& rRHS_Contribution,
				      Element::EquationIdVectorType& rEquationId,
				      ProcessInfo& rCurrentProcessInfo) override
    {
      KRATOS_TRY;

      (rCurrentElement) -> CalculateLocalSystem(rLHS_Contribution,rRHS_Contribution, rCurrentProcessInfo);

      (rCurrentElement) -> EquationIdVector(rEquationId, rCurrentProcessInfo);

      KRATOS_CATCH( "" );
    }

    /**
     * This function is designed to calculate just the RHS contribution
     * @param rCurrentElemen: The element to compute
     * @param rRHS_Contribution: The RHS vector contribution
     * @param rEquationId: The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo: The current process info instance
     */

    void Calculate_RHS_Contribution(Element::Pointer rCurrentElement,
				    LocalSystemVectorType& rRHS_Contribution,
				    Element::EquationIdVectorType& rEquationId,
				    ProcessInfo& rCurrentProcessInfo) override
    {

      KRATOS_TRY;

      // Basic operations for the element considered
      (rCurrentElement) -> CalculateRightHandSide(rRHS_Contribution, rCurrentProcessInfo);

      (rCurrentElement) -> EquationIdVector(rEquationId, rCurrentProcessInfo);

      KRATOS_CATCH( "" );
    }

    /**
     * This function is designed to calculate just the LHS contribution
     * @param rCurrentElemen: The element to compute
     * @param rLHS_Contribution: The LHS matrix contribution
     * @param rEquationId: The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo: The current process info instance
     */

    void Calculate_LHS_Contribution(Element::Pointer rCurrentElement,
				    LocalSystemMatrixType& rLHS_Contribution,
				    Element::EquationIdVectorType& rEquationId,
				    ProcessInfo& rCurrentProcessInfo) override
    {

      KRATOS_TRY;

      // Basic operations for the element considered
      (rCurrentElement) -> CalculateLeftHandSide(rLHS_Contribution, rCurrentProcessInfo);

      (rCurrentElement) -> EquationIdVector(rEquationId, rCurrentProcessInfo);

      KRATOS_CATCH( "" );
    }


    /**
     * Functions totally analogous to the precedent but applied to the "condition" objects
     * @param rCurrentCondition: The condition to compute
     * @param rLHS_Contribution: The LHS matrix contribution
     * @param rRHS_Contribution: The RHS vector contribution
     * @param rEquationId: The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo: The current process info instance
     */

    void Condition_CalculateSystemContributions(Condition::Pointer rCurrentCondition,
						LocalSystemMatrixType& rLHS_Contribution,
						LocalSystemVectorType& rRHS_Contribution,
						Element::EquationIdVectorType& rEquationId,
						ProcessInfo& rCurrentProcessInfo) override
    {
      KRATOS_TRY;

      // Basic operations for the condition considered
      (rCurrentCondition) -> CalculateLocalSystem(rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);

      (rCurrentCondition) -> EquationIdVector(rEquationId, rCurrentProcessInfo);

      KRATOS_CATCH( "" );
    }

    /**
     * Functions that calculates the RHS of a "condition" object
     * @param rCurrentCondition: The condition to compute
     * @param rRHS_Contribution: The RHS vector contribution
     * @param rEquationId: The ID's of the condition degrees of freedom
     * @param rCurrentProcessInfo: The current process info instance
     */

    void Condition_Calculate_RHS_Contribution(Condition::Pointer rCurrentCondition,
					      LocalSystemVectorType& rRHS_Contribution,
					      Element::EquationIdVectorType& rEquationId,
					      ProcessInfo& rCurrentProcessInfo) override
    {
      KRATOS_TRY;

      // Basic operations for the condition considered
      (rCurrentCondition) -> CalculateRightHandSide(rRHS_Contribution, rCurrentProcessInfo);

      (rCurrentCondition) -> EquationIdVector(rEquationId, rCurrentProcessInfo);

      KRATOS_CATCH( "" );
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

      // Check that variables are correctly allocated
      for(ModelPart::NodesContainerType::iterator it=rModelPart.NodesBegin(); it!=rModelPart.NodesEnd(); ++it)
        {
	  // Nodal data
	  KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT,(*it));

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
      
      if ( mpIntegrationMethod == NULL ) {
         ProcessInfo & rCurrentProcessInfo = rModelPart.GetProcessInfo();
         this->SetIntegrationMethod( rCurrentProcessInfo);
      }
      if ( mpIntegrationMethod == NULL ) {
	      KRATOS_ERROR << "scheme do not have a Time Integration Method " << mpIntegrationMethod << std::endl;
      }


      return ErrorCode;
      
      KRATOS_CATCH( "" );
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
        buffer << "Displacement StaticScheme";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Displacement StaticScheme";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "Displacement StaticScheme Data";     
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

    IntegrationPointerType    mpIntegrationMethod;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    virtual void SetIntegrationMethod(ProcessInfo& rCurrentProcessInfo)
    {
      this->mpIntegrationMethod = IntegrationPointerType( new StaticMethod<Variable<array_1d<double, 3> >, array_1d<double,3> > );

      // Set scheme variables
      this->mpIntegrationMethod->SetVariable(DISPLACEMENT);
      
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
  }; // Class DisplacementStaticScheme
  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{

  
  ///@}

  ///@} addtogroup block
  
}  // namespace Kratos.

#endif // KRATOS_DISPLACEMENT_STATIC_SCHEME_H_INCLUDED defined
