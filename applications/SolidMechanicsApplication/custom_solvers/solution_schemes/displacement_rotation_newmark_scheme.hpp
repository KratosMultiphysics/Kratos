//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:            November 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_DISPLACEMENT_ROTATION_NEWMARK_SCHEME_H_INCLUDED)
#define  KRATOS_DISPLACEMENT_ROTATION_NEWMARK_SCHEME_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_solvers/solution_schemes/displacement_newmark_scheme.hpp"

#include "custom_solvers/time_integration_methods/newmark_step_rotation_method.hpp"

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
  class DisplacementRotationNewmarkScheme: public DisplacementNewmarkScheme<TSparseSpace,TDenseSpace>
  {   
  public:
    
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION( DisplacementRotationNewmarkScheme );

    typedef SolutionScheme<TSparseSpace,TDenseSpace>                                    BaseType;
    typedef typename BaseType::SolutionSchemePointerType                         BasePointerType;

    typedef typename BaseType::LocalSystemVectorType                       LocalSystemVectorType; 
    typedef typename BaseType::LocalSystemMatrixType                       LocalSystemMatrixType;

    typedef DisplacementNewmarkScheme<TSparseSpace,TDenseSpace>                      DerivedType;

    typedef typename DerivedType::IntegrationPointerType                  IntegrationPointerType;
    
    typedef typename DerivedType::NodeType                                              NodeType;
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default Constructor.
    DisplacementRotationNewmarkScheme()
      :DerivedType()
    {
    }

    /// Constructor.
    DisplacementRotationNewmarkScheme(Flags& rOptions)
      :DerivedType(rOptions)
    {
    }
    
    /// Copy Constructor.
    DisplacementRotationNewmarkScheme(DisplacementRotationNewmarkScheme& rOther)
      :DerivedType(rOther)
      ,mpRotationIntegrationMethod(rOther.mpRotationIntegrationMethod)
    {
    }

    /// Clone.
    BasePointerType Clone() override
    {
      return BasePointerType( new DisplacementRotationNewmarkScheme(*this) );
    }

    /// Destructor.
    ~DisplacementRotationNewmarkScheme() override {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{


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

      int thread = OpenMPUtils::ThisThread();
      
      (rCurrentElement) -> CalculateLocalSystem(rLHS_Contribution,rRHS_Contribution, rCurrentProcessInfo);
      
      (rCurrentElement) -> CalculateSecondDerivativesContributions(this->mMatrix.M[thread],this->mVector.a[thread],rCurrentProcessInfo);

      (rCurrentElement) -> CalculateFirstDerivativesContributions(this->mMatrix.D[thread],this->mVector.v[thread],rCurrentProcessInfo);

      (rCurrentElement) -> EquationIdVector(rEquationId, rCurrentProcessInfo);

      AddDynamicsToLHS(rLHS_Contribution,this->mMatrix.D[thread],this->mMatrix.M[thread],rCurrentProcessInfo);

      AddDynamicForcesToRHS(rRHS_Contribution,this->mVector.v[thread],this->mVector.a[thread],rCurrentProcessInfo);

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

      int thread = OpenMPUtils::ThisThread();

      // Basic operations for the element considered
      (rCurrentElement) -> CalculateRightHandSide(rRHS_Contribution,rCurrentProcessInfo);

      (rCurrentElement) -> CalculateSecondDerivativesRHS(this->mVector.a[thread],rCurrentProcessInfo);

      (rCurrentElement) -> CalculateFirstDerivativesRHS(this->mVector.v[thread],rCurrentProcessInfo);

      (rCurrentElement) -> EquationIdVector(rEquationId,rCurrentProcessInfo);
      
      AddDynamicForcesToRHS(rRHS_Contribution,this->mVector.v[thread],this->mVector.a[thread],rCurrentProcessInfo);
      
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

      int thread = OpenMPUtils::ThisThread();

      // Basic operations for the condition considered
      (rCurrentCondition) -> CalculateLocalSystem(rLHS_Contribution,rRHS_Contribution,rCurrentProcessInfo);

      (rCurrentCondition) -> CalculateSecondDerivativesContributions(this->mMatrix.M[thread],this->mVector.a[thread],rCurrentProcessInfo);
	  
      (rCurrentCondition) -> CalculateFirstDerivativesContributions(this->mMatrix.D[thread],this->mVector.v[thread],rCurrentProcessInfo);

      (rCurrentCondition) -> EquationIdVector(rEquationId,rCurrentProcessInfo);

      AddDynamicsToLHS(rLHS_Contribution,this->mMatrix.D[thread],this->mMatrix.M[thread],rCurrentProcessInfo);

      AddDynamicForcesToRHS(rRHS_Contribution,this->mVector.v[thread],this->mVector.a[thread],rCurrentProcessInfo);

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

      int thread = OpenMPUtils::ThisThread();

      // Basic operations for the condition considered
      (rCurrentCondition) -> CalculateRightHandSide(rRHS_Contribution,rCurrentProcessInfo);

      (rCurrentCondition) -> CalculateSecondDerivativesRHS(this->mVector.a[thread],rCurrentProcessInfo);

      (rCurrentCondition) -> CalculateFirstDerivativesRHS(this->mVector.v[thread],rCurrentProcessInfo);

      (rCurrentCondition) -> EquationIdVector(rEquationId,rCurrentProcessInfo);

      AddDynamicForcesToRHS(rRHS_Contribution,this->mVector.v[thread],this->mVector.a[thread],rCurrentProcessInfo);
      
      KRATOS_CATCH( "" );
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
      ErrorCode  = DerivedType::Check(rModelPart);

      // Check that all required variables have been registered
      KRATOS_CHECK_VARIABLE_KEY(STEP_DISPLACEMENT);
      KRATOS_CHECK_VARIABLE_KEY(STEP_ROTATION);
      KRATOS_CHECK_VARIABLE_KEY(ROTATION);
      KRATOS_CHECK_VARIABLE_KEY(ANGULAR_VELOCITY);
      KRATOS_CHECK_VARIABLE_KEY(ANGULAR_ACCELERATION);

      // Check that variables are correctly allocated
      for(ModelPart::NodesContainerType::iterator it=rModelPart.NodesBegin(); it!=rModelPart.NodesEnd(); ++it)
        {
	  // Nodal data
	  KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(STEP_DISPLACEMENT,(*it));
	  KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(STEP_ROTATION,(*it));
	  KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ROTATION,(*it));
	  KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ANGULAR_VELOCITY,(*it));
	  KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ANGULAR_ACCELERATION,(*it));

	  // Nodal dofs
	  KRATOS_CHECK_DOF_IN_NODE(ROTATION_Z,(*it));
	  if( rModelPart.GetProcessInfo()[SPACE_DIMENSION] == 3 ){
	    KRATOS_CHECK_DOF_IN_NODE(ROTATION_X,(*it));
	    KRATOS_CHECK_DOF_IN_NODE(ROTATION_Y,(*it));
	  }
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
        buffer << "Displacement-Rotation NewmarkScheme";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Displacement-Rotation NewmarkScheme";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "Displacement-Rotation NewmarkScheme Data";     
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

    IntegrationPointerType    mpRotationIntegrationMethod;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    virtual void SetIntegrationMethod(ProcessInfo& rCurrentProcessInfo) override
    {
      
      this->mpIntegrationMethod = IntegrationPointerType( new NewmarkStepMethod<Variable<array_1d<double, 3> >, array_1d<double,3> > );

      // Set scheme variables
      this->mpIntegrationMethod->SetVariables(DISPLACEMENT,VELOCITY,ACCELERATION);

      this->mpIntegrationMethod->SetStepVariable(STEP_DISPLACEMENT);

      // Set scheme parameters
      this->mpIntegrationMethod->SetParameters(rCurrentProcessInfo);
       
      this->mpRotationIntegrationMethod = IntegrationPointerType( new NewmarkStepRotationMethod<Variable<array_1d<double, 3> >, array_1d<double,3> > );
            
      // Set rotation scheme variables
      this->mpRotationIntegrationMethod->SetVariables(ROTATION,ANGULAR_VELOCITY,ANGULAR_ACCELERATION);
      
      this->mpRotationIntegrationMethod->SetStepVariable(STEP_ROTATION);
      
      // Set scheme parameters
      this->mpRotationIntegrationMethod->SetParameters(rCurrentProcessInfo);

      // Modify ProcessInfo scheme parameters
      this->mpIntegrationMethod->SetProcessInfoParameters(rCurrentProcessInfo);
      rCurrentProcessInfo[COMPUTE_DYNAMIC_TANGENT] = true;      
    }

    virtual void IntegrationMethodUpdate(NodeType& rNode) override
    {
      this->mpIntegrationMethod->Update(rNode);
      this->mpRotationIntegrationMethod->Update(rNode);
    }

    virtual void IntegrationMethodPredict(NodeType& rNode) override
    {
      this->mpIntegrationMethod->Predict(rNode);
      this->mpRotationIntegrationMethod->Predict(rNode);
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
				  ProcessInfo& rCurrentProcessInfo) override
    {

      // Adding mass contribution to the dynamic stiffness
      if (rM.size1() != 0) // if M matrix declared
        {
	  noalias(rLHS_Contribution) += rM ;
        }

      // Adding  damping contribution
      if (rD.size1() != 0) // if D matrix declared
        {
	  noalias(rLHS_Contribution) += rD;
        }
    }


    /**
     * It adds the dynamic RHS contribution of the elements: b - M*a - D*v
     * @param rCurrentElement: The element to compute
     * @param rRHS_Contribution: The dynamic contribution for the RHS
     * @param rfv: The damping component vector
     * @param rfa: The mass component vector
     * @param rCurrentProcessInfo: The current process info instance
     */
    void AddDynamicForcesToRHS(LocalSystemVectorType& rRHS_Contribution,
			       LocalSystemVectorType& rfv ,
			       LocalSystemVectorType& rfa,
			       ProcessInfo& rCurrentProcessInfo)
    {

      // Adding inertia contribution
      if (rfa.size() != 0)
        {
	  noalias(rRHS_Contribution) -=  rfa;
        }

      // Adding damping contribution
      if (rfv.size() != 0)
        {
	  noalias(rRHS_Contribution) -=  rfv;
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
  }; // Class DisplacementRotationNewmarkScheme
  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{

  
  ///@}

  ///@} addtogroup block
  
}  // namespace Kratos.

#endif // KRATOS_DISPLACEMENT_ROTATION_NEWMARK_SCHEME_H_INCLUDED defined
