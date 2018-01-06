//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:            November 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_RESIDUAL_BASED_DISPLACEMENT_ROTATION_NEWMARK_SCHEME )
#define  KRATOS_RESIDUAL_BASED_DISPLACEMENT_ROTATION_NEWMARK_SCHEME

// System includes

// External includes

// Project includes
#include "custom_strategies/schemes/residual_based_displacement_newmark_scheme.hpp"

#include "custom_strategies/time_integration_methods/newmark_step_rotation_method.hpp"

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
  class ResidualBasedDisplacementRotationNewmarkScheme: public ResidualBasedDisplacementNewmarkScheme<TSparseSpace,TDenseSpace>
  {   
  public:
    
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedDisplacementRotationNewmarkScheme );

    typedef Scheme<TSparseSpace,TDenseSpace>                      BaseType;
    
    typedef typename BaseType::Pointer                     BaseTypePointer;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef ResidualBasedDisplacementNewmarkScheme<TSparseSpace,TDenseSpace>  DerivedType;

    typedef typename DerivedType::IntegrationTypePointer           IntegrationTypePointer;
    
    typedef typename DerivedType::NodeType                                       NodeType;
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default Constructor.
    ResidualBasedDisplacementRotationNewmarkScheme()
      :DerivedType()
    {
    }

    /// Copy Constructor.
    ResidualBasedDisplacementRotationNewmarkScheme(ResidualBasedDisplacementRotationNewmarkScheme& rOther)
      :DerivedType(rOther)
      ,mpRotationIntegrationMethod(rOther.mpRotationIntegrationMethod)
    {
    }

    /// Clone.
    BaseTypePointer Clone() override
    {
      return BaseTypePointer( new ResidualBasedDisplacementRotationNewmarkScheme(*this) );
    }

    /// Destructor.
    ~ResidualBasedDisplacementRotationNewmarkScheme() override {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{


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
      
      (rCurrentElement) -> CalculateLocalSystem(LHS_Contribution,RHS_Contribution, rCurrentProcessInfo);
      
      (rCurrentElement) -> CalculateSecondDerivativesContributions(this->mMatrix.M[thread],this->mVector.a[thread],rCurrentProcessInfo);

      (rCurrentElement) -> CalculateFirstDerivativesContributions(this->mMatrix.D[thread],this->mVector.v[thread],rCurrentProcessInfo);

      (rCurrentElement) -> EquationIdVector(EquationId, rCurrentProcessInfo);

      AddDynamicsToLHS(LHS_Contribution,this->mMatrix.D[thread],this->mMatrix.M[thread],rCurrentProcessInfo);

      AddDynamicForcesToRHS(RHS_Contribution,this->mVector.v[thread],this->mVector.a[thread],rCurrentProcessInfo);

      KRATOS_CATCH( "" );
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

      // Basic operations for the element considered
      (rCurrentElement) -> CalculateRightHandSide(RHS_Contribution,rCurrentProcessInfo);

      (rCurrentElement) -> CalculateSecondDerivativesRHS(this->mVector.a[thread],rCurrentProcessInfo);

      (rCurrentElement) -> CalculateFirstDerivativesRHS(this->mVector.v[thread],rCurrentProcessInfo);

      (rCurrentElement) -> EquationIdVector(EquationId,rCurrentProcessInfo);
      
      AddDynamicForcesToRHS(RHS_Contribution,this->mVector.v[thread],this->mVector.a[thread],rCurrentProcessInfo);
      
      KRATOS_CATCH( "" );
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

      // Basic operations for the condition considered
      (rCurrentCondition) -> CalculateLocalSystem(LHS_Contribution,RHS_Contribution,rCurrentProcessInfo);

      (rCurrentCondition) -> CalculateSecondDerivativesContributions(this->mMatrix.M[thread],this->mVector.a[thread],rCurrentProcessInfo);
	  
      (rCurrentCondition) -> CalculateFirstDerivativesContributions(this->mMatrix.D[thread],this->mVector.v[thread],rCurrentProcessInfo);

      (rCurrentCondition) -> EquationIdVector(EquationId,rCurrentProcessInfo);

      AddDynamicsToLHS(LHS_Contribution,this->mMatrix.D[thread],this->mMatrix.M[thread],rCurrentProcessInfo);

      AddDynamicForcesToRHS(RHS_Contribution,this->mVector.v[thread],this->mVector.a[thread],rCurrentProcessInfo);

      KRATOS_CATCH( "" );
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

      // Basic operations for the condition considered
      (rCurrentCondition) -> CalculateRightHandSide(RHS_Contribution,rCurrentProcessInfo);

      (rCurrentCondition) -> CalculateSecondDerivativesRHS(this->mVector.a[thread],rCurrentProcessInfo);

      (rCurrentCondition) -> CalculateFirstDerivativesRHS(this->mVector.v[thread],rCurrentProcessInfo);

      (rCurrentCondition) -> EquationIdVector(EquationId,rCurrentProcessInfo);

      AddDynamicForcesToRHS(RHS_Contribution,this->mVector.v[thread],this->mVector.a[thread],rCurrentProcessInfo);
      
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
      for(ModelPart::NodesContainerType::iterator it=rModelPart.NodesBegin(); it!=rModelPart.NodesEnd(); it++)
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

    IntegrationTypePointer    mpRotationIntegrationMethod;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    virtual void SetIntegrationMethod(ProcessInfo& rCurrentProcessInfo) override
    {
      
      this->mpIntegrationMethod = IntegrationTypePointer( new NewmarkStepMethod<Variable<array_1d<double, 3> >, array_1d<double,3> > );

      // Set scheme variables
      this->mpIntegrationMethod->SetVariables(DISPLACEMENT,VELOCITY,ACCELERATION);

      this->mpIntegrationMethod->SetStepVariable(STEP_DISPLACEMENT);

      // Set scheme parameters
      this->mpIntegrationMethod->SetParameters(rCurrentProcessInfo);
       
      this->mpRotationIntegrationMethod = IntegrationTypePointer( new NewmarkStepRotationMethod<Variable<array_1d<double, 3> >, array_1d<double,3> > );
            
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
     * @param LHS_Contribution: The dynamic contribution for the LHS
     * @param D: The damping matrix
     * @param M: The mass matrix
     * @param rCurrentProcessInfo: The current process info instance
     */

    virtual void AddDynamicsToLHS(LocalSystemMatrixType& LHS_Contribution,
				  LocalSystemMatrixType& D,
				  LocalSystemMatrixType& M,
				  ProcessInfo& rCurrentProcessInfo) override
    {

      // Adding mass contribution to the dynamic stiffness
      if (M.size1() != 0) // if M matrix declared
        {
	  noalias(LHS_Contribution) += M ;
        }

      // Adding  damping contribution
      if (D.size1() != 0) // if D matrix declared
        {
	  noalias(LHS_Contribution) += D;
        }
    }


    /**
     * It adds the dynamic RHS contribution of the elements: b - M*a - D*v
     * @param rCurrentElement: The element to compute
     * @param RHS_Contribution: The dynamic contribution for the RHS
     * @param fv: The damping component vector
     * @param fa: The mass component vector
     * @param CurrentProcessInfo: The current process info instance
     */
    void AddDynamicForcesToRHS(LocalSystemVectorType& RHS_Contribution,
			       LocalSystemVectorType& fv ,
			       LocalSystemVectorType& fa,
			       ProcessInfo& CurrentProcessInfo)
    {

      // Adding inertia contribution
      if (fa.size() != 0)
        {
	  noalias(RHS_Contribution) -=  fa;
        }

      // Adding damping contribution
      if (fv.size() != 0)
        {
	  noalias(RHS_Contribution) -=  fv;
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
  }; // Class ResidualBasedDisplacementRotationNewmarkScheme
  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{

  
  ///@}

  ///@} addtogroup block
  
}  // namespace Kratos.

#endif // KRATOS_RESIDUAL_BASED_DISPLACEMENT_ROTATION_NEWMARK_SCHEME defined
