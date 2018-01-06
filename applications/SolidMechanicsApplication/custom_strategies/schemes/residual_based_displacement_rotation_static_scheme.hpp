//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:            November 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_RESIDUAL_BASED_DISPLACEMENT_ROTATION_STATIC_SCHEME )
#define  KRATOS_RESIDUAL_BASED_DISPLACEMENT_ROTATION_STATIC_SCHEME

// System includes

// External includes

// Project includes
#include "custom_strategies/schemes/residual_based_displacement_static_scheme.hpp"

#include "custom_strategies/time_integration_methods/static_step_rotation_method.hpp"

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
  class ResidualBasedDisplacementRotationStaticScheme: public ResidualBasedDisplacementStaticScheme<TSparseSpace,TDenseSpace>
  {   
  public:
    
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedDisplacementRotationStaticScheme );

    typedef Scheme<TSparseSpace,TDenseSpace>                      BaseType;
    
    typedef typename BaseType::Pointer                     BaseTypePointer;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef ResidualBasedDisplacementStaticScheme<TSparseSpace,TDenseSpace>   DerivedType;

    typedef typename DerivedType::IntegrationTypePointer           IntegrationTypePointer;
    
    typedef typename DerivedType::NodeType                                       NodeType;
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default Constructor.
    ResidualBasedDisplacementRotationStaticScheme()
      :DerivedType()
    {
    }

    /// Copy Constructor.
    ResidualBasedDisplacementRotationStaticScheme(ResidualBasedDisplacementRotationStaticScheme& rOther)
      :DerivedType(rOther)
      ,mpRotationIntegrationMethod(rOther.mpRotationIntegrationMethod)
    {
    }

    /// Clone.
    BaseTypePointer Clone() override
    {
      return BaseTypePointer( new ResidualBasedDisplacementRotationStaticScheme(*this) );
    }

    /// Destructor.
    ~ResidualBasedDisplacementRotationStaticScheme() override {}

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
     
      (rCurrentElement) -> CalculateLocalSystem(LHS_Contribution,RHS_Contribution, rCurrentProcessInfo);
      
      (rCurrentElement) -> EquationIdVector(EquationId, rCurrentProcessInfo);

      KRATOS_CATCH( "" );
    }

    /**
     * This function is designed to calculate just the RHS contribution
     * @param rCurrentElement: The element to compute
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

      // Basic operations for the element considered
      (rCurrentElement) -> CalculateRightHandSide(RHS_Contribution,rCurrentProcessInfo);

      (rCurrentElement) -> EquationIdVector(EquationId,rCurrentProcessInfo);

      KRATOS_CATCH( "" );
    }

    /**
     * This function is designed to calculate just the LHS contribution
     * @param rCurrentElement: The element to compute
     * @param LHS_Contribution: The LHS matrix contribution
     * @param EquationId: The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo: The current process info instance
     */

    void Calculate_LHS_Contribution(Element::Pointer rCurrentElement,
				    LocalSystemMatrixType& LHS_Contribution,
				    Element::EquationIdVectorType& EquationId,
				    ProcessInfo& rCurrentProcessInfo) override
    {

      KRATOS_TRY;

      // Basic operations for the element considered
      (rCurrentElement) -> CalculateLeftHandSide(LHS_Contribution,rCurrentProcessInfo);

      (rCurrentElement) -> EquationIdVector(EquationId,rCurrentProcessInfo);

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

      // Basic operations for the condition considered
      (rCurrentCondition) -> CalculateLocalSystem(LHS_Contribution,RHS_Contribution,rCurrentProcessInfo);

      (rCurrentCondition) -> EquationIdVector(EquationId,rCurrentProcessInfo);

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

      // Basic operations for the condition considered
      (rCurrentCondition) -> CalculateRightHandSide(RHS_Contribution,rCurrentProcessInfo);

      (rCurrentCondition) -> EquationIdVector(EquationId,rCurrentProcessInfo);
      
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

      // Check that variables are correctly allocated
      for(ModelPart::NodesContainerType::iterator it=rModelPart.NodesBegin(); it!=rModelPart.NodesEnd(); it++)
        {
	  // Nodal data
	  KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(STEP_DISPLACEMENT,(*it));
	  KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(STEP_ROTATION,(*it));
	  KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ROTATION,(*it));

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
        buffer << "Displacement-Rotation StaticScheme";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Displacement-Rotation StaticScheme";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "Displacement-Rotation StaticScheme Data";     
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
      this->mpIntegrationMethod = IntegrationTypePointer( new StaticStepMethod<Variable<array_1d<double, 3> >, array_1d<double,3> > );

      // Set scheme variables
      this->mpIntegrationMethod->SetVariable(DISPLACEMENT);

      this->mpIntegrationMethod->SetStepVariable(STEP_DISPLACEMENT);

      // Set scheme parameters
      this->mpIntegrationMethod->SetParameters(rCurrentProcessInfo);
      
      this->mpRotationIntegrationMethod = IntegrationTypePointer( new StaticStepRotationMethod<Variable<array_1d<double, 3> >, array_1d<double,3> > );
     
      // Set rotation scheme variables
      this->mpRotationIntegrationMethod->SetVariable(ROTATION);
      
      this->mpRotationIntegrationMethod->SetStepVariable(STEP_ROTATION);

      // Set scheme parameters
      this->mpRotationIntegrationMethod->SetParameters(rCurrentProcessInfo);
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
  }; // Class ResidualBasedDisplacementRotationStaticScheme
  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{

  
  ///@}

  ///@} addtogroup block
  
}  // namespace Kratos.

#endif // KRATOS_RESIDUAL_BASED_DISPLACEMENT_ROTATION_STATIC_SCHEME defined
