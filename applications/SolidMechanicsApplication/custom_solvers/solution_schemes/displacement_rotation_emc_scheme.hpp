//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:            November 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_DISPLACEMENT_ROTATION_EMC_SCHEME_H_INCLUDED)
#define  KRATOS_DISPLACEMENT_ROTATION_EMC_SCHEME_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_solvers/solution_schemes/displacement_rotation_simo_scheme.hpp"

#include "custom_solvers/time_integration_methods/emc_step_rotation_method.hpp"

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

  /** @brief Bossak integration scheme (for dynamic problems)
   */
  template<class TSparseSpace,  class TDenseSpace >
  class DisplacementRotationEmcScheme: public DisplacementRotationSimoScheme<TSparseSpace,TDenseSpace>
  {   
  public:
    
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION( DisplacementRotationEmcScheme );

    typedef SolutionScheme<TSparseSpace,TDenseSpace>                                 BaseType;
    typedef typename BaseType::SolutionSchemePointerType                      BasePointerType;

    typedef typename BaseType::LocalSystemVectorType                     LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType                     LocalSystemMatrixType;

    typedef DisplacementRotationSimoScheme<TSparseSpace,TDenseSpace>               DerivedType;

    typedef typename DerivedType::IntegrationPointerType                IntegrationPointerType;
      
    typedef typename DerivedType::NodeType                                            NodeType;
    
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default Constructor.
    DisplacementRotationEmcScheme()
      :DerivedType()
    {
    }

    /// Constructor.
    DisplacementRotationEmcScheme(Flags& rOptions)
      :DerivedType(rOptions)
    {
    }
    
    /// Copy Constructor.
    DisplacementRotationEmcScheme(DisplacementRotationEmcScheme& rOther)
      :DerivedType(rOther)
    {
    }

    /// Clone.
    BasePointerType Clone() override
    {
      return BasePointerType( new DisplacementRotationEmcScheme(*this) );
    }

    /// Destructor.
    ~DisplacementRotationEmcScheme() override {}

    ///@}
    ///@name Operators
    ///@{

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
    virtual std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "Displacement-Rotation EMCScheme";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Displacement-Rotation EMCScheme";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "Displacement-Rotation EMCScheme Data";     
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

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{
    
    virtual void SetIntegrationMethod(ProcessInfo& rCurrentProcessInfo) override
    {      
      this->mpIntegrationMethod = IntegrationPointerType( new EmcStepMethod<Variable<array_1d<double, 3> >, array_1d<double,3> > );

      // Set scheme variables
      this->mpIntegrationMethod->SetVariables(DISPLACEMENT,VELOCITY,ACCELERATION);

      this->mpIntegrationMethod->SetStepVariable(STEP_DISPLACEMENT);
      
      // Set scheme parameters
      this->mpIntegrationMethod->SetParameters(rCurrentProcessInfo);
      
      this->mpRotationIntegrationMethod = IntegrationPointerType( new EmcStepRotationMethod<Variable<array_1d<double, 3> >, array_1d<double,3> > );

      // Set rotation scheme variables
      this->mpRotationIntegrationMethod->SetVariables(ROTATION,ANGULAR_VELOCITY,ANGULAR_ACCELERATION);
      
      this->mpRotationIntegrationMethod->SetStepVariable(STEP_ROTATION);
      
      // Set scheme parameters
      this->mpRotationIntegrationMethod->SetParameters(rCurrentProcessInfo);

      // Modify ProcessInfo scheme parameters
      this->mpIntegrationMethod->SetProcessInfoParameters(rCurrentProcessInfo);
      
      rCurrentProcessInfo[COMPUTE_DYNAMIC_TANGENT] = true;
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
  }; // Class DisplacementRotationEmcScheme
  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{

  
  ///@}

  ///@} addtogroup block
  
}  // namespace Kratos.

#endif // KRATOS_DISPLACEMENT_ROTATION_EMC_SCHEME_H_INCLUDED defined
