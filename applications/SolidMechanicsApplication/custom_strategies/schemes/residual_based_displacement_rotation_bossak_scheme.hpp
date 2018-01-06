//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:            November 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_RESIDUAL_BASED_DISPLACEMENT_ROTATION_BOSSAK_SCHEME )
#define  KRATOS_RESIDUAL_BASED_DISPLACEMENT_ROTATION_BOSSAK_SCHEME

// System includes

// External includes

// Project includes
#include "custom_strategies/schemes/residual_based_displacement_rotation_newmark_scheme.hpp"

#include "custom_strategies/time_integration_methods/bossak_step_method.hpp"
#include "custom_strategies/time_integration_methods/bossak_step_rotation_method.hpp"

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
  class ResidualBasedDisplacementRotationBossakScheme: public ResidualBasedDisplacementRotationNewmarkScheme<TSparseSpace,TDenseSpace>
  {   
  public:
    
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedDisplacementRotationBossakScheme );

    typedef Scheme<TSparseSpace,TDenseSpace>                      BaseType;
    
    typedef typename BaseType::Pointer                     BaseTypePointer;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef ResidualBasedDisplacementRotationNewmarkScheme<TSparseSpace,TDenseSpace>  DerivedType;

    typedef typename DerivedType::IntegrationTypePointer           IntegrationTypePointer;

    typedef typename DerivedType::NodeType                                       NodeType;
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default Constructor.
    ResidualBasedDisplacementRotationBossakScheme()
      :DerivedType()
    {
    }

    /// Copy Constructor.
    ResidualBasedDisplacementRotationBossakScheme(ResidualBasedDisplacementRotationBossakScheme& rOther)
      :DerivedType(rOther)
    {
    }

    /// Clone.
    BaseTypePointer Clone() override
    {
      return BaseTypePointer( new ResidualBasedDisplacementRotationBossakScheme(*this) );
    }

    /// Destructor.
    ~ResidualBasedDisplacementRotationBossakScheme() override {}

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

	DerivedType::Initialize(rModelPart);  

	const unsigned int NumThreads = OpenMPUtils::GetNumThreads();
	
	this->mVector.ap.resize(NumThreads);

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
        buffer << "Displacement-Rotation BossakScheme";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Displacement-Rotation BossakScheme";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "Displacement-Rotation BossakScheme Data";     
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
      this->mpIntegrationMethod = IntegrationTypePointer( new BossakStepMethod<Variable<array_1d<double, 3> >, array_1d<double,3> > );

      // Set scheme variables
      this->mpIntegrationMethod->SetVariables(DISPLACEMENT,VELOCITY,ACCELERATION);

      this->mpIntegrationMethod->SetStepVariable(STEP_DISPLACEMENT);
      
      // Set scheme parameters
      this->mpIntegrationMethod->SetParameters(rCurrentProcessInfo);

      this->mpRotationIntegrationMethod = IntegrationTypePointer( new BossakStepRotationMethod<Variable<array_1d<double, 3> >, array_1d<double,3> > );

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
  }; // Class ResidualBasedDisplacementRotationBossakScheme
  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{

  
  ///@}

  ///@} addtogroup block
  
}  // namespace Kratos.

#endif // KRATOS_RESIDUAL_BASED_DISPLACEMENT_ROTATION_BOSSAK_SCHEME defined
