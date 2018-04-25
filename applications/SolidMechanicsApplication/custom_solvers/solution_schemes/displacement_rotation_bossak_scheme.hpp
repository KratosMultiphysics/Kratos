//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:            November 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_DISPLACEMENT_ROTATION_BOSSAK_SCHEME_H_INCLUDED)
#define  KRATOS_DISPLACEMENT_ROTATION_BOSSAK_SCHEME_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_solvers/solution_schemes/displacement_rotation_newmark_scheme.hpp"

#include "custom_solvers/time_integration_methods/bossak_step_method.hpp"
#include "custom_solvers/time_integration_methods/bossak_step_rotation_method.hpp"

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
  class DisplacementRotationBossakScheme: public DisplacementRotationNewmarkScheme<TSparseSpace,TDenseSpace>
  {
  public:

    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION( DisplacementRotationBossakScheme );

    typedef SolutionScheme<TSparseSpace,TDenseSpace>                             BaseType;
    typedef typename BaseType::SolutionSchemePointerType                  BasePointerType;

    typedef typename BaseType::LocalSystemVectorType                LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType                LocalSystemMatrixType;

    typedef DisplacementRotationNewmarkScheme<TSparseSpace,TDenseSpace>       DerivedType;

    typedef typename DerivedType::IntegrationPointerType           IntegrationPointerType;

    typedef typename DerivedType::NodeType                                       NodeType;
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default Constructor.
    DisplacementRotationBossakScheme()
      :DerivedType()
    {
    }

    /// Constructor.
    DisplacementRotationBossakScheme(Flags& rOptions)
      :DerivedType(rOptions)
    {
    }

    /// Copy Constructor.
    DisplacementRotationBossakScheme(DisplacementRotationBossakScheme& rOther)
      :DerivedType(rOther)
    {
    }

    /// Clone.
    BasePointerType Clone() override
    {
      return BasePointerType( new DisplacementRotationBossakScheme(*this) );
    }

    /// Destructor.
    ~DisplacementRotationBossakScheme() override {}

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

    void SetIntegrationMethod(ProcessInfo& rCurrentProcessInfo) override
    {
      if ( this->mTimeIntegrationMethods.size() == 0 ) {
        this->mTimeIntegrationMethods.push_back(Kratos::make_shared< BossakStepMethod<Variable<array_1d<double, 3> >, array_1d<double,3> > >(DISPLACEMENT,VELOCITY,ACCELERATION));

        // Set scheme variables
        this->mTimeIntegrationMethods.front()->SetStepVariable(STEP_DISPLACEMENT);

        // Set scheme parameters
        this->mTimeIntegrationMethods.front()->SetParameters(rCurrentProcessInfo);


        this->mTimeIntegrationMethods.push_back(Kratos::make_shared< BossakStepRotationMethod<Variable<array_1d<double, 3> >, array_1d<double,3> > >(ROTATION,ANGULAR_VELOCITY,ANGULAR_ACCELERATION));


        this->mTimeIntegrationMethods.back()->SetStepVariable(STEP_ROTATION);

        // Set scheme parameters
        this->mTimeIntegrationMethods.back()->SetParameters(rCurrentProcessInfo);

        // Set parameters to process info
        this->mTimeIntegrationMethods.back()->SetProcessInfoParameters(rCurrentProcessInfo);
        
        // Modify ProcessInfo scheme parameters
        rCurrentProcessInfo[COMPUTE_DYNAMIC_TANGENT] = true;
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
  }; // Class DisplacementRotationBossakScheme
  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_DISPLACEMENT_ROTATION_BOSSAK_SCHEME_H_INCLUDED defined
