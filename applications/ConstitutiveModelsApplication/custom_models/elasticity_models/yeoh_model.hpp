ç//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_YEOH_MODEL_H_INCLUDED)
#define  KRATOS_YEOH_MODEL_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_models/elasticity_models/mooney_rivlin_model.hpp"

namespace Kratos
{
  ///@addtogroup ConstitutiveModelsApplication
  ///@{

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

  /// Short class definition.
  /** Detail class definition.
   */
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) YeohModel : public HyperElasticModel
  {
  public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of YeohModel
    KRATOS_CLASS_POINTER_DEFINITION( YeohModel );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    YeohModel() : HyperElasticModel() {}

    /// Copy constructor.
    YeohModel(YeohModel const& rOther) : HyperElasticModel(rOther) {}

    /// Assignment operator.
    YeohModel& operator=(YeohModel const& rOther)
    {
      HyperElasticModel::operator=(rOther);
      return *this;
    }

    /// Clone.
    virtual ConstitutiveModel::Pointer Clone() const override
    {
      return Kratos::make_shared<YeohModel>(*this);
    }

    /// Destructor.
    virtual ~YeohModel() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    virtual void CalculateStrainEnergy(ModelDataType& rValues, double& rDensityFunction)
    {
      KRATOS_TRY

      HyperElasticDataType Variables;
      this->CalculateStrainData(rValues, Variables);

      rDensityFunction += rVariables.GetMaterialParameters()[0] * ( Variables.Strain.Invariants.I1 - 3.0) + rVariables.GetMaterialParameters()[1] * ( Variables.Strain.Invariants.I2 - 3.0) + rVariables.GetMaterialParameters()[2] * ( Variables.Strain.Invariants.I3 - 3.0);

      KRATOS_CATCH(" ")
    }



    virtual int Check(const Properties& rProperties,
		      const ProcessInfo& rCurrentProcessInfo)
    {
      KRATOS_TRY

      HyperElasticModel::Check(rProperties,rCurrentProcessInfo);

      if( C10.Key() == 0 || rProperties[C10] <= 0.00 )
	KRATOS_ERROR << "C10 has an invalid key or value" << std::endl;

      if( C20.Key() == 0 || rProperties[C20] <= 0.00 )
	KRATOS_ERROR << "C20 has an invalid key or value" << std::endl;

      if( C30.Key() == 0 || rProperties[C30] <= 0.00 )
	KRATOS_ERROR << "C30 has an invalid key or value" << std::endl;

      if( BULK_MODULUS.Key() == 0 || rProperties[BULK_MODULUS] <= 0.00 )
	KRATOS_ERROR << "BULK_MODULUS has an invalid key or value" << std::endl;

      KRATOS_CATCH(" ")
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
        buffer << "YeohModel";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "YeohModel";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "YeohModel Data";
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

    //************// dW

    virtual double& GetFunction1stI1Derivative(HyperElasticDataType& rVariables, double& rDerivative) //dW/dI1
    {
      KRATOS_TRY

      rDerivative  = rVariables.GetMaterialParameters()[0];
      rDerivative += 2.0 * rVariables.GetMaterialParameters()[1] * (rVariables.Strain.Invariants.I1-3);
      rDerivative += 3.0 * rVariables.GetMaterialParameters()[2] * (rVariables.Strain.Invariants.I1-3) * (rVariables.Strain.Invariants.I1-3);

      return rDerivative;

      KRATOS_CATCH(" ")
    }

    virtual double& GetFunction1stI2Derivative(HyperElasticDataType& rVariables, double& rDerivative) //dW/dI2
    {
      KRATOS_TRY

      rDerivative = 0;

      return rDerivative;

      KRATOS_CATCH(" ")
    }

    virtual double& GetFunction1stI3Derivative(HyperElasticDataType& rVariables, double& rDerivative) //dW/dI3
    {
      KRATOS_TRY

      rDerivative  = 0;

      return rDerivative;

      KRATOS_CATCH(" ")
    }


    virtual double& GetFunction2ndI1Derivative(HyperElasticDataType& rVariables, double& rDerivative) //ddW/dI1dI1
    {
      KRATOS_TRY

      rDerivative = 0.0;

      return rDerivative;

      KRATOS_CATCH(" ")
    }

    virtual double& GetFunction2ndI2Derivative(HyperElasticDataType& rVariables, double& rDerivative) //ddW/dI2dI2
    {
      KRATOS_TRY

      rDerivative = 0.0;

      return rDerivative;

      KRATOS_CATCH(" ")
    }

    virtual double& GetFunction2ndI3Derivative(HyperElasticDataType& rVariables, double& rDerivative) //ddW/dI3dI3
    {
      KRATOS_TRY

      rDerivative  = 0;

      return rDerivative;

      KRATOS_CATCH(" ")
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
    friend class Serializer;


    virtual void save(Serializer& rSerializer) const override
    {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, HyperElasticModel )
    }

    virtual void load(Serializer& rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, HyperElasticModel )
    }

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

  }; // Class YeohModel

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_YEOH_MODEL_H_INCLUDED  defined
