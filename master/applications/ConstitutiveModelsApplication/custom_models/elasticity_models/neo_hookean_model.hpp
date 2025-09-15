//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_NEO_HOOKEAN_MODEL_H_INCLUDED)
#define  KRATOS_NEO_HOOKEAN_MODEL_H_INCLUDED

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
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) NeoHookeanModel : public MooneyRivlinModel
  {
  public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of NeoHookeanModel
    KRATOS_CLASS_POINTER_DEFINITION(NeoHookeanModel);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    NeoHookeanModel() : MooneyRivlinModel() {}

    /// Copy constructor.
    NeoHookeanModel(NeoHookeanModel const& rOther) : MooneyRivlinModel(rOther) {}

    /// Assignment operator.
    NeoHookeanModel& operator=(NeoHookeanModel const& rOther)
    {
      MooneyRivlinModel::operator=(rOther);
      return *this;
    }

    /// Clone.
    ConstitutiveModel::Pointer Clone() const override
    {
      return Kratos::make_shared<NeoHookeanModel>(*this);
    }

    /// Destructor.
    ~NeoHookeanModel() override {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void CalculateStrainEnergy(ModelDataType& rValues, double& rDensityFunction) override
    {
      KRATOS_TRY

      HyperElasticDataType Variables;
      this->CalculateStrainData(rValues, Variables);

      this->CalculateAndAddVolumetricStrainEnergy( Variables, rDensityFunction );

      rDensityFunction += rValues.MaterialParameters.GetModelParameters()[0] * ( Variables.Strain.Invariants.I1 - 3.0);

      KRATOS_CATCH(" ")
    }



    int Check(const Properties& rProperties, const ProcessInfo& rCurrentProcessInfo) override
    {
      KRATOS_TRY

      MooneyRivlinModel::Check(rProperties,rCurrentProcessInfo);

      if( C10.Key() == 0 || rProperties[C10] <= 0.00 )
	KRATOS_ERROR << "C10 has an invalid key or value" << std::endl;
      if( BULK_MODULUS.Key() == 0 || rProperties[BULK_MODULUS] <= 0.00 )
	KRATOS_ERROR << "BULK_MODULUS has an invalid key or value" << std::endl;

      return 0;

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
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "NeoHookeanModel";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "NeoHookeanModel";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "NeoHookeanModel Data";
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


    // specialized methods:

    virtual void CalculateVolumetricFactor(HyperElasticDataType& rVariables, double& rFactor)
    {
      KRATOS_TRY

      rFactor = 0.5 * (rVariables.Strain.Invariants.I3-1.0);

      KRATOS_CATCH(" ")
    }

    virtual void CalculateConstitutiveMatrixFactor(HyperElasticDataType& rVariables, double& rFactor)
    {
      KRATOS_TRY

      rFactor = rVariables.Strain.Invariants.I3;

      KRATOS_CATCH(" ")
    }


    //**** Particular calculation for the basic neo-hookean models ****//

    // virtual void CalculateAndAddStressTensor(HyperElasticDataType& rVariables, MatrixType& rStressMatrix) override
    // {
    //   KRATOS_TRY

    //   const ModelDataType&  rModelData        = rVariables.GetModelData();
    //   const StressMeasureType& rStressMeasure = rModelData.GetStressMeasure();

    //   MatrixType StressMatrix;

    //   const MaterialDataType& rMaterial = rVariables.GetMaterialParameters();

    //   double Factor = 0;
    //   this->CalculateVolumetricFactor(rVariables,Factor);

    //   if( rStressMeasure == ConstitutiveModelData::StressMeasureType::StressMeasure_PK2 ){ //mStrainMatrix = RightCauchyGreen (C)

    // 	StressMatrix  = rMaterial.GetLameLambda() * Factor * rVariables.Strain.InverseMatrix;
    //     StressMatrix += rMaterial.GetLameMu() * ( msIdentityMatrix - rVariables.Strain.InverseMatrix);

    // 	rStressMatrix += StressMatrix;

    //   }
    //   else if( rStressMeasure == ConstitutiveModelData::StressMeasureType::StressMeasure_Kirchhoff ){ //mStrainMatrix = LeftCauchyGreen (b)

    // 	StressMatrix  = rMaterial.GetLameLambda() * Factor * msIdentityMatrix;
    //     StressMatrix += rMaterial.GetLameMu() * ( rVariables.Strain.Matrix - msIdentityMatrix );

    // 	rStressMatrix += StressMatrix;
    //   }

    //   rVariables.State().Set(ConstitutiveModelData::STRESS_COMPUTED);

    //   KRATOS_CATCH(" ")
    // }


    // virtual double& AddConstitutiveComponent(HyperElasticDataType& rVariables, double &rCabcd,
    // 					     const unsigned int& a, const unsigned int& b,
    // 					     const unsigned int& c, const unsigned int& d) override
    // {
    //   KRATOS_TRY


    //   double Cabcd = 0;

    //   const MaterialDataType& rMaterial = rVariables.GetMaterialParameters();

    //   const ModelDataType&  rModelData        = rVariables.GetModelData();
    //   const StressMeasureType& rStressMeasure = rModelData.GetStressMeasure();

    //   double FactorA = 0;
    //   this->CalculateConstitutiveMatrixFactor(rVariables,FactorA);

    //   double FactorB = 0;
    //   this->CalculateVolumetricFactor(rVariables,FactorB);

    //   if( rStressMeasure == ConstitutiveModelData::StressMeasureType::StressMeasure_PK2 ){ //mStrainMatrix = RightCauchyGreen (C)
    // 	Cabcd  = FactorA * rMaterial.GetLameLambda() * (rVariables.Strain.InverseMatrix(a,b)*rVariables.Strain.InverseMatrix(c,d));
    // 	Cabcd += (rMaterial.GetLameMu() - rMaterial.GetLameLambda() * FactorB) * (rVariables.Strain.InverseMatrix(a,c)*rVariables.Strain.InverseMatrix(b,d)+rVariables.Strain.InverseMatrix(a,d)*rVariables.Strain.InverseMatrix(b,c));
    //   }
    //   else if( rStressMeasure == ConstitutiveModelData::StressMeasureType::StressMeasure_Kirchhoff ){ //mStrainMatrix = LeftCauchyGreen (b)
    // 	Cabcd  = FactorA * rMaterial.GetLameLambda() * (msIdentityMatrix(a,b)*msIdentityMatrix(c,d));
    // 	Cabcd += (rMaterial.GetLameMu() - rMaterial.GetLameLambda() * FactorB) * (msIdentityMatrix(a,c)*msIdentityMatrix(b,d)+msIdentityMatrix(a,d)*msIdentityMatrix(b,c));
    //   }

    //   rCabcd += Cabcd;

    //   rVariables.State().Set(ConstitutiveModelData::CONSTITUTIVE_MATRIX_COMPUTED);

    //   return rCabcd;

    //   KRATOS_CATCH(" ")
    // }

    //**** Particular calculation for the basic neo-hookean models ****//

    //************// W

    //option: g(J) = (lambda/4)*(J²-1) - (lamda/2)*lnJ - (mu)*lnJ  (default implemented)

    void CalculateAndAddVolumetricStrainEnergy(HyperElasticDataType& rVariables, double& rVolumetricDensityFunction) override
    {
      KRATOS_TRY

      const MaterialDataType& rMaterial = rVariables.GetMaterialParameters();

      //g(J) = (lambda/4)*(J²-1) - (lamda/2)*lnJ - (mu)*lnJ
      rVolumetricDensityFunction  = rMaterial.GetLameLambda() * 0.25 * ( rVariables.Strain.Invariants.I3 - 1.0);
      rVolumetricDensityFunction -= rMaterial.GetLameLambda() * 0.5 * std::log( rVariables.Strain.Invariants.J );
      rVolumetricDensityFunction -= rMaterial.GetLameMu() * std::log( rVariables.Strain.Invariants.J );

      KRATOS_CATCH(" ")
    }

    //************// dW

    double& GetFunction1stI1Derivative(HyperElasticDataType& rVariables, double& rDerivative) override //dW/dI1
    {
      KRATOS_TRY

      const MaterialDataType& rMaterial = rVariables.GetMaterialParameters();
      rDerivative = rMaterial.GetModelParameters()[0];

      return rDerivative;

      KRATOS_CATCH(" ")
    }

    double& GetFunction1stI2Derivative(HyperElasticDataType& rVariables, double& rDerivative) override //dW/dI2
    {
      KRATOS_TRY

      rDerivative =  0.0;

      return rDerivative;

      KRATOS_CATCH(" ")
    }

    double& GetFunction1stI3Derivative(HyperElasticDataType& rVariables, double& rDerivative) override //dW/dI3
    {
      KRATOS_TRY

      const MaterialDataType&  rMaterial = rVariables.GetMaterialParameters();

      //derivative of "g(J) = (lambda/4)*(J²-1) - (lamda/2)*lnJ - (mu)*lnJ"
      //dg(J)/dI3 = (lambda/4) - (lambda/2)*(1/J²)/2 - (mu)*(1/J²)/2
      rDerivative  = rMaterial.GetLameLambda();
      rDerivative += 2.0 * rMaterial.GetLameMu();
      rDerivative /= -rVariables.Strain.Invariants.I3;
      rDerivative += rMaterial.GetLameLambda();
      rDerivative *= 0.25;

      return rDerivative;

      KRATOS_CATCH(" ")
    }


    double& GetFunction2ndI1Derivative(HyperElasticDataType& rVariables, double& rDerivative) override //ddW/dI1dI1
    {
      KRATOS_TRY

      rDerivative = 0.0;

      return rDerivative;

      KRATOS_CATCH(" ")
    }

    double& GetFunction2ndI2Derivative(HyperElasticDataType& rVariables, double& rDerivative) override //ddW/dI2dI2
    {
      KRATOS_TRY

      rDerivative = 0.0;

      return rDerivative;

      KRATOS_CATCH(" ")
    }

    double& GetFunction2ndI3Derivative(HyperElasticDataType& rVariables, double& rDerivative) override //ddW/dI3dI3
    {
      KRATOS_TRY

      const MaterialDataType&   rMaterial = rVariables.GetMaterialParameters();

      //ddg(J)/dI3dI3 = (lambda/4)*(1/J⁴) + (mu/2)*(1/J⁴)
      rDerivative  = 0.25 * rMaterial.GetLameLambda();
      rDerivative += 0.5  * rMaterial.GetLameMu();
      rDerivative /= (rVariables.Strain.Invariants.I3 * rVariables.Strain.Invariants.I3);

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


    void save(Serializer& rSerializer) const override
    {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, MooneyRivlinModel )
    }

    void load(Serializer& rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, MooneyRivlinModel )
    }

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

  }; // Class NeoHookeanModel

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_NEO_HOOKEAN_MODEL_H_INCLUDED  defined
