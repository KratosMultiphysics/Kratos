//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_ISOCHORIC_NEO_HOOKEAN_MODEL_H_INCLUDED )
#define  KRATOS_ISOCHORIC_NEO_HOOKEAN_MODEL_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_models/elasticity_models/isochoric_mooney_rivlin_model.hpp"

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
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) IsochoricNeoHookeanModel : public IsochoricMooneyRivlinModel
  {
  public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of IsochoricNeoHookeanModel
    KRATOS_CLASS_POINTER_DEFINITION( IsochoricNeoHookeanModel );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    IsochoricNeoHookeanModel() : IsochoricMooneyRivlinModel() {}

    /// Copy constructor.
    IsochoricNeoHookeanModel(IsochoricNeoHookeanModel const& rOther) : IsochoricMooneyRivlinModel(rOther) {}

    /// Assignment operator.
    IsochoricNeoHookeanModel& operator=(IsochoricNeoHookeanModel const& rOther)
    {
      IsochoricMooneyRivlinModel::operator=(rOther);
      return *this;
    }

    /// Clone.
    ConstitutiveModel::Pointer Clone() const override
    {
      return Kratos::make_shared<IsochoricNeoHookeanModel>(*this);
    }

    /// Destructor.
    ~IsochoricNeoHookeanModel() override {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    // Simplyfied methods must be implemented for performance purposes

    /**
     * Calculate Stresses
     */



    /**
     * Calculate Constitutive Components
     */



    /**
     * Check
     */

    int Check(const Properties& rMaterialProperties, const ProcessInfo& rCurrentProcessInfo) override
    {
      KRATOS_TRY

      HyperElasticModel::Check(rMaterialProperties,rCurrentProcessInfo);

      if( C10.Key() == 0 || rMaterialProperties[C10] <= 0.00 )
	KRATOS_ERROR << "C10 has an invalid key or value" << std::endl;

      if( BULK_MODULUS.Key() == 0 || rMaterialProperties[BULK_MODULUS] <= 0.00 )
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
        buffer << "IsochoricNeoHookeanModel";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "IsochoricNeoHookeanModel";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "IsochoricNeoHookeanModel Data";
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

    //specialized methods:

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

    // virtual void CalculateAndAddIsochoricStressTensor(HyperElasticDataType& rVariables, MatrixType& rStressMatrix) override
    // {
    //   KRATOS_TRY

    //   const ModelDataType&  rModelData        = rVariables.GetModelData();
    //   const StressMeasureType& rStressMeasure = rModelData.GetStressMeasure();

    //   MatrixType StressMatrix;
    //   const MaterialDataType& rMaterial = rVariables.GetMaterialParameters();

    //   if( rStressMeasure == ConstitutiveModelData::StressMeasure_PK2 ){ //Variables.Strain.Matrix = RightCauchyGreen (C)

    // 	StressMatrix  = msIdentityMatrix;
    // 	StressMatrix -= 1.0/3.0 * ( rVariables.Strain.Matrix(0,0) + rVariables.Strain.Matrix(1,1) + rVariables.Strain.Matrix(2,2) ) * rVariables.Strain.InverseMatrix;

    // 	StressMatrix *= rMaterial.GetLameMu() * rVariables.Strain.Invariants.J_13 * rVariables.Strain.Invariants.J_13;

    // 	rStressMatrix += StressMatrix;
    //   }
    //   else if( rStressMeasure == ConstitutiveModelData::StressMeasure_Kirchhoff ){ //Variables.Strain.Matrix = LeftCauchyGreen (b)

    // 	StressMatrix  = rVariables.Strain.Matrix;
    // 	StressMatrix -= 1.0/3.0 * ( rVariables.Strain.Matrix(0,0) + rVariables.Strain.Matrix(1,1) + rVariables.Strain.Matrix(2,2) ) * msIdentityMatrix;
    // 	StressMatrix *= rMaterial.GetLameMu() * rVariables.Strain.Invariants.J_13 * rVariables.Strain.Invariants.J_13;

    // 	rStressMatrix += StressMatrix;
    //   }


    //   KRATOS_CATCH(" ")
    // }

    // virtual void CalculateAndAddVolumetricStressTensor(HyperElasticDataType& rVariables, MatrixType& rStressMatrix) override
    // {
    //   KRATOS_TRY

    //   const ModelDataType&  rModelData        = rVariables.GetModelData();
    //   const StressMeasureType& rStressMeasure = rModelData.GetStressMeasure();


    //   MatrixType StressMatrix;
    //   const MaterialDataType& rMaterial = rVariables.GetMaterialParameters();

    //   double Factor = 0;
    //   this->CalculateVolumetricFactor(rVariables,Factor);


    //   if( rStressMeasure == ConstitutiveModelData::StressMeasure_PK2 ){ //Variables.Strain.Matrix = RightCauchyGreen (C)

    // 	StressMatrix = rMaterial.GetBulkModulus() * Factor * rVariables.Strain.InverseMatrix;

    // 	rStressMatrix += StressMatrix;
    //   }
    //   else if( rStressMeasure == ConstitutiveModelData::StressMeasure_Kirchhoff ){ //Variables.Strain.Matrix = LeftCauchyGreen (b)

    // 	StressMatrix = rMaterial.GetBulkModulus() * Factor * msIdentityMatrix;

    // 	rStressMatrix += StressMatrix;
    //   }


    //   KRATOS_CATCH(" ")
    // }


    // virtual double& AddIsochoricConstitutiveComponent(HyperElasticDataType& rVariables, double &rCabcd,
    // 						      const unsigned int& a, const unsigned int& b,
    // 						      const unsigned int& c, const unsigned int& d) override
    // {
    //   KRATOS_TRY


    //   double Cabcd = 0;

    //   const MaterialDataType& rMaterial = rVariables.GetMaterialParameters();

    //   const ModelDataType&  rModelData         = rVariables.GetModelData();
    //   const StressMeasureType& rStressMeasure  = rModelData.GetStressMeasure();
    //   const MatrixType& rIsochoricStressMatrix = rModelData.GetStressMatrix();


    //   if( rStressMeasure == ConstitutiveModelData::StressMeasure_PK2 ){ //mStrainMatrix = RightCauchyGreen (C)

    // 	Cabcd  = (1.0/3.0) * (rVariables.Strain.InverseMatrix(a,b)*rVariables.Strain.InverseMatrix(d,c));

    // 	Cabcd -= 0.5 * (rVariables.Strain.InverseMatrix(a,c)*rVariables.Strain.InverseMatrix(b,d)+rVariables.Strain.InverseMatrix(a,d)*rVariables.Strain.InverseMatrix(b,c));

    // 	Cabcd *= rMaterial.GetLameMu() * ( rVariables.Strain.Matrix(0,0) + rVariables.Strain.Matrix(1,1) + rVariables.Strain.Matrix(2,2) ) * rVariables.Strain.Invariants.J_13 * rVariables.Strain.Invariants.J_13;

    // 	Cabcd += (rVariables.Strain.InverseMatrix(c,d)*rIsochoricStressMatrix(a,b)+rIsochoricStressMatrix(c,d)*rVariables.Strain.InverseMatrix(a,b));

    // 	Cabcd *= (-2.0/3.0);

    //   }
    //   else if( rStressMeasure == ConstitutiveModelData::StressMeasure_Kirchhoff ){ //mStrainMatrix = LeftCauchyGreen (b)

    // 	Cabcd  = (1.0/3.0) * (msIdentityMatrix(a,b)*msIdentityMatrix(c,d));


    // 	Cabcd -= 0.5 * (msIdentityMatrix(a,c)*msIdentityMatrix(b,d)+msIdentityMatrix(a,d)*msIdentityMatrix(b,c));


    // 	Cabcd *= rMaterial.GetLameMu() * ( rVariables.Strain.Matrix(0,0) + rVariables.Strain.Matrix(1,1) + rVariables.Strain.Matrix(2,2) ) * rVariables.Strain.Invariants.J_13 * rVariables.Strain.Invariants.J_13;


    // 	Cabcd +=  (msIdentityMatrix(c,d)*rIsochoricStressMatrix(a,b)+rIsochoricStressMatrix(c,d)*msIdentityMatrix(a,b));

    // 	Cabcd *= (-2.0/3.0);

    //   }

    //   rCabcd += Cabcd;

    //   rVariables.State().Set(ConstitutiveModelData::CONSTITUTIVE_MATRIX_COMPUTED);

    //   return rCabcd;

    //   KRATOS_CATCH(" ")
    // }


    // virtual double& AddVolumetricConstitutiveComponent(HyperElasticDataType& rVariables, double &rCabcd,
    // 						       const unsigned int& a, const unsigned int& b,
    // 						       const unsigned int& c, const unsigned int& d) override
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


    //   if( rStressMeasure == ConstitutiveModelData::StressMeasure_PK2 ){ //mStrainMatrix = RightCauchyGreen (C)

    // 	Cabcd  = FactorA * (rVariables.Strain.InverseMatrix(a,b)*rVariables.Strain.InverseMatrix(c,d));

    // 	Cabcd -= FactorB * (rVariables.Strain.InverseMatrix(a,c)*rVariables.Strain.InverseMatrix(b,d)+rVariables.Strain.InverseMatrix(a,d)*rVariables.Strain.InverseMatrix(b,c));

    // 	Cabcd *= rMaterial.GetBulkModulus();

    //   }
    //   else if( rStressMeasure == ConstitutiveModelData::StressMeasure_Kirchhoff ){ //mStrainMatrix = LeftCauchyGreen (b)

    // 	Cabcd  = FactorA * (msIdentityMatrix(a,b)*msIdentityMatrix(c,d));

    // 	Cabcd -= FactorB * (msIdentityMatrix(a,c)*msIdentityMatrix(b,d)+msIdentityMatrix(a,d)*msIdentityMatrix(b,c));

    // 	Cabcd *= rMaterial.GetBulkModulus();

    //   }

    //   rCabcd += Cabcd;

    //   rVariables.State().Set(ConstitutiveModelData::CONSTITUTIVE_MATRIX_COMPUTED);

    //   return rCabcd;

    //   KRATOS_CATCH(" ")
    // }

    //************// W

    void CalculateAndAddIsochoricStrainEnergy(HyperElasticDataType& rVariables, double& rIsochoricDensityFunction) override
    {
      KRATOS_TRY

      const MaterialDataType& rMaterial = rVariables.GetMaterialParameters();

      rIsochoricDensityFunction += rMaterial.GetModelParameters()[0] * ( rVariables.Strain.Invariants.J_13 * rVariables.Strain.Invariants.I1 - 3.0);

      KRATOS_CATCH(" ")
    }


    void CalculateAndAddVolumetricStrainEnergy(HyperElasticDataType& rVariables, double& rVolumetricDensityFunction) override
    {
      KRATOS_TRY

      const MaterialDataType& rMaterial = rVariables.GetMaterialParameters();

      //energy function "U(J) = (K/4)*(J²-1) - (K/2)*lnJ"
      rVolumetricDensityFunction += rMaterial.GetBulkModulus() * 0.25 * ( rVariables.Strain.Invariants.J * rVariables.Strain.Invariants.J - 1.0);
      rVolumetricDensityFunction -= rMaterial.GetBulkModulus() * 0.5 * std::log( rVariables.Strain.Invariants.J );

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

      rDerivative = 0.0;

      return rDerivative;

      KRATOS_CATCH(" ")
    }

    double& GetFunction1stI3Derivative(HyperElasticDataType& rVariables, double& rDerivative) override //dW/dI3
    {
      KRATOS_TRY

      rDerivative = 0.0;

      return rDerivative;

      KRATOS_CATCH(" ")
    }


    double& GetVolumetricFunction1stJDerivative(HyperElasticDataType& rVariables, double& rDerivative) override //dU/dJ
    {
      KRATOS_TRY

      const MaterialDataType& rMaterial = rVariables.GetMaterialParameters();

      //derivative of "U(J) = (K/4)*(J²-1) - (K/2)*lnJ"
      //dU(J)/dJ = (K/2)*(J-1/J)
      rDerivative = 0.5 * rMaterial.GetBulkModulus() * ( rVariables.Strain.Invariants.J * rVariables.Strain.Invariants.J - 1.0 );

      rDerivative /= rVariables.Strain.Invariants.J;

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

      rDerivative = 0.0;

      return rDerivative;

      KRATOS_CATCH(" ")
    }


    double& GetVolumetricFunction2ndJDerivative(HyperElasticDataType& rVariables, double& rDerivative) override //ddU/dJdJ
    {
      KRATOS_TRY

      const MaterialDataType& rMaterial = rVariables.GetMaterialParameters();

      //derivative of "dU(J)/dJ = (K/2)*(J-1/J)"
      //ddU(J)/dJdJ = (K/2)*(1-1/J²)
      rDerivative = 0.5 * rMaterial.GetBulkModulus() * (rVariables.Strain.Invariants.J * rVariables.Strain.Invariants.J + 1.0 );

      rDerivative /= rVariables.Strain.Invariants.J * rVariables.Strain.Invariants.J ;

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


    void save(Serializer& rSerializer) const  override
    {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, IsochoricMooneyRivlinModel )
    }

    void load(Serializer& rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, IsochoricMooneyRivlinModel )
    }

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

  }; // Class IsochoricNeoHookeanModel

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_ISOCHORIC_NEO_HOOKEAN_MODEL_H_INCLUDED  defined
