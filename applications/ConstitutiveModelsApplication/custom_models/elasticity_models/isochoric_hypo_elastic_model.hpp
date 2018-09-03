//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2018 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_ISOCHORIC_HYPO_ELASTIC_MODEL_H_INCLUDED )
#define  KRATOS_ISOCHORIC_HYPO_ELASTIC_MODEL_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_models/elasticity_models/hypo_elastic_model.hpp"

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
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) IsochoricHypoElasticModel : public HypoElasticModel
  {
  public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of IsochoricHypoElasticModel
    KRATOS_CLASS_POINTER_DEFINITION( IsochoricHypoElasticModel );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    IsochoricHypoElasticModel() : HypoElasticModel() {}

    /// Copy constructor.
    IsochoricHypoElasticModel(IsochoricHypoElasticModel const& rOther) : HypoElasticModel(rOther) {}

    /// Assignment operator.
    IsochoricHypoElasticModel& operator=(IsochoricHypoElasticModel const& rOther)
    {
	HypoElasticModel::operator=(rOther);
	return *this;
    }

    /// Clone.
    ConstitutiveModel::Pointer Clone() const override
    {
      return Kratos::make_shared<IsochoricHypoElasticModel>(*this);
    }

    /// Destructor.
    ~IsochoricHypoElasticModel() override {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    void CalculateStrainEnergy(ModelDataType& rValues, double& rDensityFunction) override
    {
      KRATOS_TRY

      ElasticDataType Variables;
      this->InitializeElasticData(rValues,Variables);

      rDensityFunction = 0;
      this->CalculateAndAddIsochoricStrainEnergy( Variables, rDensityFunction );
      this->CalculateAndAddVolumetricStrainEnergy( Variables, rDensityFunction );


      KRATOS_CATCH(" ")
    }


    void CalculateStressTensor(ModelDataType& rValues, MatrixType& rStressMatrix) override
    {
      KRATOS_TRY

      ElasticDataType Variables;
      this->InitializeElasticData(rValues,Variables);

      VectorType StrainVector;
      ConstitutiveModelUtilities::StrainTensorToVector(Variables.StrainMatrix, StrainVector);

      this->CalculateAndAddConstitutiveTensor(Variables);

      VectorType StressVector;
      noalias(StressVector) = ZeroVector(6);
      this->CalculateAndAddIsochoricStressTensor(Variables,StrainVector,StressVector);

      rValues.StressMatrix = ConstitutiveModelUtilities::VectorToSymmetricTensor(StressVector,rValues.StressMatrix); //store isochoric stress matrix as StressMatrix

      this->CalculateAndAddVolumetricStressTensor(Variables,StrainVector,StressVector);

      rStressMatrix = ConstitutiveModelUtilities::VectorToSymmetricTensor(StressVector,rStressMatrix);

      // Rate to stress value
      rStressMatrix *= rValues.GetProcessInfo()[DELTA_TIME];

      // Isochoric stress stored in the HistoryVector
      this->AddHistoricalStress(rValues, rStressMatrix);

      Variables.State().Set(ConstitutiveModelData::STRESS_COMPUTED);

      KRATOS_CATCH(" ")
    }


    void CalculateStressAndConstitutiveTensors(ModelDataType& rValues, MatrixType& rStressMatrix, Matrix& rConstitutiveMatrix) override
    {
      KRATOS_TRY

      ElasticDataType Variables;
      this->InitializeElasticData(rValues,Variables);

      VectorType StrainVector;
      ConstitutiveModelUtilities::StrainTensorToVector(Variables.StrainMatrix, StrainVector);

      //Calculate Constitutive Matrix
      this->CalculateAndAddConstitutiveTensor(Variables,rConstitutiveMatrix);

      //Calculate Stress Matrix
      VectorType StressVector;
      noalias(StressVector) = ZeroVector(6);
      this->CalculateAndAddIsochoricStressTensor(Variables,StrainVector,StressVector);

      rValues.StressMatrix = ConstitutiveModelUtilities::VectorToSymmetricTensor(StressVector,rValues.StressMatrix); //store isochoric stress matrix as StressMatrix

      this->CalculateAndAddVolumetricStressTensor(Variables,StrainVector,StressVector);

      rStressMatrix = ConstitutiveModelUtilities::VectorToSymmetricTensor(StressVector,rStressMatrix);

      // Rate to stress value
      rStressMatrix *= rValues.GetProcessInfo()[DELTA_TIME];

      // Isochoric stress stored in the HistoryVector
      this->AddHistoricalStress(rValues, rStressMatrix);

      Variables.State().Set(ConstitutiveModelData::STRESS_COMPUTED);

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
        buffer << "IsochoricHypoElasticModel";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "IsochoricHypoElasticModel";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "IsochoricHypoElasticModel Data";
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


    void CalculateAndAddIsochoricStressTensor(ElasticDataType& rVariables, VectorType& rStrainVector, VectorType& rStressVector) override
    {
      KRATOS_TRY

      //total stress
      VectorType StressVector;
      this->CalculateAndAddStressTensor(rVariables,rStrainVector,StressVector);

      //deviatoric stress
      double MeanStress = (1.0/3.0) * (StressVector[0]+StressVector[1]+StressVector[2]);
      for (unsigned int i = 0; i < 3; i++)
        StressVector[i] -= MeanStress;

      rStressVector += StressVector;

      KRATOS_CATCH(" ")
    }


    void CalculateAndAddVolumetricStressTensor(ElasticDataType& rVariables, VectorType& rStrainVector, VectorType& rStressVector) override
    {
      KRATOS_TRY

      //total stress
      VectorType StressVector;
      this->CalculateAndAddStressTensor(rVariables,rStrainVector,StressVector);

      //volumetric stress
      double MeanStress = (1.0/3.0) * (StressVector[0]+StressVector[1]+StressVector[2]);
      for (unsigned int i = 0; i < 3; i++)
        rStressVector[i] += MeanStress;

      KRATOS_CATCH(" ")
    }


    void CalculateAndAddIsochoricStrainEnergy(ElasticDataType& rVariables, double& rIsochoricDensityFunction) override
    {
      KRATOS_TRY

      KRATOS_ERROR << "calling the class function in IsochoricHypoElasticModel ... illegal operation" << std::endl;

      KRATOS_CATCH(" ")
    }

    void CalculateAndAddVolumetricStrainEnergy(ElasticDataType& rVariables, double& rVolumetricDensityFunction) override
    {
      KRATOS_TRY

      KRATOS_ERROR << "calling the class function in IsochoricHypoElasticModel ... illegal operation" << std::endl;

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
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, HypoElasticModel )
    }

    void load(Serializer& rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, HypoElasticModel )
    }

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

  }; // Class IsochoricHypoElasticModel

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_ISOCHORIC_HYPO_ELASTIC_MODEL_H_INCLUDED  defined
