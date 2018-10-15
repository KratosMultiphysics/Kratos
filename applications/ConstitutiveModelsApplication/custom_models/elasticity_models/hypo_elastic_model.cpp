//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2018 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_models/elasticity_models/hypo_elastic_model.hpp"

namespace Kratos
{

  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

  HypoElasticModel::HypoElasticModel()
    : ConstitutiveModel()
  {
  }

  //******************************COPY CONSTRUCTOR**************************************
  //************************************************************************************

  HypoElasticModel::HypoElasticModel(const HypoElasticModel& rOther)
    : ConstitutiveModel(rOther)
  {
  }

  //********************************CLONE***********************************************
  //************************************************************************************

  ConstitutiveModel::Pointer HypoElasticModel::Clone() const
  {
    return Kratos::make_shared<HypoElasticModel>(*this);
  }

  //********************************ASSIGNMENT******************************************
  //************************************************************************************
  HypoElasticModel& HypoElasticModel::operator=(HypoElasticModel const& rOther)
  {
    ConstitutiveModel::operator=(rOther);
    return *this;
  }

  //*******************************DESTRUCTOR*******************************************
  //************************************************************************************

  HypoElasticModel::~HypoElasticModel()
  {
  }

  //***********************PUBLIC OPERATIONS FROM BASE CLASS****************************
  //************************************************************************************

  void HypoElasticModel::InitializeModel(ModelDataType& rValues)
  {
    KRATOS_TRY


    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************

  void HypoElasticModel::FinalizeModel(ModelDataType& rValues)
  {
    KRATOS_TRY


    //update total stress measure
    ConstitutiveModelUtilities::SymmetricTensorToVector(rValues.StressMatrix, this->mHistoryVector);

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************

  void HypoElasticModel::InitializeElasticData(ModelDataType& rValues, ElasticDataType& rVariables)
  {
    KRATOS_TRY

    //set model data pointer
    rVariables.SetModelData(rValues);
    rVariables.SetState(rValues.State);

    // ConstitutiveModelData::StressMeasureType::StressMeasure_Kirchhoff  allowed only

    // symmetric spatial velocity gradient
    noalias(rVariables.StrainMatrix) = 0.5 * (rValues.StrainMatrix + trans(rValues.StrainMatrix)); // spatial velocity gradient is rValues.StrainMatrix

    rValues.SetStrainMeasure(ConstitutiveModelData::StrainMeasureType::CauchyGreen_None);
    rValues.MaterialParameters.LameMuBar = rValues.MaterialParameters.LameMu;

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************

  void HypoElasticModel::CalculateStrainEnergy(ModelDataType& rValues, double& rDensityFunction)
  {
    KRATOS_TRY

    KRATOS_ERROR << "calling the base class function in HypoElasticModel ... illegal operation" << std::endl;

    KRATOS_CATCH(" ")
  }

  //************// W

  void HypoElasticModel::CalculateAndAddIsochoricStrainEnergy(ElasticDataType& rVariables, double& rIsochoricDensityFunction)
  {
    KRATOS_TRY

    KRATOS_ERROR << "calling the base class function in HypoElasticModel ... illegal operation" << std::endl;

    KRATOS_CATCH(" ")
  }


  void HypoElasticModel::CalculateAndAddVolumetricStrainEnergy(ElasticDataType& rVariables, double& rVolumetricDensityFunction)
  {
    KRATOS_TRY

    KRATOS_ERROR << "calling the base class function in HypoElasticModel ... illegal operation" << std::endl;

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************

  void HypoElasticModel::AddHistoricalStress(ModelDataType& rValues, MatrixType& rStressMatrix)
  {
    KRATOS_TRY


    // Add historical stress
    const MatrixType& rSpatialVelocityGradient = rValues.GetStrainMatrix(); // spatial velocity gradient is rValues.StrainMatrix
    const double& rDeltaTime                   = rValues.GetProcessInfo()[DELTA_TIME];

    // Skewsymmetric spatial velocity gradient
    Matrix W = 0.5 * rDeltaTime * (rSpatialVelocityGradient - trans(rSpatialVelocityGradient));

    // Exponential map using quaternions
    Quaternion<double> QuaternionValue = Quaternion<double>::FromRotationMatrix( W );
    QuaternionValue.ToRotationMatrix(W);


    MatrixType PreviousStressMatrix;
    PreviousStressMatrix = ConstitutiveModelUtilities::VectorToSymmetricTensor(this->mHistoryVector, PreviousStressMatrix);

    PreviousStressMatrix = prod( W, PreviousStressMatrix );
    PreviousStressMatrix = prod( PreviousStressMatrix, trans(W) );

    rStressMatrix += PreviousStressMatrix;

    // To store when FinalizeModel (can be the total or the isochoric stress)
    rValues.StressMatrix = rStressMatrix;

    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************

  void HypoElasticModel::CalculateStressTensor(ModelDataType& rValues, MatrixType& rStressMatrix)
  {
    KRATOS_TRY

    ElasticDataType Variables;
    this->InitializeElasticData(rValues,Variables);

    VectorType StrainVector;
    ConstitutiveModelUtilities::StrainTensorToVector(Variables.StrainMatrix, StrainVector);

    this->CalculateAndAddConstitutiveTensor(Variables);

    VectorType StressVector;
    this->CalculateAndAddStressTensor(Variables,StrainVector,StressVector);

    rStressMatrix = ConstitutiveModelUtilities::VectorToSymmetricTensor(StressVector,rStressMatrix);

    // Rate to stress value
    rStressMatrix *= rValues.GetProcessInfo()[DELTA_TIME];

    // Total stress stored in the HistoryVector
    this->AddHistoricalStress(rValues, rStressMatrix);

    Variables.State().Set(ConstitutiveModelData::STRESS_COMPUTED);

    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************

  void HypoElasticModel::CalculateAndAddStressTensor(ElasticDataType& rVariables, VectorType& rStrainVector, VectorType& rStressVector)
  {
    KRATOS_TRY

    noalias(rStressVector) = prod(rVariables.ConstitutiveTensor,rStrainVector);

    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************

  void HypoElasticModel::CalculateIsochoricStressTensor(ModelDataType& rValues, MatrixType& rStressMatrix)
  {
    KRATOS_TRY

    ElasticDataType Variables;
    this->InitializeElasticData(rValues,Variables);

    VectorType StrainVector;
    ConstitutiveModelUtilities::StrainTensorToVector(Variables.StrainMatrix, StrainVector);

    this->CalculateAndAddConstitutiveTensor(Variables);

    VectorType StressVector;
    StressVector.clear();
    this->CalculateAndAddIsochoricStressTensor(Variables,StrainVector,StressVector);

    rStressMatrix = ConstitutiveModelUtilities::VectorToSymmetricTensor(StressVector,rStressMatrix);

    // Rate to stress value
    rStressMatrix *= rValues.GetProcessInfo()[DELTA_TIME];

    // Isochoric stress stored in the HistoryVector
    this->AddHistoricalStress(rValues, rStressMatrix);

    Variables.State().Set(ConstitutiveModelData::STRESS_COMPUTED);

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************

  void HypoElasticModel::CalculateAndAddIsochoricStressTensor(ElasticDataType& rVariables, VectorType& rStrainVector, VectorType& rStressVector)
  {
    KRATOS_TRY

    KRATOS_ERROR << "calling the class function in HypoElasticModel ... illegal operation" << std::endl;

    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************

  void HypoElasticModel::CalculateVolumetricStressTensor(ModelDataType& rValues, MatrixType& rStressMatrix)
  {
    KRATOS_TRY

    KRATOS_ERROR << "calling the class function in HypoElasticModel ... illegal operation" << std::endl;

    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************

  void HypoElasticModel::CalculateAndAddVolumetricStressTensor(ElasticDataType& rVariables, VectorType& rStrainVector, VectorType& rStressVector)
  {
    KRATOS_TRY

    KRATOS_ERROR << "calling the base class function in HypoElasticModel ... illegal operation" << std::endl;

    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************

  void HypoElasticModel::CalculateConstitutiveTensor(ModelDataType& rValues, Matrix& rConstitutiveMatrix)
  {
    KRATOS_TRY

    ElasticDataType Variables;
    this->InitializeElasticData(rValues,Variables);

    //Set constitutive matrix to zero before adding
    rConstitutiveMatrix.clear();

    this->CalculateAndAddConstitutiveTensor(Variables,rConstitutiveMatrix);

    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************

  void HypoElasticModel::CalculateAndAddConstitutiveTensor(ElasticDataType& rVariables, Matrix& rConstitutiveMatrix)
  {
    KRATOS_TRY

    this->CalculateAndAddConstitutiveTensor(rVariables);

    rConstitutiveMatrix = ConstitutiveModelUtilities::ConstitutiveTensorToMatrix(rVariables.ConstitutiveTensor,rConstitutiveMatrix);

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************

  void HypoElasticModel::CalculateAndAddConstitutiveTensor(ElasticDataType& rVariables)
  {
    KRATOS_TRY

    //Calculate Elastic ConstitutiveMatrix
    const ModelDataType&  rModelData  = rVariables.GetModelData();
    const MaterialDataType& rMaterial = rModelData.GetMaterialParameters();

    rVariables.ConstitutiveTensor.clear();

    // Lame constants
    const double& rYoungModulus       = rMaterial.GetYoungModulus();
    const double& rPoissonCoefficient = rMaterial.GetPoissonCoefficient();

    // 3D linear elastic constitutive matrix
    rVariables.ConstitutiveTensor ( 0 , 0 ) = (rYoungModulus*(1.0-rPoissonCoefficient)/((1.0+rPoissonCoefficient)*(1.0-2.0*rPoissonCoefficient)));
    rVariables.ConstitutiveTensor ( 1 , 1 ) = rVariables.ConstitutiveTensor ( 0 , 0 );
    rVariables.ConstitutiveTensor ( 2 , 2 ) = rVariables.ConstitutiveTensor ( 0 , 0 );

    rVariables.ConstitutiveTensor ( 3 , 3 ) = rVariables.ConstitutiveTensor ( 0 , 0 )*(1.0-2.0*rPoissonCoefficient)/(2.0*(1.0-rPoissonCoefficient));
    rVariables.ConstitutiveTensor ( 4 , 4 ) = rVariables.ConstitutiveTensor ( 3 , 3 );
    rVariables.ConstitutiveTensor ( 5 , 5 ) = rVariables.ConstitutiveTensor ( 3 , 3 );

    rVariables.ConstitutiveTensor ( 0 , 1 ) = rVariables.ConstitutiveTensor ( 0 , 0 )*rPoissonCoefficient/(1.0-rPoissonCoefficient);
    rVariables.ConstitutiveTensor ( 1 , 0 ) = rVariables.ConstitutiveTensor ( 0 , 1 );

    rVariables.ConstitutiveTensor ( 0 , 2 ) = rVariables.ConstitutiveTensor ( 0 , 1 );
    rVariables.ConstitutiveTensor ( 2 , 0 ) = rVariables.ConstitutiveTensor ( 0 , 1 );

    rVariables.ConstitutiveTensor ( 1 , 2 ) = rVariables.ConstitutiveTensor ( 0 , 1 );
    rVariables.ConstitutiveTensor ( 2 , 1 ) = rVariables.ConstitutiveTensor ( 0 , 1 );


    rVariables.State().Set(ConstitutiveModelData::CONSTITUTIVE_MATRIX_COMPUTED);


    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************

  void HypoElasticModel::CalculateStressAndConstitutiveTensors(ModelDataType& rValues, MatrixType& rStressMatrix, Matrix& rConstitutiveMatrix)
  {
    KRATOS_TRY

    ElasticDataType Variables;
    this->InitializeElasticData(rValues,Variables);

    VectorType StrainVector;
    ConstitutiveModelUtilities::StrainTensorToVector(Variables.StrainMatrix, StrainVector);

    //Set constitutive matrix to zero before adding
    rConstitutiveMatrix.clear();

    this->CalculateAndAddConstitutiveTensor(Variables, rConstitutiveMatrix);

    VectorType StressVector;
    this->CalculateAndAddStressTensor(Variables,StrainVector,StressVector);

    rStressMatrix = ConstitutiveModelUtilities::VectorToSymmetricTensor(StressVector,rStressMatrix);

    // Rate to stress value
    rStressMatrix *= rValues.GetProcessInfo()[DELTA_TIME];

    // Total stress stored in the HistoryVector
    this->AddHistoricalStress(rValues, rStressMatrix);

    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************

  void HypoElasticModel::CalculateIsochoricStressAndConstitutiveTensors(ModelDataType& rValues, MatrixType& rStressMatrix, Matrix& rConstitutiveMatrix)
  {
    KRATOS_TRY

    KRATOS_ERROR << "calling the base class function in HypoElasticModel ... illegal operation" << std::endl;

    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************

  void HypoElasticModel::CalculateVolumetricStressAndConstitutiveTensors(ModelDataType& rValues, MatrixType& rStressMatrix, Matrix& rConstitutiveMatrix)
  {
    KRATOS_TRY

    KRATOS_ERROR << "calling the base class function in HypoElasticModel ... illegal operation" << std::endl;

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************

  void HypoElasticModel::CalculateIsochoricConstitutiveTensor(ModelDataType& rValues, Matrix& rConstitutiveMatrix)
  {
    KRATOS_TRY

    ElasticDataType Variables;
    this->InitializeElasticData(rValues,Variables);
    this->CalculateAndAddIsochoricConstitutiveTensor(Variables,rConstitutiveMatrix);

    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************

  void HypoElasticModel::CalculateAndAddIsochoricConstitutiveTensor(ElasticDataType& rVariables, Matrix& rConstitutiveMatrix)
  {
    KRATOS_TRY

    this->CalculateAndAddIsochoricConstitutiveTensor(rVariables);

    rConstitutiveMatrix = ConstitutiveModelUtilities::ConstitutiveTensorToMatrix(rVariables.ConstitutiveTensor, rConstitutiveMatrix);

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************

  void HypoElasticModel::CalculateAndAddIsochoricConstitutiveTensor(ElasticDataType& rVariables)
  {
    KRATOS_TRY

    KRATOS_ERROR << "calling the base class function in HypoElasticModel ... illegal operation" << std::endl;

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************

  void HypoElasticModel::CalculateVolumetricConstitutiveTensor(ModelDataType& rValues, Matrix& rConstitutiveMatrix)
  {
    KRATOS_TRY

    ElasticDataType Variables;
    this->InitializeElasticData(rValues,Variables);
    this->CalculateAndAddVolumetricConstitutiveTensor(Variables,rConstitutiveMatrix);

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************

  void HypoElasticModel::CalculateAndAddVolumetricConstitutiveTensor(ElasticDataType& rVariables, Matrix& rConstitutiveMatrix)
  {
    KRATOS_TRY

    this->CalculateAndAddVolumetricConstitutiveTensor(rVariables);

    rConstitutiveMatrix = ConstitutiveModelUtilities::ConstitutiveTensorToMatrix(rVariables.ConstitutiveTensor, rConstitutiveMatrix);

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************

  void HypoElasticModel::CalculateAndAddVolumetricConstitutiveTensor(ElasticDataType& rVariables)
  {
    KRATOS_TRY

    KRATOS_ERROR << "calling the base class function in HypoElasticModel ... illegal operation" << std::endl;

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************

  int HypoElasticModel::Check(const Properties& rProperties, const ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    if(YOUNG_MODULUS.Key() == 0 || rProperties[YOUNG_MODULUS] <= 0.00)
      KRATOS_ERROR << "YOUNG_MODULUS has Key zero or invalid value" << std::endl;

    if(POISSON_RATIO.Key() == 0){
      KRATOS_ERROR << "POISSON_RATIO has Key zero invalid value" << std::endl;
    }
    else{
      const double& nu = rProperties[POISSON_RATIO];
      if( (nu > 0.499 && nu < 0.501) || (nu < -0.999 && nu > -1.01) )
	KRATOS_ERROR << "POISSON_RATIO has an invalid value" << std::endl;
    }

    return 0;


    KRATOS_CATCH(" ")
  }



} // Namespace Kratos
