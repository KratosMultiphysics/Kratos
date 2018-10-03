//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_models/elasticity_models/linear_elastic_model.hpp"


namespace Kratos
{

  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

  LinearElasticModel::LinearElasticModel()
    : ConstitutiveModel()
  {
  }

  //******************************COPY CONSTRUCTOR**************************************
  //************************************************************************************

  LinearElasticModel::LinearElasticModel(const LinearElasticModel& rOther)
    : ConstitutiveModel(rOther)
  {
  }

  //********************************CLONE***********************************************
  //************************************************************************************

  ConstitutiveModel::Pointer LinearElasticModel::Clone() const
  {
    return Kratos::make_shared<LinearElasticModel>(*this);
  }

  //********************************ASSIGNMENT******************************************
  //************************************************************************************
  LinearElasticModel& LinearElasticModel::operator=(LinearElasticModel const& rOther)
  {
    ConstitutiveModel::operator=(rOther);
    return *this;
  }

  //*******************************DESTRUCTOR*******************************************
  //************************************************************************************

  LinearElasticModel::~LinearElasticModel()
  {
  }

  //***********************PUBLIC OPERATIONS FROM BASE CLASS****************************
  //************************************************************************************

  void LinearElasticModel::InitializeModel(ModelDataType& rValues)
  {
    KRATOS_TRY


    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************

  void LinearElasticModel::FinalizeModel(ModelDataType& rValues)
  {
    KRATOS_TRY


    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************

  void LinearElasticModel::InitializeElasticData(ModelDataType& rValues, ElasticDataType& rVariables)
  {
    KRATOS_TRY

    //set model data pointer
    rVariables.SetModelData(rValues);
    rVariables.SetState(rValues.State);

    //add initial strain
    if(this->mOptions.Is(ConstitutiveModel::ADD_HISTORY_VECTOR) && this->mOptions.Is(ConstitutiveModel::HISTORY_STRAIN_MEASURE) ){
       VectorType StrainVector;
       ConstitutiveModelUtilities::StrainTensorToVector(rValues.StrainMatrix, StrainVector);
       for(unsigned int i=0; i<StrainVector.size(); i++)
       {
          StrainVector[i] += this->mHistoryVector[i];
       }
       rValues.StrainMatrix = ConstitutiveModelUtilities::StrainVectorToTensor(StrainVector, rValues.StrainMatrix);
    }

    rValues.SetStrainMeasure( ConstitutiveModelData::CauchyGreen_None);
    rValues.MaterialParameters.LameMuBar = rValues.MaterialParameters.LameMu;

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************

  void LinearElasticModel::CalculateStrainEnergy(ModelDataType& rValues, double& rDensityFunction)
  {
    KRATOS_TRY

    KRATOS_ERROR << "calling the base class function in LinearElasticModel ... illegal operation" << std::endl;

    KRATOS_CATCH(" ")
  }



  //************************************************************************************
  //************************************************************************************

  void LinearElasticModel::CalculateStressTensor(ModelDataType& rValues, MatrixType& rStressMatrix)
  {
    KRATOS_TRY

    ElasticDataType Variables;
    this->InitializeElasticData(rValues,Variables);

    VectorType StrainVector;
    ConstitutiveModelUtilities::StrainTensorToVector(rValues.StrainMatrix, StrainVector);

    this->CalculateAndAddConstitutiveTensor(Variables);

    VectorType StressVector;
    this->CalculateAndAddStressTensor(Variables,StrainVector,StressVector);

    rStressMatrix = ConstitutiveModelUtilities::VectorToSymmetricTensor(StressVector,rStressMatrix);

    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************

  void LinearElasticModel::CalculateAndAddStressTensor(ElasticDataType& rVariables, VectorType& rStrainVector, VectorType& rStressVector)
  {
    KRATOS_TRY

    noalias(rStressVector) = prod(rVariables.ConstitutiveTensor,rStrainVector);

    rVariables.State().Set(ConstitutiveModelData::STRESS_COMPUTED);

    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************

  void LinearElasticModel::CalculateIsochoricStressTensor(ModelDataType& rValues, MatrixType& rStressMatrix)
  {
    KRATOS_TRY

    ElasticDataType Variables;
    this->InitializeElasticData(rValues,Variables);

    VectorType StrainVector;
    ConstitutiveModelUtilities::StrainTensorToVector(rValues.StrainMatrix, StrainVector);

    this->CalculateAndAddConstitutiveTensor(Variables);

    VectorType StressVector;
    this->CalculateAndAddIsochoricStressTensor(Variables,StrainVector,StressVector);

    rStressMatrix = ConstitutiveModelUtilities::VectorToSymmetricTensor(StressVector,rStressMatrix);

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************

  void LinearElasticModel::CalculateAndAddIsochoricStressTensor(ElasticDataType& rVariables, VectorType& rStrainVector, VectorType& rStressVector)
  {
    KRATOS_TRY

    //total stress
    noalias(rStressVector) = prod(rVariables.ConstitutiveTensor,rStrainVector);

    //deviatoric stress
    double MeanStress = (1.0/3.0) * (rStressVector[0]+rStressVector[1]+rStressVector[2]);
    rStressVector[0] -= MeanStress;
    rStressVector[1] -= MeanStress;
    rStressVector[2] -= MeanStress;

    rVariables.State().Set(ConstitutiveModelData::STRESS_COMPUTED);

    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************

  void LinearElasticModel::CalculateVolumetricStressTensor(ModelDataType& rValues, MatrixType& rStressMatrix)
  {
    KRATOS_TRY

    this->CalculateStressTensor( rValues, rStressMatrix);
    double MeanPressure = 0;
    for (unsigned int i = 0; i < 3; i++)
       MeanPressure += rStressMatrix(i,i)/3.0;
    rStressMatrix = IdentityMatrix(3);
    rStressMatrix *= MeanPressure;


    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************

  void LinearElasticModel::CalculateAndAddVolumetricStressTensor(ElasticDataType& rVariables, VectorType& rStrainVector, VectorType& rStressVector)
  {
    KRATOS_TRY

    KRATOS_ERROR << "calling the base class function in LinearElasticModel ... illegal operation" << std::endl;

    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************

  void LinearElasticModel::CalculateConstitutiveTensor(ModelDataType& rValues, Matrix& rConstitutiveMatrix)
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

  void LinearElasticModel::CalculateAndAddConstitutiveTensor(ElasticDataType& rVariables, Matrix& rConstitutiveMatrix)
  {
    KRATOS_TRY

    this->CalculateAndAddConstitutiveTensor(rVariables);

    rConstitutiveMatrix = ConstitutiveModelUtilities::ConstitutiveTensorToMatrix(rVariables.ConstitutiveTensor,rConstitutiveMatrix);

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************

  void LinearElasticModel::CalculateAndAddConstitutiveTensor(ElasticDataType& rVariables)
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

  void LinearElasticModel::CalculateStressAndConstitutiveTensors(ModelDataType& rValues, MatrixType& rStressMatrix, Matrix& rConstitutiveMatrix)
  {
    KRATOS_TRY

    ElasticDataType Variables;
    this->InitializeElasticData(rValues,Variables);

    VectorType StrainVector;
    ConstitutiveModelUtilities::StrainTensorToVector(rValues.StrainMatrix, StrainVector);

    //Set constitutive matrix to zero before adding
    rConstitutiveMatrix.clear();

    this->CalculateAndAddConstitutiveTensor(Variables, rConstitutiveMatrix);

    VectorType StressVector;
    this->CalculateAndAddStressTensor(Variables,StrainVector,StressVector);

    rStressMatrix = ConstitutiveModelUtilities::VectorToSymmetricTensor(StressVector,rStressMatrix);


    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************

  void LinearElasticModel::CalculateIsochoricStressAndConstitutiveTensors(ModelDataType& rValues, MatrixType& rStressMatrix, Matrix& rConstitutiveMatrix)
  {
    KRATOS_TRY

    KRATOS_ERROR << "calling the base class function in LinearElasticModel ... illegal operation" << std::endl;

    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************

  void LinearElasticModel::CalculateVolumetricStressAndConstitutiveTensors(ModelDataType& rValues, MatrixType& rStressMatrix, Matrix& rConstitutiveMatrix)
  {
    KRATOS_TRY

    KRATOS_ERROR << "calling the base class function in LinearElasticModel ... illegal operation" << std::endl;

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************

  void LinearElasticModel::CalculateIsochoricConstitutiveTensor(ModelDataType& rValues, Matrix& rConstitutiveMatrix)
  {
    KRATOS_TRY

    ElasticDataType Variables;
    this->InitializeElasticData(rValues,Variables);
    this->CalculateAndAddIsochoricConstitutiveTensor(Variables,rConstitutiveMatrix);

    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************

  void LinearElasticModel::CalculateAndAddIsochoricConstitutiveTensor(ElasticDataType& rVariables, Matrix& rConstitutiveMatrix)
  {
    KRATOS_TRY

    this->CalculateAndAddIsochoricConstitutiveTensor(rVariables);

    rConstitutiveMatrix = ConstitutiveModelUtilities::ConstitutiveTensorToMatrix(rVariables.ConstitutiveTensor, rConstitutiveMatrix);

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************

  void LinearElasticModel::CalculateAndAddIsochoricConstitutiveTensor(ElasticDataType& rVariables)
  {
    KRATOS_TRY

    KRATOS_ERROR << "calling the base class function in LinearElasticModel ... illegal operation" << std::endl;

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************

  void LinearElasticModel::CalculateVolumetricConstitutiveTensor(ModelDataType& rValues, Matrix& rConstitutiveMatrix)
  {
    KRATOS_TRY

    ElasticDataType Variables;
    this->InitializeElasticData(rValues,Variables);
    this->CalculateAndAddVolumetricConstitutiveTensor(Variables,rConstitutiveMatrix);

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************

  void LinearElasticModel::CalculateAndAddVolumetricConstitutiveTensor(ElasticDataType& rVariables, Matrix& rConstitutiveMatrix)
  {
    KRATOS_TRY

    this->CalculateAndAddVolumetricConstitutiveTensor(rVariables);

    rConstitutiveMatrix = ConstitutiveModelUtilities::ConstitutiveTensorToMatrix(rVariables.ConstitutiveTensor, rConstitutiveMatrix);

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************

  void LinearElasticModel::CalculateAndAddVolumetricConstitutiveTensor(ElasticDataType& rVariables)
  {
    KRATOS_TRY

    KRATOS_ERROR << "calling the base class function in LinearElasticModel ... illegal operation" << std::endl;

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************

  int LinearElasticModel::Check(const Properties& rMaterialProperties, const ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    if(YOUNG_MODULUS.Key() == 0 || rMaterialProperties[YOUNG_MODULUS] <= 0.00)
      KRATOS_ERROR << "YOUNG_MODULUS has Key zero or invalid value" << std::endl;

    if(POISSON_RATIO.Key() == 0){
      KRATOS_ERROR << "POISSON_RATIO has Key zero invalid value" << std::endl;
    }
    else{
      const double& nu = rMaterialProperties[POISSON_RATIO];
      if( (nu > 0.499 && nu < 0.501) || (nu < -0.999 && nu > -1.01) )
	KRATOS_ERROR << "POISSON_RATIO has an invalid value" << std::endl;
    }

    return 0;


    KRATOS_CATCH(" ")
  }



} // Namespace Kratos
