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
#include "custom_models/elasticity_models/hyper_elastic_model.hpp"


namespace Kratos
{

  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

  HyperElasticModel::HyperElasticModel()
    : ConstitutiveModel()
    , msIdentityMatrix( identity_matrix<double>(3) )
  {
    KRATOS_TRY

    noalias(this->mHistoryVector) = ZeroVector(6);

    this->mHistoryVector[0] = 1.0;
    this->mHistoryVector[1] = 1.0;
    this->mHistoryVector[2] = 1.0;

    KRATOS_CATCH(" ")
  }

  //******************************COPY CONSTRUCTOR**************************************
  //************************************************************************************

  HyperElasticModel::HyperElasticModel(const HyperElasticModel& rOther)
    : ConstitutiveModel(rOther)
    , msIdentityMatrix(rOther.msIdentityMatrix)
  {
  }

  //********************************CLONE***********************************************
  //************************************************************************************

  ConstitutiveModel::Pointer HyperElasticModel::Clone() const
  {
    return Kratos::make_shared<HyperElasticModel>(*this);
  }

  //********************************ASSIGNMENT******************************************
  //************************************************************************************
  HyperElasticModel& HyperElasticModel::operator=(HyperElasticModel const& rOther)
  {
    ConstitutiveModel::operator=(rOther);
    return *this;
  }

  //*******************************DESTRUCTOR*******************************************
  //************************************************************************************

  HyperElasticModel::~HyperElasticModel()
  {
  }


  //***********************PUBLIC OPERATIONS FROM BASE CLASS****************************
  //************************************************************************************

  void HyperElasticModel::InitializeModel(ModelDataType& rValues)
  {
    KRATOS_TRY


    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************

  void HyperElasticModel::FinalizeModel(ModelDataType& rValues)
  {
    KRATOS_TRY

    //update total strain measure
    ConstitutiveModelUtilities::SymmetricTensorToVector(rValues.StrainMatrix, this->mHistoryVector);

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************

  void HyperElasticModel::CalculateStrainEnergy(ModelDataType& rValues, double& rDensityFunction)
  {
    KRATOS_TRY

    KRATOS_ERROR << "calling the base class function in HyperElasticModel ... illegal operation" << std::endl;

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************

  void HyperElasticModel::CalculateAndAddIsochoricStrainEnergy(HyperElasticDataType& rVariables, double& rIsochoricDensityFunction)
  {
    KRATOS_TRY

    KRATOS_ERROR << "calling the base class function in HyperElasticModel ... illegal operation" << std::endl;

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************

  void HyperElasticModel::CalculateAndAddVolumetricStrainEnergy(HyperElasticDataType& rVariables, double& rVolumetricDensityFunction)
  {
    KRATOS_TRY

    KRATOS_ERROR << "calling the base class function in HyperElasticModel ... illegal operation" << std::endl;

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************

  void HyperElasticModel::CalculateStressTensor(ModelDataType& rValues, MatrixType& rStressMatrix)
  {
    KRATOS_TRY

    HyperElasticDataType Variables;
    this->CalculateStrainData(rValues,Variables);

    this->CalculateAndAddStressTensor(Variables,rStressMatrix);

    rValues.StressMatrix = rStressMatrix; //store total stress as StressMatrix

    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************

  void HyperElasticModel::CalculateAndAddStressTensor(HyperElasticDataType& rVariables, MatrixType& rStressMatrix)
  {
    KRATOS_TRY

    KRATOS_ERROR << "calling the base class function in HyperElasticModel ... illegal operation" << std::endl;

    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************

  void HyperElasticModel::CalculateIsochoricStressTensor(ModelDataType& rValues, MatrixType& rStressMatrix)
  {
    KRATOS_TRY

    HyperElasticDataType Variables;
    this->CalculateStrainData(rValues,Variables);
    this->CalculateAndAddIsochoricStressTensor(Variables,rStressMatrix);

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************

  void HyperElasticModel::CalculateAndAddIsochoricStressTensor(HyperElasticDataType& rVariables, MatrixType& rStressMatrix)
  {
    KRATOS_TRY

    KRATOS_ERROR << "calling the base class function in HyperElasticModel ... illegal operation" << std::endl;

    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************

  void HyperElasticModel::CalculateVolumetricStressTensor(ModelDataType& rValues, MatrixType& rStressMatrix)
  {
    KRATOS_TRY

    HyperElasticDataType Variables;
    this->CalculateStrainData(rValues,Variables);
    this->CalculateAndAddVolumetricStressTensor(Variables,rStressMatrix);

    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************

  void HyperElasticModel::CalculateAndAddVolumetricStressTensor(HyperElasticDataType& rVariables, MatrixType& rStressMatrix)
  {
    KRATOS_TRY

    KRATOS_ERROR << "calling the base class function in HyperElasticModel ... illegal operation" << std::endl;

    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************

  void HyperElasticModel::CalculateConstitutiveTensor(ModelDataType& rValues, Matrix& rConstitutiveMatrix)
  {
    KRATOS_TRY

    //Initialize ConstitutiveMatrix
    HyperElasticDataType Variables;
    this->CalculateStrainData(rValues,Variables);

    //Set constitutive matrix to zero before adding
    rConstitutiveMatrix.clear();

    //Calculate Constitutive Matrix
    this->CalculateAndAddConstitutiveTensor(Variables,rConstitutiveMatrix);

    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************

  void HyperElasticModel::CalculateAndAddConstitutiveTensor(HyperElasticDataType& rVariables, Matrix& rConstitutiveMatrix)
  {
    KRATOS_TRY

    //Calculate HyperElastic ConstitutiveMatrix
    const ModelDataType&  rModelData        = rVariables.GetModelData();
    const SizeType&       rVoigtSize        = rModelData.GetVoigtSize();
    const VoigtIndexType& rIndexVoigtTensor = rModelData.GetVoigtIndexTensor();


    for(SizeType i=0; i<rVoigtSize; i++)
      {
	for(SizeType j=0; j<rVoigtSize; j++)
	  {
	    rConstitutiveMatrix(i,j) = this->AddConstitutiveComponent(rVariables,rConstitutiveMatrix(i,j),
								      rIndexVoigtTensor[i][0],rIndexVoigtTensor[i][1],
								      rIndexVoigtTensor[j][0],rIndexVoigtTensor[j][1]);
	  }

      }

    rVariables.State().Set(ConstitutiveModelData::CONSTITUTIVE_MATRIX_COMPUTED,true);

    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************

  void HyperElasticModel::CalculateStressAndConstitutiveTensors(ModelDataType& rValues, MatrixType& rStressMatrix, Matrix& rConstitutiveMatrix)
  {
    KRATOS_TRY

    HyperElasticDataType Variables;
    this->CalculateStrainData(rValues,Variables);

    //Calculate Stress Matrix
    this->CalculateAndAddStressTensor(Variables,rStressMatrix);

    rValues.StressMatrix = rStressMatrix; //store total stress as StressMatrix

    //Set constitutive matrix to zero before adding
    rConstitutiveMatrix.clear();

    //Calculate Constitutive Matrix
    this->CalculateAndAddConstitutiveTensor(Variables,rConstitutiveMatrix);

    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************

  void HyperElasticModel::CalculateIsochoricStressAndConstitutiveTensors(ModelDataType& rValues, MatrixType& rStressMatrix, Matrix& rConstitutiveMatrix)
  {
    KRATOS_TRY

    HyperElasticDataType Variables;
    this->CalculateStrainData(rValues,Variables);

    //Calculate Stress Matrix
    this->CalculateAndAddIsochoricStressTensor(Variables,rStressMatrix);

    rValues.StressMatrix = rStressMatrix; //store isochoric stress as StressMatrix

    //Set constitutive matrix to zero before adding
    rConstitutiveMatrix.clear();

    //Calculate Constitutive Matrix
    this->CalculateAndAddIsochoricConstitutiveTensor(Variables,rConstitutiveMatrix);

    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************

  void HyperElasticModel::CalculateVolumetricStressAndConstitutiveTensors(ModelDataType& rValues, MatrixType& rStressMatrix, Matrix& rConstitutiveMatrix)
  {
    KRATOS_TRY

    HyperElasticDataType Variables;
    this->CalculateStrainData(rValues,Variables);

    //Calculate Stress Matrix
    this->CalculateAndAddVolumetricStressTensor(Variables,rStressMatrix);

    rValues.StressMatrix = rStressMatrix; //store volumetric stress as StressMatrix

    //Set constitutive matrix to zero before adding
    rConstitutiveMatrix.clear();

    //Calculate Constitutive Matrix
    this->CalculateAndAddVolumetricConstitutiveTensor(Variables,rConstitutiveMatrix);

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************

  void HyperElasticModel::CalculateIsochoricConstitutiveTensor(ModelDataType& rValues, Matrix& rConstitutiveMatrix)
  {
    KRATOS_TRY

    HyperElasticDataType Variables;
    this->CalculateStrainData(rValues,Variables);

    //Calculate Isochoric Constitutive Matrix
    this->CalculateAndAddIsochoricConstitutiveTensor(Variables,rConstitutiveMatrix);

    KRATOS_CATCH(" ")
  }


  //***********************PROTECTED OPERATIONS FROM BASE CLASS*************************
  //************************************************************************************

  void HyperElasticModel::CalculateStrainData(ModelDataType& rValues, HyperElasticDataType& rVariables)
  {
    KRATOS_TRY

    KRATOS_ERROR << "calling the base class function in HyperElasticModel ... illegal operation" << std::endl;

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************

  void HyperElasticModel::CalculateAndAddIsochoricConstitutiveTensor(HyperElasticDataType& rVariables, Matrix& rConstitutiveMatrix)
  {
    KRATOS_TRY

    //Calculate HyperElastic ConstitutiveMatrix
    const ModelDataType&  rModelData        = rVariables.GetModelData();
    const SizeType&       rVoigtSize        = rModelData.GetVoigtSize();
    const VoigtIndexType& rIndexVoigtTensor = rModelData.GetVoigtIndexTensor();

    for(SizeType i=0; i<rVoigtSize; i++)
      {
	for(SizeType j=0; j<rVoigtSize; j++)
	  {

	    rConstitutiveMatrix(i,j) = this->AddIsochoricConstitutiveComponent(rVariables,rConstitutiveMatrix(i,j),
									       rIndexVoigtTensor[i][0],rIndexVoigtTensor[i][1],
									       rIndexVoigtTensor[j][0],rIndexVoigtTensor[j][1]);
	  }

      }

    rVariables.State().Set(ConstitutiveModelData::CONSTITUTIVE_MATRIX_COMPUTED,true);

    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************

  void HyperElasticModel::CalculateVolumetricConstitutiveTensor(ModelDataType& rValues, Matrix& rConstitutiveMatrix)
  {
    KRATOS_TRY

    HyperElasticDataType Variables;
    this->CalculateStrainData(rValues,Variables);

    //Calculate Volumetric Constitutive Matrix
    this->CalculateAndAddVolumetricConstitutiveTensor(Variables,rConstitutiveMatrix);


    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************

  void HyperElasticModel::CalculateAndAddVolumetricConstitutiveTensor(HyperElasticDataType& rVariables, Matrix& rConstitutiveMatrix)
  {
    KRATOS_TRY

    //Calculate HyperElastic ConstitutiveMatrix
    const ModelDataType&  rModelData        = rVariables.GetModelData();
    const SizeType&       rVoigtSize        = rModelData.GetVoigtSize();
    const VoigtIndexType& rIndexVoigtTensor = rModelData.GetVoigtIndexTensor();

    for(SizeType i=0; i<rVoigtSize; i++)
      {
	for(SizeType j=0; j<rVoigtSize; j++)
	  {

	    rConstitutiveMatrix(i,j) = this->AddVolumetricConstitutiveComponent(rVariables,rConstitutiveMatrix(i,j),
										rIndexVoigtTensor[i][0],rIndexVoigtTensor[i][1],
										rIndexVoigtTensor[j][0],rIndexVoigtTensor[j][1]);
	  }

      }

    rVariables.State().Set(ConstitutiveModelData::CONSTITUTIVE_MATRIX_COMPUTED,true);

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************

  double& HyperElasticModel::AddConstitutiveComponent(HyperElasticDataType& rVariables, double &rCabcd,
						      const unsigned int& a, const unsigned int& b,
						      const unsigned int& c, const unsigned int& d)
  {
    KRATOS_TRY

    KRATOS_ERROR << "calling the base class function in HyperElasticModel ... illegal operation" << std::endl;

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************

  double& HyperElasticModel::AddIsochoricConstitutiveComponent(HyperElasticDataType& rVariables, double &rCabcd,
							       const unsigned int& a, const unsigned int& b,
							       const unsigned int& c, const unsigned int& d)
  {
    KRATOS_TRY

    KRATOS_ERROR << "calling the base class function in HyperElasticModel ... illegal operation" << std::endl;

    KRATOS_CATCH(" ")
  }



  //************************************************************************************
  //************************************************************************************

  double& HyperElasticModel::AddVolumetricConstitutiveComponent(HyperElasticDataType& rVariables, double &rCabcd,
								const unsigned int& a, const unsigned int& b,
								const unsigned int& c, const unsigned int& d)
  {
    KRATOS_TRY

    KRATOS_ERROR << "calling the base class function in HyperElasticModel ... illegal operation" << std::endl;

    KRATOS_CATCH(" ")
  }

  //***********************PROTECTED OPERATIONS FROM BASE CLASS*************************
  //************************************************************************************


  void HyperElasticModel::CalculateStrainInvariants(const MatrixType& rStrainMatrix, double& rI1, double& rI2, double& rI3)
  {
    KRATOS_TRY

    rI1 = rStrainMatrix(0,0) + rStrainMatrix(1,1) + rStrainMatrix(2,2);

    rI2 = (  rStrainMatrix(0,0)*rStrainMatrix(1,1)
	     + rStrainMatrix(1,1)*rStrainMatrix(2,2)
	     + rStrainMatrix(0,0)*rStrainMatrix(2,2)
	     - rStrainMatrix(0,1)*rStrainMatrix(1,0)
	     - rStrainMatrix(1,2)*rStrainMatrix(2,1)
	     - rStrainMatrix(0,2)*rStrainMatrix(2,0) );

    rI3 = (  rStrainMatrix(0,0)*rStrainMatrix(1,1)*rStrainMatrix(2,2)
	     + rStrainMatrix(0,1)*rStrainMatrix(1,2)*rStrainMatrix(2,0)
	     + rStrainMatrix(1,0)*rStrainMatrix(2,1)*rStrainMatrix(0,2)
	     - rStrainMatrix(2,0)*rStrainMatrix(1,1)*rStrainMatrix(0,2)
	     - rStrainMatrix(1,0)*rStrainMatrix(0,1)*rStrainMatrix(2,2)
	     - rStrainMatrix(0,0)*rStrainMatrix(2,1)*rStrainMatrix(1,2) );

    //std::cout<<" I1: "<<rI1<<" I2: "<<rI2<<" I3: "<<rI3<<std::endl;

    KRATOS_CATCH(" ")
  }



  //************************************************************************************
  //************************************************************************************

  void HyperElasticModel::CalculateInvariants(HyperElasticDataType& rVariables)
  {
    KRATOS_TRY


    //invariants
    this->CalculateStrainInvariants( rVariables.Strain.Matrix, rVariables.Strain.Invariants.I1, rVariables.Strain.Invariants.I2, rVariables.Strain.Invariants.I3 );

    //jacobian
    rVariables.Strain.Invariants.J    = rVariables.GetModelData().GetTotalDeformationDet();
    rVariables.Strain.Invariants.J_13 = pow(rVariables.Strain.Invariants.J,(-1.0/3.0));

    //std::cout<<" Strain.Invariants [I1:"<<rVariables.Strain.Invariants.I1<<" I2:"<<rVariables.Strain.Invariants.I2<<" I3:"<<rVariables.Strain.Invariants.I3<<"] J:"<<rVariables.Strain.Invariants.J<<std::endl;

    rVariables.Strain.Invariants.I3 = rVariables.Strain.Invariants.J * rVariables.Strain.Invariants.J; //for volumetric consistency


    KRATOS_CATCH(" ")
  }


  void HyperElasticModel::CalculateScalingFactors(HyperElasticDataType& rVariables)
  {
    KRATOS_TRY

    KRATOS_ERROR << "calling the base class function in HyperElasticModel ... illegal operation" << std::endl;

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************


  //************// right cauchy green: C

  HyperElasticModel::MatrixType& HyperElasticModel::GetJRightCauchyGreenDerivative(const StrainData& rStrain, MatrixType& rDerivative) //dJ/dC
  {
    KRATOS_TRY

    noalias(rDerivative) = rStrain.InverseMatrix;
    rDerivative *= rStrain.Invariants.J * 0.5;

    return rDerivative;

    KRATOS_CATCH(" ")
  }


  double& HyperElasticModel::GetJRightCauchyGreen1stDerivative(const StrainData& rStrain,
							       double& rDerivative,
							       const double& a,
							       const double& b) ///dJ/dC
  {
    KRATOS_TRY

    rDerivative = 0.5 * rStrain.Invariants.J * rStrain.InverseMatrix(a,b);

    return rDerivative;

    KRATOS_CATCH(" ")
  }

  double& HyperElasticModel::GetJRightCauchyGreenSquare1stDerivative(const StrainData& rStrain,
								     double& rDerivative,
								     const double& a,
								     const double& b,
								     const double& c,
								     const double& d) //dJ/dC * dJ/dC
  {
    KRATOS_TRY

    rDerivative  = 0.5 * rStrain.Invariants.J * rStrain.InverseMatrix(a,b);
    rDerivative *= (0.5 * rStrain.Invariants.J * rStrain.InverseMatrix(c,d));

    return rDerivative;

    KRATOS_CATCH(" ")
  }


  double& HyperElasticModel::GetJRightCauchyGreen2ndDerivative(const StrainData& rStrain,
							       double& rDerivative,
							       const double& a,
							       const double& b,
							       const double& c,
							       const double& d) //ddJ/dCdC
  {
    KRATOS_TRY

    rDerivative  = (-1) * ConstitutiveModelUtilities::CalculateFourthOrderTensor(rStrain.InverseMatrix,rDerivative,a,b,c,d);
    rDerivative += 0.5 * rStrain.InverseMatrix(a,b)*rStrain.InverseMatrix(c,d);
    rDerivative *= 0.5 * rStrain.Invariants.J;

    return rDerivative;

    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************


  //************// left cauchy green : b

  HyperElasticModel::MatrixType& HyperElasticModel::GetJLeftCauchyGreenDerivative(const StrainData& rStrain, MatrixType& rDerivative) //dJ/db
  {
    KRATOS_TRY

    noalias(rDerivative)  = this->msIdentityMatrix;
    rDerivative *= rStrain.Invariants.J * 0.5;

    return rDerivative;

    KRATOS_CATCH(" ")
  }

  double& HyperElasticModel::GetJLeftCauchyGreen1stDerivative(const StrainData& rStrain,
							      double& rDerivative,
							      const double& a,
							      const double& b) //dJ/db
  {
    KRATOS_TRY

    rDerivative  = 0.5 * rStrain.Invariants.J * this->msIdentityMatrix(a,b);

    return rDerivative;

    KRATOS_CATCH(" ")
  }


  double& HyperElasticModel::GetJLeftCauchyGreenSquare1stDerivative(const StrainData& rStrain,
								    double& rDerivative,
								    const double& a,
								    const double& b,
								    const double& c,
								    const double& d) //dJ/db * dJ/db
  {
    KRATOS_TRY

    rDerivative  = 0.5 * rStrain.Invariants.J * this->msIdentityMatrix(a,b);
    rDerivative *= (0.5 * rStrain.Invariants.J * this->msIdentityMatrix(c,d));

    return rDerivative;

    KRATOS_CATCH(" ")
  }

  double& HyperElasticModel::GetJLeftCauchyGreen2ndDerivative(const StrainData& rStrain,
							      double& rDerivative,
							      const double& a,
							      const double& b,
							      const double& c,
							      const double& d) //ddJ/dbdb
  {
    KRATOS_TRY

    rDerivative  = (-1.0) * ConstitutiveModelUtilities::CalculateFourthOrderUnitTensor(this->msIdentityMatrix,rDerivative,a,b,c,d);
    rDerivative += 0.5 * this->msIdentityMatrix(a,b)*this->msIdentityMatrix(c,d);
    rDerivative *= 0.5 * rStrain.Invariants.J;

    return rDerivative;

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************

  //isochoric volumetric slit

  double& HyperElasticModel::GetVolumetricFunction1stJDerivative(HyperElasticDataType& rVariables, double& rDerivative) //dU/dJ
  {
    KRATOS_TRY

    KRATOS_ERROR << "calling the base class function in HyperElasticModel ... illegal operation" << std::endl;

    return rDerivative;

    KRATOS_CATCH(" ")
  }

  double& HyperElasticModel::GetVolumetricFunction2ndJDerivative(HyperElasticDataType& rVariables, double& rDerivative) //ddU/dJdJ
  {
    KRATOS_TRY

    KRATOS_ERROR << "calling the base class function in HyperElasticModel ... illegal operation" << std::endl;

    return rDerivative;

    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************

  int HyperElasticModel::Check(const Properties& rMaterialProperties, const ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    return 0;

    KRATOS_CATCH(" ")
  }



} // Namespace Kratos
