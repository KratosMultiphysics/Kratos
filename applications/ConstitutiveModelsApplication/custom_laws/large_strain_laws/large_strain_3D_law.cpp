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
#include "custom_laws/large_strain_laws/large_strain_3D_law.hpp"

namespace Kratos
{

  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

  LargeStrain3DLaw::LargeStrain3DLaw()
    : Constitutive3DLaw()
  {
    KRATOS_TRY

    //member variables initialization
    mTotalDeformationDet = 1.0;

    MatrixType Identity = IdentityMatrix(3);
    noalias(mInverseTotalDeformationMatrix) = Identity;


    KRATOS_CATCH(" ")
  }

  //******************************CONSTRUCTOR WITH THE MODEL****************************
  //************************************************************************************

  LargeStrain3DLaw::LargeStrain3DLaw(ModelTypePointer pModel)
    : Constitutive3DLaw()
  {
    KRATOS_TRY

    //model
    mpModel = pModel->Clone();

    //member variables initialization
    mTotalDeformationDet = 1.0;

    MatrixType Identity = IdentityMatrix(3);
    noalias(mInverseTotalDeformationMatrix) = Identity;

    KRATOS_CATCH(" ")
  }

  //******************************COPY CONSTRUCTOR**************************************
  //************************************************************************************

  LargeStrain3DLaw::LargeStrain3DLaw(const LargeStrain3DLaw& rOther)
    : Constitutive3DLaw(rOther)
    ,mTotalDeformationDet(rOther.mTotalDeformationDet)
    ,mInverseTotalDeformationMatrix(rOther.mInverseTotalDeformationMatrix)
  {
    mpModel = rOther.mpModel->Clone();
  }

  //*******************************ASSIGMENT OPERATOR***********************************
  //************************************************************************************

  LargeStrain3DLaw& LargeStrain3DLaw::operator=(const LargeStrain3DLaw& rOther)
  {
    Constitutive3DLaw::operator=(rOther);
    mpModel = rOther.mpModel->Clone();
    mTotalDeformationDet = rOther.mTotalDeformationDet;
    mInverseTotalDeformationMatrix = rOther.mInverseTotalDeformationMatrix;
    return *this;
  }

  //********************************CLONE***********************************************
  //************************************************************************************

  ConstitutiveLaw::Pointer LargeStrain3DLaw::Clone() const
  {
    return Kratos::make_shared<LargeStrain3DLaw>(*this);
  }


  //*******************************DESTRUCTOR*******************************************
  //************************************************************************************

  LargeStrain3DLaw::~LargeStrain3DLaw()
  {
  }

  //*******************************OPERATIONS FROM BASE CLASS***************************
  //************************************************************************************

  //***********************HAS : DOUBLE - VECTOR - MATRIX*******************************
  //************************************************************************************

  bool LargeStrain3DLaw::Has( const Variable<double>& rThisVariable )
  {
    KRATOS_TRY

    if(rThisVariable == DETERMINANT_F)
      return true;

    return mpModel->Has(rThisVariable);


    KRATOS_CATCH(" ")
  }


  //***********************SET VALUE: DOUBLE - VECTOR - MATRIX**************************
  //************************************************************************************

  void LargeStrain3DLaw::SetValue( const Variable<double>& rThisVariable, const double& rValue,
				    const ProcessInfo& rCurrentProcessInfo )
  {
    KRATOS_TRY

    mpModel->SetValue(rThisVariable,rValue, rCurrentProcessInfo);

    if(rThisVariable == DETERMINANT_F){
      mTotalDeformationDet = rValue;
    }

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************

  void LargeStrain3DLaw::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
				    const ProcessInfo& rCurrentProcessInfo )
  {
    KRATOS_TRY

    mpModel->SetValue(rThisVariable,rValue, rCurrentProcessInfo);

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************

  void LargeStrain3DLaw::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
				   const ProcessInfo& rCurrentProcessInfo )
  {
    KRATOS_TRY

    mpModel->SetValue(rThisVariable,rValue, rCurrentProcessInfo);

    KRATOS_CATCH(" ")
  }


  //***********************GET VALUE: DOUBLE - VECTOR - MATRIX**************************
  //************************************************************************************

  double& LargeStrain3DLaw::GetValue( const Variable<double>& rThisVariable, double& rValue )
  {
    KRATOS_TRY

    rValue = mpModel->GetValue(rThisVariable,rValue);

    if(rThisVariable == DETERMINANT_F){
      rValue = mTotalDeformationDet;
    } else {
       rValue = mpModel->GetValue(rThisVariable,rValue);
    }

    return rValue;

    KRATOS_CATCH(" ")
  }

  //************* STARTING - ENDING  METHODS
  //************************************************************************************
  //************************************************************************************


  //************* COMPUTING  METHODS
  //************************************************************************************
  //************************************************************************************
  void LargeStrain3DLaw::InitializeMaterial( const Properties& rProperties,
					     const GeometryType& rElementGeometry,
					     const Vector& rShapeFunctionsValues )
  {
    KRATOS_TRY

    ConstitutiveLaw::InitializeMaterial(rProperties,rElementGeometry,rShapeFunctionsValues);

    //member variables initialization
    mTotalDeformationDet = 1.0;

    MatrixType Identity = IdentityMatrix(3);
    noalias(mInverseTotalDeformationMatrix) = Identity;

    mpModel->InitializeMaterial(rProperties);

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************

  void LargeStrain3DLaw::InitializeModelData(Parameters& rValues,ModelDataType& rModelValues)
  {
    KRATOS_TRY

    if(rValues.GetMaterialProperties().Has(PROPERTIES_LAYOUT))
    {
      PropertiesLayout::Pointer pPropertiesLayout = rValues.GetMaterialProperties()[PROPERTIES_LAYOUT].Clone();
      pPropertiesLayout->Configure(rValues.GetMaterialProperties(),rValues.GetElementGeometry(),rValues.GetShapeFunctionsValues());
      rModelValues.SetPropertiesLayout(pPropertiesLayout);
    }

    rModelValues.SetOptions(rValues.GetOptions());
    rModelValues.SetProperties(rValues.GetMaterialProperties());
    rModelValues.SetProcessInfo(rValues.GetProcessInfo());
    rModelValues.SetVoigtSize(this->GetStrainSize());
    rModelValues.SetVoigtIndexTensor(this->GetVoigtIndexTensor());

    LawDataType& rVariables = rModelValues.rConstitutiveLawData();

    //a.- Calculate incremental deformation gradient determinant
    rVariables.TotalDeformationDet = rValues.GetDeterminantF();
    rVariables.DeltaDeformationDet = rVariables.TotalDeformationDet/mTotalDeformationDet; //determinant incremental F

    //b.- Calculate incremental deformation gradient
    const MatrixType& rTotalDeformationMatrix = rValues.GetDeformationGradientF();

    rVariables.TotalDeformationMatrix = ConstitutiveModelUtilities::DeformationGradientTo3D(rTotalDeformationMatrix, rVariables.TotalDeformationMatrix);
    noalias(rVariables.DeltaDeformationMatrix) = prod(rVariables.TotalDeformationMatrix, mInverseTotalDeformationMatrix); //incremental F

    if( rValues.GetOptions().Is(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE) )
      rModelValues.State.Set(ConstitutiveModelData::UPDATE_INTERNAL_VARIABLES);

    //initialize model
    mpModel->InitializeModel(rModelValues);


    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************

  void LargeStrain3DLaw::FinalizeModelData(Parameters& rValues,ModelDataType& rModelValues)
  {
    KRATOS_TRY

    //Finalize Material response
    if(rValues.GetOptions().Is(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE)){

      const Matrix& rDeformationMatrix = rValues.GetDeformationGradientF();
      const double& rDeformationDet    = rValues.GetDeterminantF();

      //update total deformation gradient
      MatrixType TotalDeformationMatrix;
      TotalDeformationMatrix = ConstitutiveModelUtilities::DeformationGradientTo3D(rDeformationMatrix,TotalDeformationMatrix);
      ConstitutiveModelUtilities::InvertMatrix3( TotalDeformationMatrix, mInverseTotalDeformationMatrix, mTotalDeformationDet);
      mTotalDeformationDet = rDeformationDet; //special treatment of the determinant

      //finalize model (update total strain measure)
      mpModel->FinalizeModel(rModelValues);

    }

    KRATOS_CATCH(" ")
  }


  //*****************************MATERIAL RESPONSES*************************************
  //************************************************************************************

  void LargeStrain3DLaw::CalculateMaterialResponsePK2(Parameters& rValues)
  {
    KRATOS_TRY

    ModelDataType ModelValues;

    this->CalculateMaterialResponsePK2(rValues,ModelValues);

    KRATOS_CATCH(" ")
  }


  void LargeStrain3DLaw::CalculateMaterialResponsePK2(Parameters& rValues, ModelDataType& rModelValues)
  {
    KRATOS_TRY

    //0.- Check if the constitutive parameters are passed correctly to the law calculation
    //CheckParameters(rValues);

    const Flags& rOptions = rValues.GetOptions();

    //1.- Initialize hyperelastic model parameters
    ModelDataType rModelValues;

    LawDataType& rVariables = rModelValues.rConstitutiveLawData();
    rVariables.StressMeasure = ConstitutiveModelData::StressMeasureType::StressMeasure_PK2; //required stress measure

    this->InitializeModelData(rValues, rModelValues);

    //2.-Calculate domain variables (Temperature, Pressure, Size) and calculate material parameters
    this->CalculateDomainVariables(rValues, rModelValues);

    ConstitutiveModelData::CalculateMaterialParameters(rModelValues);

    //3.-Calculate Total PK2 stress and  Constitutive Matrix related to Total PK2 stress
    if(rOptions.Is(ConstitutiveLaw::COMPUTE_STRESS) && rOptions.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)){

      Vector& rStressVector       = rValues.GetStressVector();
      Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix();

      this->CalculateStressVectorAndConstitutiveMatrix(rModelValues, rStressVector, rConstitutiveMatrix);

    }
    else{

      //4.-Calculate Total PK2 stress
      if(rOptions.Is(ConstitutiveLaw::COMPUTE_STRESS)){

	Vector& rStressVector       = rValues.GetStressVector();
	this->CalculateStressVector(rModelValues, rStressVector);

      }

      //5.-Calculate Constitutive Matrix related to Total PK2 stress
      if(rOptions.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)){

      	Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix();
	this->CalculateConstitutiveMatrix(rModelValues, rConstitutiveMatrix);

      }

    }

    if(rOptions.Is(ConstitutiveLaw::COMPUTE_STRAIN_ENERGY))
      {

      }


    //6.- Finalize hyperelastic model parameters
    this->FinalizeModelData(rValues,rModelValues);

    // std::cout<<" StrainVector "<<rValues.GetStrainVector()<<std::endl;
    // std::cout<<" StressVector "<<rValues.GetStressVector()<<std::endl;
    // std::cout<<" ConstitutiveMatrix "<<rValues.GetConstitutiveMatrix()<<std::endl;

    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************

  void LargeStrain3DLaw::CalculateMaterialResponseKirchhoff(Parameters& rValues)
  {
    KRATOS_TRY

    ModelDataType ModelValues;

    this->CalculateMaterialResponseKirchhoff(rValues,ModelValues);

    KRATOS_CATCH(" ")
  }

  void LargeStrain3DLaw::CalculateMaterialResponseKirchhoff(Parameters& rValues, ModelDataType& rModelValues)
  {
    KRATOS_TRY

    //0.- Check if the constitutive parameters are passed correctly to the law calculation
    //CheckParameters(rValues);

    const Flags& rOptions = rValues.GetOptions();

    //1.- Initialize hyperelastic model parameters
    LawDataType& rVariables = rModelValues.rConstitutiveLawData();
    rVariables.StressMeasure = ConstitutiveModelData::StressMeasureType::StressMeasure_Kirchhoff; //set required stress measure

    this->InitializeModelData(rValues, rModelValues);

    //2.-Calculate domain variables (Temperature, Pressure, Size) and calculate material parameters
    this->CalculateDomainVariables(rValues, rModelValues);

    ConstitutiveModelData::CalculateMaterialParameters(rModelValues);

    //3.-Calculate Total kirchhoff stress and  Constitutive Matrix related to Total Kirchhoff stress

    if(rOptions.Is(ConstitutiveLaw::COMPUTE_STRESS) && rOptions.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)){

      Vector& rStressVector       = rValues.GetStressVector();
      Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix();

      this->CalculateStressVectorAndConstitutiveMatrix(rModelValues, rStressVector, rConstitutiveMatrix);

    }
    else{

      //4.-Calculate Total Kirchhoff stress

      if(rOptions.Is(ConstitutiveLaw::COMPUTE_STRESS)){

	Vector& rStressVector       = rValues.GetStressVector();
	this->CalculateStressVector(rModelValues, rStressVector);

      }

      //5.-Calculate Constitutive Matrix related to Total Kirchhoff stress

      if(rOptions.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)){

      	Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix();
	this->CalculateConstitutiveMatrix(rModelValues, rConstitutiveMatrix);

      }

    }

    if(rOptions.Is(ConstitutiveLaw::COMPUTE_STRAIN_ENERGY))
      {

      }


    //6.- Finalize hyperelastic model parameters
    this->FinalizeModelData(rValues,rModelValues);


    // std::cout<<" StrainVector "<<rValues.GetStrainVector()<<std::endl;
    // std::cout<<" StressVector "<<rValues.GetStressVector()<<std::endl;
    // std::cout<<" ConstitutiveMatrix "<<rValues.GetConstitutiveMatrix()<<std::endl;


    KRATOS_CATCH(" ")

  }


  //*******************************COMPUTE STRESS VECTOR********************************
  //************************************************************************************

  void LargeStrain3DLaw::CalculateStressVector(ModelDataType& rModelValues, Vector& rStressVector)
  {
    KRATOS_TRY

    MatrixType StressMatrix;
    StressMatrix.clear();

    if(rModelValues.GetOptions().Is(ConstitutiveLaw::ISOCHORIC_TENSOR_ONLY)){
      mpModel->CalculateIsochoricStressTensor(rModelValues, StressMatrix);
    }
    else if(rModelValues.GetOptions().Is(ConstitutiveLaw::VOLUMETRIC_TENSOR_ONLY)){
      mpModel->CalculateVolumetricStressTensor(rModelValues, StressMatrix);
    }
    else{
      mpModel->CalculateStressTensor(rModelValues, StressMatrix);
    }

    rStressVector = ConstitutiveModelUtilities::StressTensorToVector(StressMatrix, rStressVector);

    KRATOS_CATCH(" ")
  }

  //***********************COMPUTE ALGORITHMIC CONSTITUTIVE MATRIX**********************
  //************************************************************************************

  void LargeStrain3DLaw::CalculateConstitutiveMatrix(ModelDataType& rModelValues, Matrix& rConstitutiveMatrix)
  {
    KRATOS_TRY

    //Calculate ConstitutiveMatrix
    if(rModelValues.GetOptions().Is(ConstitutiveLaw::ISOCHORIC_TENSOR_ONLY)){

      mpModel->CalculateIsochoricConstitutiveTensor(rModelValues, rConstitutiveMatrix);
    }
    else if(rModelValues.GetOptions().Is(ConstitutiveLaw::VOLUMETRIC_TENSOR_ONLY)){

      mpModel->CalculateVolumetricConstitutiveTensor(rModelValues, rConstitutiveMatrix);
    }
    else{

      mpModel->CalculateConstitutiveTensor(rModelValues, rConstitutiveMatrix);
    }

    KRATOS_CATCH(" ")
  }

  //******************COMPUTE STRESS AND ALGORITHMIC CONSTITUTIVE MATRIX****************
  //************************************************************************************

  void LargeStrain3DLaw::CalculateStressVectorAndConstitutiveMatrix(ModelDataType& rModelValues, Vector& rStressVector, Matrix& rConstitutiveMatrix)
  {
    KRATOS_TRY

    MatrixType StressMatrix;
    StressMatrix.clear();

    //Calculate Stress and ConstitutiveMatrix
    if(rModelValues.GetOptions().Is(ConstitutiveLaw::ISOCHORIC_TENSOR_ONLY)){

      mpModel->CalculateIsochoricStressAndConstitutiveTensors(rModelValues, StressMatrix, rConstitutiveMatrix);
    }
    else if(rModelValues.GetOptions().Is(ConstitutiveLaw::VOLUMETRIC_TENSOR_ONLY)){

      mpModel->CalculateVolumetricStressAndConstitutiveTensors(rModelValues, StressMatrix, rConstitutiveMatrix);
    }
    else{

      mpModel->CalculateStressAndConstitutiveTensors(rModelValues, StressMatrix, rConstitutiveMatrix);
    }

    rStressVector = ConstitutiveModelUtilities::StressTensorToVector(StressMatrix, rStressVector);


    KRATOS_CATCH(" ")
  }

  //*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
  //************************************************************************************

  void LargeStrain3DLaw::GetLawFeatures(Features& rFeatures)
  {
    KRATOS_TRY

    //Set the type of law
    rFeatures.mOptions.Set( THREE_DIMENSIONAL_LAW );
    rFeatures.mOptions.Set( FINITE_STRAINS );
    rFeatures.mOptions.Set( ISOTROPIC );

    //Get model features
    GetModelFeatures(rFeatures);

    //Set strain measure required by the consitutive law
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

    //Set the strain size
    rFeatures.mStrainSize = GetStrainSize();

    //Set the spacedimension
    rFeatures.mSpaceDimension = WorkingSpaceDimension();

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************

  void LargeStrain3DLaw::GetModelFeatures(Features& rFeatures)
  {
    KRATOS_TRY

    //Get model variables and set law characteristics
    if( mpModel != NULL ){

      std::vector<Variable<double> > ScalarVariables;
      std::vector<Variable<array_1d<double,3> > > ComponentVariables;

      mpModel->GetDomainVariablesList(ScalarVariables, ComponentVariables);

      for(std::vector<Variable<array_1d<double,3> > >::iterator cv_it=ComponentVariables.begin(); cv_it != ComponentVariables.end(); ++cv_it)
	{
	  if( *cv_it == DISPLACEMENT ){
	    for(std::vector<Variable<double> >::iterator sv_it=ScalarVariables.begin(); sv_it != ScalarVariables.end(); ++sv_it)
	      {
		if( *sv_it == PRESSURE )
		  rFeatures.mOptions.Set( U_P_LAW );
	      }
	  }
	  // if( *cv_it == VELOCITY ){
	  //   for(std::vector<Variables<double> >::iterator sv_it=ScalarVariables.begin(); sv_it != ScalarVariables.end(); ++sv_it)
	  //     {
	  // 	if( *sv_it == PRESSURE )
	  // 	  rFeatures.mOptions.Set( V_P_LAW );
	  //     }
	  // }
	}

      //...
    }


    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************

  int LargeStrain3DLaw::Check(const Properties& rProperties,
			       const GeometryType& rElementGeometry,
			       const ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY


    mpModel->Check(rProperties,rCurrentProcessInfo);

    return 0;

    KRATOS_CATCH(" ")
  }


} // Namespace Kratos
