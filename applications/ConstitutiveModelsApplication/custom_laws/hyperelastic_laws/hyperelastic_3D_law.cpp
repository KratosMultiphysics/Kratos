//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                December 2016 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_laws/hyperelastic_laws/hyperelastic_3D_law.hpp"

namespace Kratos
{
  
  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

  HyperElastic3DLaw::HyperElastic3DLaw()
    : Elastic3DLaw()
  {
    KRATOS_TRY

    //member variables initialization
    mDeterminantF0 = 1.0;

    MatrixType Identity = identity_matrix<double>(3);   
    noalias(mInverseDeformationGradientF0) = Identity;
    mCauchyGreenVector.clear();
    mCauchyGreenVector[0] = 1.0;
    mCauchyGreenVector[1] = 1.0;    
    mCauchyGreenVector[2] = 1.0;
    
    KRATOS_CATCH(" ")    
  }
  
  //******************************CONSTRUCTOR WITH THE MODEL****************************
  //************************************************************************************

  HyperElastic3DLaw::HyperElastic3DLaw(ModelTypePointer pModel)
    : Elastic3DLaw()
    , mpModel(pModel)
  {
    KRATOS_TRY

    //member variables initialization
    mDeterminantF0 = 1.0;

    MatrixType Identity = identity_matrix<double>(3);    
    noalias(mInverseDeformationGradientF0) = Identity;
    mCauchyGreenVector.clear();
    mCauchyGreenVector[0] = 1.0;
    mCauchyGreenVector[1] = 1.0;    
    mCauchyGreenVector[2] = 1.0;
    
    KRATOS_CATCH(" ")    
  }
  
  //******************************COPY CONSTRUCTOR**************************************
  //************************************************************************************

  HyperElastic3DLaw::HyperElastic3DLaw(const HyperElastic3DLaw& rOther)
    : Elastic3DLaw(rOther)
    , mpModel(rOther.mpModel)
    , mDeterminantF0(rOther.mDeterminantF0)  
    , mInverseDeformationGradientF0(rOther.mInverseDeformationGradientF0)    
    , mCauchyGreenVector(rOther.mCauchyGreenVector)
  {
  }

  //********************************CLONE***********************************************
  //************************************************************************************

  ConstitutiveLaw::Pointer HyperElastic3DLaw::Clone() const
  {
    return ( HyperElastic3DLaw::Pointer(new HyperElastic3DLaw(*this)) );
  }

  
  //*******************************DESTRUCTOR*******************************************
  //************************************************************************************

  HyperElastic3DLaw::~HyperElastic3DLaw()
  {
  }

  //*******************************OPERATIONS FROM BASE CLASS***************************
  //************************************************************************************

  //***********************HAS : DOUBLE - VECTOR - MATRIX*******************************
  //************************************************************************************
  
    bool HyperElastic3DLaw::Has( const Variable<double>& rThisVariable )
  {
    KRATOS_TRY
      
    if(rThisVariable == DETERMINANT_F)
      return true;
    
    return false;
    
    KRATOS_CATCH(" ")
  }
  
 
  //***********************SET VALUE: DOUBLE - VECTOR - MATRIX**************************
  //************************************************************************************

  void HyperElastic3DLaw::SetValue( const Variable<double>& rThisVariable, const double& rValue,
				    const ProcessInfo& rCurrentProcessInfo )
  {
    KRATOS_TRY

    if(rThisVariable == DETERMINANT_F){
      mDeterminantF0 = rValue;
    }

    KRATOS_CATCH(" ")
  }
  
  //************************************************************************************
  //************************************************************************************
  
  void HyperElastic3DLaw::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
				    const ProcessInfo& rCurrentProcessInfo )
  {
    KRATOS_TRY

    // A method to compute the initial linear strain from the stress is needed
    //if(rThisVariable == INITIAL_STRESS_VECTOR)

    // A method to compute the initial linear strain from the stress is needed
    // if(rThisVariable == INITIAL_STRAIN_VECTOR){
    //   mCauchyGreenVector = rValue;
    // }
         
    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************
  
  void HyperElastic3DLaw::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
				    const ProcessInfo& rCurrentProcessInfo )
  {
    KRATOS_TRY

    // A method to compute the initial linear strain from the stress is needed
    //if(rThisVariable == INITIAL_STRESS_VECTOR)

    // A method to compute the initial linear strain from the stress is needed
    //if(rThisVariable == INITIAL_STRAIN_MATRIX){
         
    KRATOS_CATCH(" ")
  }

  
  //***********************GET VALUE: DOUBLE - VECTOR - MATRIX**************************
  //************************************************************************************

  double& HyperElastic3DLaw::GetValue( const Variable<double>& rThisVariable, double& rValue )
  {
    KRATOS_TRY
      
    if(rThisVariable == DETERMINANT_F){
      rValue = mDeterminantF0;
    }
      
    if (rThisVariable==PLASTIC_STRAIN){
    }
  
    if (rThisVariable==DELTA_PLASTIC_STRAIN){
    }

    return( rValue );
    
    KRATOS_CATCH(" ")   
  }

  //************* STARTING - ENDING  METHODS
  //************************************************************************************
  //************************************************************************************


  //************* COMPUTING  METHODS
  //************************************************************************************
  //************************************************************************************
  
  void  HyperElastic3DLaw::InitializeModelData(Parameters& rValues,ModelDataType& rModelValues)
  {
    KRATOS_TRY

    rModelValues.SetMaterialProperties(rValues.GetMaterialProperties());
    rModelValues.SetProcessInfo(rValues.GetProcessInfo());
    rModelValues.SetVoigtSize(this->GetStrainSize());
    rModelValues.SetVoigtIndexTensor(this->GetVoigtIndexTensor());

    ConstitutiveLawDataType& rVariables = rModelValues.rConstitutiveLawData();

    // if there is no initial strain and no plasticity
    // rVariables.StressMeasure = ConstitutiveModelData::StressMeasure_PK2;        //required stress measure
    // rVariables.StrainMeasure = ConstitutiveModelData::CauchyGreen_None;         //provided strain measure

    // const Matrix& rDeformationGradientF = rValues.GetDeformationGradientF();   //total deformation gradient    
    // rVariables.DeformationGradientF = ConstitutiveLawUtilities::DeformationGradientTo3D(rVariables.DeformationGradientF, rDeformationGradientF);
    // rVariables.DeterminantF  = rValues.GetDeterminantF();
    

    rVariables.StressMeasure = ConstitutiveModelData::StressMeasure_Kirchhoff; //required stress measure
    rVariables.StrainMeasure = ConstitutiveModelData::CauchyGreen_Left;        //provided strain measure
   
    //a.- Calculate incremental deformation gradient determinant
    rVariables.DeterminantF = rValues.GetDeterminantF();
    
    rVariables.DeterminantF /= mDeterminantF0; //determinant incremental F
        
    //b.- Calculate incremental deformation gradient
    const MatrixType& rDeformationGradientF = rValues.GetDeformationGradientF(); 
    rVariables.DeformationGradientF = ConstitutiveLawUtilities::DeformationGradientTo3D(rVariables.DeformationGradientF, rDeformationGradientF);

    rVariables.DeformationGradientF = prod(rVariables.DeformationGradientF, mInverseDeformationGradientF0); //incremental F
    
    //c.- Calculate incremental left cauchy green tensor
    rModelValues.StrainMatrix = ConstitutiveLawUtilities::VectorToSymmetricTensor(mCauchyGreenVector, rModelValues.StrainMatrix);
    
    rModelValues.StrainMatrix = prod(rModelValues.StrainMatrix,trans(rVariables.DeformationGradientF));
    rModelValues.StrainMatrix = prod(rVariables.DeformationGradientF,rModelValues.StrainMatrix);

    //d.- Set Total DeterminantF and DeformationGracientF
    rVariables.DeterminantF         = rValues.GetDeterminantF();
    rVariables.DeformationGradientF = rValues.GetDeformationGradientF();

    
    if( rValues.GetOptions().Is(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE) )
      rModelValues.State.Set(ConstitutiveModelData::UPDATE_INTERNAL_VARIABLES);
    
    KRATOS_CATCH(" ")      
  }

  //************************************************************************************
  //************************************************************************************

  void  HyperElastic3DLaw::FinalizeModelData(Parameters& rValues,ModelDataType& rModelValues)
  {
    KRATOS_TRY
      
    //Finalize Material response
    if(rValues.GetOptions().Is(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE)){
      const Matrix& rDeformationGradientF    = rValues.GetDeformationGradientF();
      const double& rDeterminantF            = rValues.GetDeterminantF();
      
      //update total deformation gradient
      MatrixType DeformationGradientF0;
      noalias(DeformationGradientF0) = ConstitutiveLawUtilities::DeformationGradientTo3D(DeformationGradientF0,rDeformationGradientF);
      ConstitutiveLawUtilities::InvertMatrix3( DeformationGradientF0, mInverseDeformationGradientF0, mDeterminantF0);
      mDeterminantF0 = rDeterminantF; //special treatment of the determinant
	
      //update total strain measure
      mCauchyGreenVector = ConstitutiveLawUtilities::SymmetricTensorToVector(rModelValues.StrainMatrix, mCauchyGreenVector);
    }
    
    KRATOS_CATCH(" ")      
  }
  
  
  //*****************************MATERIAL RESPONSES*************************************
  //************************************************************************************


  void  HyperElastic3DLaw::CalculateMaterialResponsePK2(Parameters& rValues)
  {
    KRATOS_TRY
 
    //0.- Check if the constitutive parameters are passed correctly to the law calculation
    CheckParameters(rValues);

    const Flags& rOptions = rValues.GetOptions();
    
    //1.- Initialize hyperelastic model parameters    
    ModelDataType ModelValues;
    this->InitializeModelData(rValues, ModelValues);

    //2.-Get problem variables (Temperature, Pressure, Size) and calculate material parameters
    this->CalculateDomainVariables(rValues, ModelValues);

    ConstitutiveModelData::CalculateMaterialParameters(ModelValues);
    
    //3.-Calculate Total Strain
    
    if(rOptions.Is(ConstitutiveLaw::COMPUTE_STRAIN)) //large strains
      {
	Vector& rStrainVector   = rValues.GetStrainVector();
	
	//E= 0.5*(C-1) Green-Lagrange Strain
	ConstitutiveLawUtilities::RightCauchyToGreenLagrangeStrain(ModelValues.StrainMatrix,rStrainVector);

	//ConstitutiveLawDataType& rVariables = ModelValues.rConstitutiveLawData();
	//E= 0.5*(FT*F-1) Green-Lagrange Strain
	//ConstitutiveLawUtilities::CalculateGreenLagrangeStrain(rVariables.DeformationGradientF,rStrainVector);
      }

    //4.-Calculate Total PK2 stress and  Constitutive Matrix related to Total PK2 stress
    if(rOptions.Is(ConstitutiveLaw::COMPUTE_STRESS) && rOptions.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)){

      Vector& rStressVector       = rValues.GetStressVector();
      Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix();

      this->CalculateStressVectorAndConstitutiveMatrix(ModelValues, rStressVector, rConstitutiveMatrix);

    }
    else{

      //5.-Calculate Total PK2 stress
      if(rOptions.Is(ConstitutiveLaw::COMPUTE_STRESS)){
	
	Vector& rStressVector       = rValues.GetStressVector();
	this->CalculateStressVector(ModelValues, rStressVector);
	
      }

      //6.-Calculate Constitutive Matrix related to Total PK2 stress
      if(rOptions.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)){
	
      	Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix();
	this->CalculateConstitutiveMatrix(ModelValues, rConstitutiveMatrix);
	
      }
 
    } 
    
    //7.- Finalize hyperelastic model parameters    
    this->FinalizeModelData(rValues,ModelValues);

    if(rOptions.Is(ConstitutiveLaw::COMPUTE_STRAIN_ENERGY))
      {
	     
      }  

    KRATOS_CATCH(" ")
      
  }


  //************************************************************************************
  //************************************************************************************

  void HyperElastic3DLaw::CalculateMaterialResponseKirchhoff(Parameters& rValues)
  {
    KRATOS_TRY
 
    //0.- Check if the constitutive parameters are passed correctly to the law calculation
    CheckParameters(rValues);

    const Flags& rOptions = rValues.GetOptions();
    
    //1.- Initialize hyperelastic model parameters    
    ModelDataType ModelValues;
    this->InitializeModelData(rValues, ModelValues);

    //2.-Calculate domain variables (Temperature, Pressure, Size) and calculate material parameters
    this->CalculateDomainVariables(rValues, ModelValues);

    ConstitutiveModelData::CalculateMaterialParameters(ModelValues);
    
    //3.-Calculate Total Strain
    
    if(rOptions.Is(ConstitutiveLaw::COMPUTE_STRAIN)) //large strains
      {
        Vector& rStrainVector = rValues.GetStrainVector();
	
	// e= 0.5*(1-inv(b))
        ConstitutiveLawUtilities::InverseLeftCauchyToAlmansiStrain(ModelValues.StrainMatrix,rStrainVector);

	//ConstitutiveLawDataType& rVariables = ModelValues.rConstitutiveLawData();
        //e= 0.5*(1-invFT*invF) Almansi Strain
        //ConstitutiveLawUtilities::CalculateAlmansiStrain(rVariables.DeformationGradientF,StrainVector);
      }

    //4.-Calculate Total kirchhoff stress and  Constitutive Matrix related to Total Kirchhoff stress

    if(rOptions.Is(ConstitutiveLaw::COMPUTE_STRESS) && rOptions.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)){

      Vector& rStressVector       = rValues.GetStressVector();
      Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix();

      this->CalculateStressVectorAndConstitutiveMatrix(ModelValues, rStressVector, rConstitutiveMatrix);

    }
    else{

      //5.-Calculate Total Kirchhoff stress

      if(rOptions.Is(ConstitutiveLaw::COMPUTE_STRESS)){
	
	Vector& rStressVector       = rValues.GetStressVector();
	this->CalculateStressVector(ModelValues, rStressVector);
	
      }

      //6.-Calculate Constitutive Matrix related to Total Kirchhoff stress

      if(rOptions.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)){
	
      	Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix();
	this->CalculateConstitutiveMatrix(ModelValues, rConstitutiveMatrix);
	
      }
 
    }
    
    if(rOptions.Is(ConstitutiveLaw::COMPUTE_STRAIN_ENERGY))
      {
	     
      }  
    

    //7.- Finalize hyperelastic model parameters    
    this->FinalizeModelData(rValues,ModelValues);

    
    KRATOS_CATCH(" ")
      
  }


  //*******************************COMPUTE STRESS VECTOR********************************
  //************************************************************************************

  void HyperElastic3DLaw::CalculateStressVector(ModelDataType& rModelValues, Vector& rStressVector)
  {
    KRATOS_TRY

    MatrixType StressMatrix;
      
    if(rModelValues.GetOptions().Is(ConstitutiveLaw::ISOCHORIC_TENSOR_ONLY)){
      mpModel->CalculateIsochoricStressTensor(rModelValues, StressMatrix);
    }
    else if(rModelValues.GetOptions().Is(ConstitutiveLaw::VOLUMETRIC_TENSOR_ONLY)){     
      mpModel->CalculateVolumetricStressTensor(rModelValues, StressMatrix);
    }
    else{      
      mpModel->CalculateStressTensor(rModelValues, StressMatrix);
    }

    rStressVector = MathUtils<double>::StressTensorToVector(StressMatrix, rStressVector.size());
      
    KRATOS_CATCH(" ")
  }
  
  //***********************COMPUTE ALGORITHMIC CONSTITUTIVE MATRIX**********************
  //************************************************************************************

  void HyperElastic3DLaw::CalculateConstitutiveMatrix(ModelDataType& rModelValues, Matrix& rConstitutiveMatrix)
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

  void HyperElastic3DLaw::CalculateStressVectorAndConstitutiveMatrix(ModelDataType& rModelValues, Vector& rStressVector, Matrix& rConstitutiveMatrix)
  {
    KRATOS_TRY

    MatrixType StressMatrix;
    
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
        
    KRATOS_CATCH(" ")
  }
  
  //*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
  //************************************************************************************

  void HyperElastic3DLaw::GetLawFeatures(Features& rFeatures)
  {
    KRATOS_TRY
    
    //Set the type of law
    rFeatures.mOptions.Set( THREE_DIMENSIONAL_LAW );
    rFeatures.mOptions.Set( FINITE_STRAINS );
    rFeatures.mOptions.Set( ISOTROPIC );

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

  int HyperElastic3DLaw::Check(const Properties& rMaterialProperties,
			       const GeometryType& rElementGeometry,
			       const ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY
      
    if(YOUNG_MODULUS.Key() == 0 || rMaterialProperties[YOUNG_MODULUS]<= 0.00)
      KRATOS_THROW_ERROR( std::invalid_argument,"YOUNG_MODULUS has Key zero or invalid value ", "" )

    const double& nu = rMaterialProperties[POISSON_RATIO];
    const bool check = bool( (nu >0.499 && nu<0.501 ) || (nu < -0.999 && nu > -1.01 ) );

    if(POISSON_RATIO.Key() == 0 || check==true)
      KRATOS_THROW_ERROR( std::invalid_argument,"POISSON_RATIO has Key zero invalid value ", "" )


    if(DENSITY.Key() == 0 || rMaterialProperties[DENSITY]<0.00)
      KRATOS_THROW_ERROR( std::invalid_argument,"DENSITY has Key zero or invalid value ", "" )

    mpModel->Check(rMaterialProperties,rCurrentProcessInfo);
    
    return 0;
    
    KRATOS_CATCH(" ")
  }

  
} // Namespace Kratos
