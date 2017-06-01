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
    mDeterminantF0 = 1.0;

    MatrixType Identity = identity_matrix<double>(3);   
    noalias(mInverseDeformationGradientF0) = Identity;

    
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
    mDeterminantF0 = 1.0;

    MatrixType Identity = identity_matrix<double>(3);    
    noalias(mInverseDeformationGradientF0) = Identity;
    
    KRATOS_CATCH(" ")    
  }
  
  //******************************COPY CONSTRUCTOR**************************************
  //************************************************************************************

  LargeStrain3DLaw::LargeStrain3DLaw(const LargeStrain3DLaw& rOther)
    : Constitutive3DLaw(rOther)
    ,mDeterminantF0(rOther.mDeterminantF0)
    ,mInverseDeformationGradientF0(rOther.mInverseDeformationGradientF0)
  {
    mpModel = rOther.mpModel->Clone();
  }

  //*******************************ASSIGMENT OPERATOR***********************************
  //************************************************************************************

  LargeStrain3DLaw& LargeStrain3DLaw::operator=(const LargeStrain3DLaw& rOther)
  {
    Constitutive3DLaw::operator=(rOther);
    mpModel = rOther.mpModel->Clone();
    mDeterminantF0 = rOther.mDeterminantF0;
    mInverseDeformationGradientF0 = rOther.mInverseDeformationGradientF0;
    return *this;
  } 
  
  //********************************CLONE***********************************************
  //************************************************************************************

  ConstitutiveLaw::Pointer LargeStrain3DLaw::Clone() const
  {
    return ( LargeStrain3DLaw::Pointer(new LargeStrain3DLaw(*this)) );
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

    if(rThisVariable == DETERMINANT_F){
      mDeterminantF0 = rValue;
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
      rValue = mDeterminantF0;
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
  
  void LargeStrain3DLaw::InitializeModelData(Parameters& rValues,ModelDataType& rModelValues)
  {
    KRATOS_TRY

    rModelValues.SetOptions(rValues.GetOptions());
    rModelValues.SetMaterialProperties(rValues.GetMaterialProperties());
    rModelValues.SetProcessInfo(rValues.GetProcessInfo());
    rModelValues.SetVoigtSize(this->GetStrainSize());
    rModelValues.SetVoigtIndexTensor(this->GetVoigtIndexTensor());

    LawDataType& rVariables = rModelValues.rConstitutiveLawData();
       
    //a.- Calculate incremental deformation gradient determinant
    rVariables.DeterminantF0 = rValues.GetDeterminantF();    
    rVariables.DeterminantF  = rVariables.DeterminantF0/mDeterminantF0; //determinant incremental F
        
    //b.- Calculate incremental deformation gradient
    const MatrixType& rDeformationGradientF0 = rValues.GetDeformationGradientF();

    noalias(rVariables.DeformationGradientF0) = ConstitutiveModelUtilities::DeformationGradientTo3D(rDeformationGradientF0, rVariables.DeformationGradientF0);
    rVariables.DeformationGradientF = prod(rVariables.DeformationGradientF0, mInverseDeformationGradientF0); //incremental F
        
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
      
      const Matrix& rDeformationGradientF  = rValues.GetDeformationGradientF();
      const double& rDeterminantF          = rValues.GetDeterminantF();
            
      //update total deformation gradient
      MatrixType DeformationGradientF0;
      noalias(DeformationGradientF0) = ConstitutiveModelUtilities::DeformationGradientTo3D(rDeformationGradientF,DeformationGradientF0);
      ConstitutiveModelUtilities::InvertMatrix3( DeformationGradientF0, mInverseDeformationGradientF0, mDeterminantF0);
      mDeterminantF0 = rDeterminantF; //special treatment of the determinant
	
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

    //0.- Check if the constitutive parameters are passed correctly to the law calculation
    CheckParameters(rValues);
    
    const Flags& rOptions = rValues.GetOptions();
    
    //1.- Initialize hyperelastic model parameters    
    ModelDataType ModelValues;

    LawDataType& rVariables = ModelValues.rConstitutiveLawData();
    rVariables.StressMeasure = ConstitutiveModelData::StressMeasure_PK2; //required stress measure

    this->InitializeModelData(rValues, ModelValues);
    
    //2.-Get problem variables (Temperature, Pressure, Size) and calculate material parameters
    this->CalculateDomainVariables(rValues, ModelValues);

    ConstitutiveModelData::CalculateMaterialParameters(ModelValues);    

    //3.-Calculate Total PK2 stress and  Constitutive Matrix related to Total PK2 stress
    if(rOptions.Is(ConstitutiveLaw::COMPUTE_STRESS) && rOptions.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)){

      Vector& rStressVector       = rValues.GetStressVector();
      Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix();
      
      this->CalculateStressVectorAndConstitutiveMatrix(ModelValues, rStressVector, rConstitutiveMatrix);

    }
    else{

      //4.-Calculate Total PK2 stress
      if(rOptions.Is(ConstitutiveLaw::COMPUTE_STRESS)){
	
	Vector& rStressVector       = rValues.GetStressVector();
	this->CalculateStressVector(ModelValues, rStressVector);
	
      }

      //5.-Calculate Constitutive Matrix related to Total PK2 stress
      if(rOptions.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)){
	
      	Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix();
	this->CalculateConstitutiveMatrix(ModelValues, rConstitutiveMatrix);
	
      }
 
    } 
    
    //6.- Finalize hyperelastic model parameters    
    this->FinalizeModelData(rValues,ModelValues);

    if(rOptions.Is(ConstitutiveLaw::COMPUTE_STRAIN_ENERGY))
      {
	     
      }  

    KRATOS_CATCH(" ")
      
  }


  //************************************************************************************
  //************************************************************************************

  void LargeStrain3DLaw::CalculateMaterialResponseKirchhoff(Parameters& rValues)
  {
    KRATOS_TRY
 
    //0.- Check if the constitutive parameters are passed correctly to the law calculation
    CheckParameters(rValues);

    const Flags& rOptions = rValues.GetOptions();
    
    //1.- Initialize hyperelastic model parameters    
    ModelDataType ModelValues;

    LawDataType& rVariables = ModelValues.rConstitutiveLawData();
    rVariables.StressMeasure = ConstitutiveModelData::StressMeasure_Kirchhoff; //set required stress measure
    
    this->InitializeModelData(rValues, ModelValues);

    //2.-Calculate domain variables (Temperature, Pressure, Size) and calculate material parameters
    this->CalculateDomainVariables(rValues, ModelValues);

    ConstitutiveModelData::CalculateMaterialParameters(ModelValues);
    
    //3.-Calculate Total kirchhoff stress and  Constitutive Matrix related to Total Kirchhoff stress

    if(rOptions.Is(ConstitutiveLaw::COMPUTE_STRESS) && rOptions.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)){

      Vector& rStressVector       = rValues.GetStressVector();
      Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix();

      this->CalculateStressVectorAndConstitutiveMatrix(ModelValues, rStressVector, rConstitutiveMatrix);

    }
    else{

      //4.-Calculate Total Kirchhoff stress

      if(rOptions.Is(ConstitutiveLaw::COMPUTE_STRESS)){
	
	Vector& rStressVector       = rValues.GetStressVector();
	this->CalculateStressVector(ModelValues, rStressVector);
	
      }

      //5.-Calculate Constitutive Matrix related to Total Kirchhoff stress

      if(rOptions.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)){
	
      	Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix();
	this->CalculateConstitutiveMatrix(ModelValues, rConstitutiveMatrix);
	
      }
 
    }
    
    if(rOptions.Is(ConstitutiveLaw::COMPUTE_STRAIN_ENERGY))
      {
	     
      }  
    

    //6.- Finalize hyperelastic model parameters    
    this->FinalizeModelData(rValues,ModelValues);


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

    rStressVector.clear();
    rStressVector = ConstitutiveModelUtilities::StressTensorToVector(StressMatrix, rStressVector);
        
    KRATOS_CATCH(" ")
  }
  
  //***********************COMPUTE ALGORITHMIC CONSTITUTIVE MATRIX**********************
  //************************************************************************************

  void LargeStrain3DLaw::CalculateConstitutiveMatrix(ModelDataType& rModelValues, Matrix& rConstitutiveMatrix)
  {
    KRATOS_TRY
      
    rConstitutiveMatrix.clear();
    
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
    rConstitutiveMatrix.clear();
    
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
      
      for(std::vector<Variable<array_1d<double,3> > >::iterator cv_it=ComponentVariables.begin(); cv_it != ComponentVariables.end(); cv_it++)
	{
	  if( *cv_it == DISPLACEMENT ){
	    for(std::vector<Variable<double> >::iterator sv_it=ScalarVariables.begin(); sv_it != ScalarVariables.end(); sv_it++)
	      {
		if( *sv_it == PRESSURE )
		  rFeatures.mOptions.Set( U_P_LAW );
	      }
	  }
	  // if( *cv_it == VELOCITY ){
	  //   for(std::vector<Variables<double> >::iterator sv_it=ScalarVariables.begin(); sv_it != ScalarVariables.end(); )
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

  int LargeStrain3DLaw::Check(const Properties& rMaterialProperties,
			       const GeometryType& rElementGeometry,
			       const ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY
      
    if(YOUNG_MODULUS.Key() == 0 || rMaterialProperties[YOUNG_MODULUS]<= 0.00)
      KRATOS_ERROR << "YOUNG_MODULUS has Key zero or invalid value" << std::endl;

    const double& nu = rMaterialProperties[POISSON_RATIO];
    const bool check = bool( (nu >0.499 && nu<0.501 ) || (nu < -0.999 && nu > -1.01 ) );

    if(POISSON_RATIO.Key() == 0 || check==true)
      KRATOS_ERROR << "POISSON_RATIO has Key zero invalid value" << std::endl;


    if(DENSITY.Key() == 0 || rMaterialProperties[DENSITY]<0.00)
      KRATOS_ERROR << "DENSITY has Key zero or invalid value" << std::endl;

    mpModel->Check(rMaterialProperties,rCurrentProcessInfo);
    
    return 0;
    
    KRATOS_CATCH(" ")
  }

  
} // Namespace Kratos
