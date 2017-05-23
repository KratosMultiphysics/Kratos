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
    mStrainVector.clear();
    mStrainVector[0] = 1.0;
    mStrainVector[1] = 1.0;    
    mStrainVector[2] = 1.0;
    
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
    mStrainVector.clear();
    mStrainVector[0] = 1.0;
    mStrainVector[1] = 1.0;    
    mStrainVector[2] = 1.0;
    
    KRATOS_CATCH(" ")    
  }
  
  //******************************COPY CONSTRUCTOR**************************************
  //************************************************************************************

  LargeStrain3DLaw::LargeStrain3DLaw(const LargeStrain3DLaw& rOther)
    : Constitutive3DLaw(rOther)
    ,mDeterminantF0(rOther.mDeterminantF0)
    ,mInverseDeformationGradientF0(rOther.mInverseDeformationGradientF0)
    ,mStrainVector(rOther.mStrainVector)
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
    mStrainVector = rOther.mStrainVector;    
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

    // A method to compute the initial linear strain from the stress is needed
    //if(rThisVariable == INITIAL_STRESS_VECTOR)

    // A method to compute the initial linear strain from the stress is needed
    // if(rThisVariable == INITIAL_STRAIN_VECTOR){
    //   mStrainVector = rValue;
    // }
         
    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************
  
  void LargeStrain3DLaw::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
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

    // if there is no initial strain and no plasticity
    // rVariables.StressMeasure = ConstitutiveModelData::StressMeasure_PK2;        //required stress measure
    // rVariables.StrainMeasure = ConstitutiveModelData::CauchyGreen_None;         //provided strain measure

    // const Matrix& rDeformationGradientF = rValues.GetDeformationGradientF();   //total deformation gradient    
    // noalias(rVariables.DeformationGradientF) = ConstitutiveModelUtilities::DeformationGradientTo3D(rDeformationGradientF,rVariables.DeformationGradientF);
    // rVariables.DeterminantF  = rValues.GetDeterminantF();
    

    rVariables.StressMeasure = ConstitutiveModelData::StressMeasure_Kirchhoff; //required stress measure
    rVariables.StrainMeasure = ConstitutiveModelData::CauchyGreen_Left;        //provided strain measure
   
    //a.- Calculate incremental deformation gradient determinant
    rVariables.DeterminantF  = rValues.GetDeterminantF();    
    rVariables.DeterminantF /= mDeterminantF0; //determinant incremental F
        
    //b.- Calculate incremental deformation gradient
    const MatrixType& rDeformationGradientF = rValues.GetDeformationGradientF();

    noalias(rVariables.DeformationGradientF) = ConstitutiveModelUtilities::DeformationGradientTo3D(rDeformationGradientF, rVariables.DeformationGradientF);
    rVariables.DeformationGradientF = prod(rVariables.DeformationGradientF, mInverseDeformationGradientF0); //incremental F
    noalias(rVariables.IncrementalDeformationGradientF) = rVariables.DeformationGradientF;
    
    if( rValues.GetOptions().Is(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE) )
      rModelValues.State.Set(ConstitutiveModelData::UPDATE_INTERNAL_VARIABLES);
    
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
	
      //update total strain measure
      mStrainVector = ConstitutiveModelUtilities::SymmetricTensorToVector(rModelValues.StrainMatrix, mStrainVector);

      // //update total strain measure ( in the UpdateInternalVariables method of the plasticity model )
      // mStrainVector[0] = rModelValues.StressMatrix(0,0);
      // mStrainVector[1] = rModelValues.StressMatrix(1,1);
      // mStrainVector[2] = rModelValues.StressMatrix(2,2);
      
      // mStrainVector[3] = rModelValues.StressMatrix(0,1);
      // mStrainVector[4] = rModelValues.StressMatrix(1,2);
      // mStrainVector[5] = rModelValues.StressMatrix(2,0);
      
      // mStrainVector   *=  ( 1.0 / rModelValues.MaterialParameters.LameMu );
      
      // double VolumetricPart = (rModelValues.StrainMatrix(0,0)+rModelValues.StrainMatrix(1,1)+rModelValues.StrainMatrix(2,2))/3.0;
      
      // mStrainVector[0] += VolumetricPart;
      // mStrainVector[1] += VolumetricPart;
      // mStrainVector[2] += VolumetricPart;

      
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
    this->InitializeModelData(rValues, ModelValues);

    
    // Calculate incremental right cauchy green tensor
    LawDataType& rVariables = ModelValues.rConstitutiveLawData();
    rVariables.StressMeasure = ConstitutiveModelData::StressMeasure_PK2; //required stress measure
    rVariables.StrainMeasure = ConstitutiveModelData::CauchyGreen_Right;  //provided strain measure

    ModelValues.StrainMatrix = ConstitutiveModelUtilities::VectorToSymmetricTensor(mStrainVector, ModelValues.StrainMatrix);

    ModelValues.StrainMatrix = prod(ModelValues.StrainMatrix,rVariables.DeformationGradientF);
    ModelValues.StrainMatrix = prod(trans(rVariables.DeformationGradientF),ModelValues.StrainMatrix);

    // Set Total DeterminantF and DeformationGradientF
    rVariables.DeterminantF         = rValues.GetDeterminantF();
    rVariables.DeformationGradientF = rValues.GetDeformationGradientF();   
    
    //2.-Get problem variables (Temperature, Pressure, Size) and calculate material parameters
    this->CalculateDomainVariables(rValues, ModelValues);

    ConstitutiveModelData::CalculateMaterialParameters(ModelValues);
    
    //3.-Calculate Total Strain
    
    if(rOptions.Is(ConstitutiveLaw::COMPUTE_STRAIN)) //large strains
      {
	Vector& rStrainVector   = rValues.GetStrainVector();
	
	//E= 0.5*(C-1) Green-Lagrange Strain
	ConstitutiveModelUtilities::RightCauchyToGreenLagrangeStrain(ModelValues.StrainMatrix,rStrainVector);

	//LawDataType& rVariables = ModelValues.rConstitutiveLawData();
	//E= 0.5*(FT*F-1) Green-Lagrange Strain
	//ConstitutiveModelUtilities::CalculateGreenLagrangeStrain(rVariables.DeformationGradientF,rStrainVector);
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

  void LargeStrain3DLaw::CalculateMaterialResponseKirchhoff(Parameters& rValues)
  {
    KRATOS_TRY
 
    //0.- Check if the constitutive parameters are passed correctly to the law calculation
    CheckParameters(rValues);

    const Flags& rOptions = rValues.GetOptions();
    
    //1.- Initialize hyperelastic model parameters    
    ModelDataType ModelValues;
    this->InitializeModelData(rValues, ModelValues);

    // Calculate incremental left cauchy green tensor
    LawDataType& rVariables = ModelValues.rConstitutiveLawData();
    rVariables.StressMeasure = ConstitutiveModelData::StressMeasure_Kirchhoff; //required stress measure
    rVariables.StrainMeasure = ConstitutiveModelData::CauchyGreen_Left;  //provided strain measure

    ModelValues.StrainMatrix = ConstitutiveModelUtilities::VectorToSymmetricTensor(mStrainVector,ModelValues.StrainMatrix);
    
    ModelValues.StrainMatrix = prod(ModelValues.StrainMatrix,trans(rVariables.DeformationGradientF));
    ModelValues.StrainMatrix = prod(rVariables.DeformationGradientF,ModelValues.StrainMatrix);

    
    // Set Total DeterminantF and DeformationGradientF
    rVariables.DeterminantF         = rValues.GetDeterminantF();
    rVariables.DeformationGradientF = rValues.GetDeformationGradientF();   
    
    //2.-Calculate domain variables (Temperature, Pressure, Size) and calculate material parameters
    this->CalculateDomainVariables(rValues, ModelValues);

    ConstitutiveModelData::CalculateMaterialParameters(ModelValues);
    
    //3.-Calculate Total Strain
    
    if(rOptions.Is(ConstitutiveLaw::COMPUTE_STRAIN)) //large strains
      {
        Vector& rStrainVector = rValues.GetStrainVector();
	
	// e= 0.5*(1-inv(b))
        MatrixType TotalLeftCauchyGreen = prod( rVariables.DeformationGradientF, trans( rVariables.DeformationGradientF ));
        ConstitutiveModelUtilities::LeftCauchyToAlmansiStrain( TotalLeftCauchyGreen, rStrainVector);

	//LawDataType& rVariables = ModelValues.rConstitutiveLawData();
        //e= 0.5*(1-invFT*invF) Almansi Strain
        //ConstitutiveModelUtilities::CalculateAlmansiStrain(rVariables.DeformationGradientF,StrainVector);

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
