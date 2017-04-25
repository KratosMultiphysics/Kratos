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
#include "custom_models/elasticity_models/hyperelastic_model.hpp"


namespace Kratos
{

  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

  HyperElasticModel::HyperElasticModel()
    : ConstitutiveModel()
    , msIdentityMatrix( identity_matrix<double>(3) )
  {
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
    return ( HyperElasticModel::Pointer(new HyperElasticModel(*this)) );
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

  void HyperElasticModel::CalculateStrainData(ModelDataType& rValues, HyperElasticDataType& rVariables)
  {
    KRATOS_TRY

    //set model data pointer
    rVariables.SetModelData(rValues);
    rVariables.SetState(rValues.State);
    
    //cauchy green tensor
    const MatrixType& rStrainMatrix         = rValues.GetStrainMatrix();
    const MatrixType& rDeformationGradientF = rValues.GetDeformationGradientF();

    const StrainMeasureType& rStrainMeasure = rValues.GetStrainMeasure();
    const StressMeasureType& rStressMeasure = rValues.GetStressMeasure();
    
    if( rStressMeasure == ConstitutiveModelData::StressMeasure_PK2 ){ //mCauchyGreenMatrix = RightCauchyGreen (C=FT*F)  C^-1=(FT*F)^-1=F^-1*FT^-1
      
      if( rStrainMeasure == ConstitutiveModelData::CauchyGreen_Right ){
	noalias(rVariables.Strain.CauchyGreenMatrix) = rStrainMatrix;
	ConstitutiveModelUtilities::InvertMatrix3( rVariables.Strain.CauchyGreenMatrix, rVariables.Strain.InverseCauchyGreenMatrix, rVariables.Strain.Invariants.I3 );
	rValues.State.Set(ConstitutiveModelData::COMPUTED_STRAIN);
      }
      else if( rStrainMeasure == ConstitutiveModelData::CauchyGreen_None ){
	noalias(rVariables.Strain.CauchyGreenMatrix) = prod(trans(rDeformationGradientF), rDeformationGradientF);
	ConstitutiveModelUtilities::InvertMatrix3( rVariables.Strain.CauchyGreenMatrix, rVariables.Strain.InverseCauchyGreenMatrix, rVariables.Strain.Invariants.I3 );

	rValues.StrainMatrix = rVariables.Strain.CauchyGreenMatrix;
	rValues.State.Set(ConstitutiveModelData::COMPUTED_STRAIN);
      }
      else{
	KRATOS_ERROR << "calling initialize HyperElasticModel .. StrainMeasure provided is inconsistent" << std::endl;
      }
      
    }
    else if( rStressMeasure == ConstitutiveModelData::StressMeasure_Kirchhoff ){ //mCauchyGreenMatrix = LeftCauchyGreen (b=F*FT)

      if( rStrainMeasure == ConstitutiveModelData::CauchyGreen_Left ){
	noalias(rVariables.Strain.CauchyGreenMatrix) = rStrainMatrix;
	rValues.State.Set(ConstitutiveModelData::COMPUTED_STRAIN);
      }
      else if( rStrainMeasure == ConstitutiveModelData::CauchyGreen_None ){
	noalias(rVariables.Strain.CauchyGreenMatrix) = prod(rDeformationGradientF,trans(rDeformationGradientF));

	rValues.StrainMatrix = rVariables.Strain.CauchyGreenMatrix;
	rValues.State.Set(ConstitutiveModelData::COMPUTED_STRAIN);
      }
      else{
	KRATOS_ERROR << "calling initialize HyperElasticModel .. StrainMeasure provided is inconsistent" << std::endl;
      }

      rValues.MaterialParameters.LameMuBar =  rValues.MaterialParameters.LameMu * ( rVariables.Strain.CauchyGreenMatrix(0,0) + rVariables.Strain.CauchyGreenMatrix(1,1) + rVariables.Strain.CauchyGreenMatrix(2,2) ) * 1.0/3.0 * pow(rValues.GetDeterminantF(),(-2.0/3.0));
      
    }
    else{
      KRATOS_ERROR << "calling initialize HyperElasticModel .. StressMeasure required is inconsistent"  << std::endl;
    }

    
    this->CalculateInvariants(rVariables);

    //Algorithmic moduli factors
    this->CalculateScalingFactors(rVariables);

      
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

    const ModelDataType&  rModelData        = rVariables.GetModelData();
    const StressMeasureType& rStressMeasure = rModelData.GetStressMeasure();
    
    MatrixType StressPartMatrix;
    MatrixType StressMatrix;
	    
    if( rStressMeasure == ConstitutiveModelData::StressMeasure_PK2 ){ //mCauchyGreenMatrix = RightCauchyGreen (C)

      StressPartMatrix = GetI1RightCauchyGreenDerivative(rVariables.Strain,StressPartMatrix);
      noalias(StressMatrix)  = rVariables.Factors.Alpha1 * StressPartMatrix;

      StressPartMatrix = GetI2RightCauchyGreenDerivative(rVariables.Strain,StressPartMatrix);
      noalias(StressMatrix) += rVariables.Factors.Alpha2 * StressPartMatrix;

      StressPartMatrix = GetI3RightCauchyGreenDerivative(rVariables.Strain,StressPartMatrix);
      noalias(StressMatrix) += rVariables.Factors.Alpha3 * StressPartMatrix;

      StressMatrix *= 2.0;
      
      rStressMatrix += StressMatrix;
      
    }
    else if( rStressMeasure == ConstitutiveModelData::StressMeasure_Kirchhoff ){ //mCauchyGreenMatrix = LeftCauchyGreen (b)
     
      StressPartMatrix = GetI1LeftCauchyGreenDerivative(rVariables.Strain,StressPartMatrix);
      noalias(StressMatrix)  = rVariables.Factors.Alpha1 * StressPartMatrix;
     
      StressPartMatrix = GetI2LeftCauchyGreenDerivative(rVariables.Strain,StressPartMatrix);
      noalias(StressMatrix) += rVariables.Factors.Alpha2 * StressPartMatrix;
      
      StressPartMatrix = GetI3LeftCauchyGreenDerivative(rVariables.Strain,StressPartMatrix);      
      noalias(StressMatrix) += rVariables.Factors.Alpha3 * StressPartMatrix;

      StressMatrix *= 2.0;

      rStressMatrix += StressMatrix;
    }

    
    rVariables.State().Set(ConstitutiveModelData::COMPUTED_STRESS);
    
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
    rConstitutiveMatrix.clear();
  
    HyperElasticDataType Variables;
    this->CalculateStrainData(rValues,Variables);

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

    rVariables.State().Set(ConstitutiveModelData::COMPUTED_CONSTITUTIVE_MATRIX);
    
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

    //rVariables.State().Set(ConstitutiveModelData::COMPUTED_CONSTITUTIVE_MATRIX);
    
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

    //rVariables.State().Set(ConstitutiveModelData::COMPUTED_CONSTITUTIVE_MATRIX);
    
    KRATOS_CATCH(" ")
  }
  
  //************************************************************************************
  //************************************************************************************
  
  double& HyperElasticModel::AddConstitutiveComponent(HyperElasticDataType& rVariables, double &rCabcd,
						      const unsigned int& a, const unsigned int& b,
						      const unsigned int& c, const unsigned int& d)
  {
    KRATOS_TRY

    //this is a simplified version, more terms of the derivatives are needed to be general:
      
    const ModelDataType&  rModelData        = rVariables.GetModelData();
    const StressMeasureType& rStressMeasure = rModelData.GetStressMeasure();
    
    double Cabcd = 0;
    double Dabcd = 0;

    if( rStressMeasure == ConstitutiveModelData::StressMeasure_PK2 ){ //mCauchyGreenMatrix = LeftCauchyGreen (C)

      //2nd derivatives
      Dabcd = GetI1RightCauchyGreen2ndDerivative(rVariables.Strain,Dabcd,a,b,c,d);
      Cabcd += rVariables.Factors.Alpha1 * Dabcd;

      Dabcd = GetI2RightCauchyGreen2ndDerivative(rVariables.Strain,Dabcd,a,b,c,d);
      Cabcd += rVariables.Factors.Alpha2 * Dabcd;

      Dabcd = GetI3RightCauchyGreen2ndDerivative(rVariables.Strain,Dabcd,a,b,c,d);
      Cabcd += rVariables.Factors.Alpha3 * Dabcd;

      //1st derivatives
      Dabcd = GetI1RightCauchyGreenSquare1stDerivative(rVariables.Strain,Dabcd,a,b,c,d);
      Cabcd += rVariables.Factors.Beta1 * Dabcd;
      
      Dabcd = GetI2RightCauchyGreenSquare1stDerivative(rVariables.Strain,Dabcd,a,b,c,d);
      Cabcd += rVariables.Factors.Beta2 * Dabcd;

      Dabcd = GetI3RightCauchyGreenSquare1stDerivative(rVariables.Strain,Dabcd,a,b,c,d);
      Cabcd += rVariables.Factors.Beta3 * Dabcd;

      Cabcd *= 4.0;


    }
    else if( rStressMeasure == ConstitutiveModelData::StressMeasure_Kirchhoff ){ //mCauchyGreenMatrix = LeftCauchyGreen (b)
     
      //2nd derivatives
      Dabcd = GetI1LeftCauchyGreen2ndDerivative(rVariables.Strain,Dabcd,a,b,c,d);
      Cabcd += rVariables.Factors.Alpha1 * Dabcd;
      
      Dabcd = GetI2LeftCauchyGreen2ndDerivative(rVariables.Strain,Dabcd,a,b,c,d);
      Cabcd += rVariables.Factors.Alpha2 * Dabcd;

      Dabcd = GetI3LeftCauchyGreen2ndDerivative(rVariables.Strain,Dabcd,a,b,c,d);
      Cabcd += rVariables.Factors.Alpha3 * Dabcd;
            
      //1st derivatives
      Dabcd = GetI1LeftCauchyGreenSquare1stDerivative(rVariables.Strain,Dabcd,a,b,c,d);
      Cabcd += rVariables.Factors.Beta1 * Dabcd;

      Dabcd = GetI2LeftCauchyGreenSquare1stDerivative(rVariables.Strain,Dabcd,a,b,c,d);
      Cabcd += rVariables.Factors.Beta2 * Dabcd;

      Dabcd = GetI3LeftCauchyGreenSquare1stDerivative(rVariables.Strain,Dabcd,a,b,c,d);
      Cabcd += rVariables.Factors.Beta3 * Dabcd;

      
      Cabcd *= 4.0;

    }

    rCabcd += Cabcd;
    
    rVariables.State().Set(ConstitutiveModelData::COMPUTED_CONSTITUTIVE_MATRIX);
    
    return rCabcd;

    
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

  
  void HyperElasticModel::CalculateStrainInvariants(const MatrixType& rCauchyGreenMatrix, double& rI1, double& rI2, double& rI3)
  {
    KRATOS_TRY
      
    rI1 = rCauchyGreenMatrix(0,0) + rCauchyGreenMatrix(1,1) + rCauchyGreenMatrix(2,2);

    rI2 = (  rCauchyGreenMatrix(1,1)*rCauchyGreenMatrix(2,2)	     
	     + rCauchyGreenMatrix(2,2)*rCauchyGreenMatrix(3,3)
	     + rCauchyGreenMatrix(1,1)*rCauchyGreenMatrix(3,3)
	     - rCauchyGreenMatrix(1,2)*rCauchyGreenMatrix(2,1)
	     - rCauchyGreenMatrix(2,3)*rCauchyGreenMatrix(3,2)
	     - rCauchyGreenMatrix(1,3)*rCauchyGreenMatrix(3,1) );
	
    rI3 = (  rCauchyGreenMatrix(0,0)*rCauchyGreenMatrix(1,1)*rCauchyGreenMatrix(2,2)
	     + rCauchyGreenMatrix(0,1)*rCauchyGreenMatrix(1,2)*rCauchyGreenMatrix(2,0)
	     + rCauchyGreenMatrix(1,0)*rCauchyGreenMatrix(2,1)*rCauchyGreenMatrix(0,2)
	     - rCauchyGreenMatrix(2,0)*rCauchyGreenMatrix(1,1)*rCauchyGreenMatrix(0,2)
	     - rCauchyGreenMatrix(1,0)*rCauchyGreenMatrix(0,1)*rCauchyGreenMatrix(2,2)
	     - rCauchyGreenMatrix(0,0)*rCauchyGreenMatrix(2,1)*rCauchyGreenMatrix(1,2) );

    //std::cout<<" I1: "<<rI1<<" I2: "<<rI2<<" I3: "<<rI3<<std::endl;
    
    KRATOS_CATCH(" ")
  }

  

  //************************************************************************************
  //************************************************************************************
  
  void HyperElasticModel::CalculateInvariants(HyperElasticDataType& rVariables)
  {
    KRATOS_TRY

      
    //invariants
    this->CalculateStrainInvariants( rVariables.Strain.CauchyGreenMatrix, rVariables.Strain.Invariants.I1, rVariables.Strain.Invariants.I2, rVariables.Strain.Invariants.I3 );
 
    //jacobian
    rVariables.Strain.Invariants.J    = rVariables.GetModelData().GetDeterminantF();
    rVariables.Strain.Invariants.J_13 = pow(rVariables.Strain.Invariants.J,(-1.0/3.0));

    //std::cout<<" Strain.Invariants [I1:"<<rVariables.Strain.Invariants.I1<<" I2:"<<rVariables.Strain.Invariants.I2<<" I3:"<<rVariables.Strain.Invariants.I3<<"] J:"<<rVariables.Strain.Invariants.J<<std::endl;
    
    rVariables.Strain.Invariants.I3 = rVariables.Strain.Invariants.J * rVariables.Strain.Invariants.J; //for volumetric consistency
    
            
    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************
   
  void HyperElasticModel::CalculateScalingFactors(HyperElasticDataType& rVariables)
  {
    KRATOS_TRY

    rVariables.Factors.Alpha1 = this->GetFunction1stI1Derivative(rVariables, rVariables.Factors.Alpha1);
    rVariables.Factors.Alpha2 = this->GetFunction1stI2Derivative(rVariables, rVariables.Factors.Alpha2);
    rVariables.Factors.Alpha3 = this->GetFunction1stI3Derivative(rVariables, rVariables.Factors.Alpha3);

    rVariables.Factors.Beta1  = this->GetFunction2ndI1Derivative(rVariables, rVariables.Factors.Beta1);
    rVariables.Factors.Beta2  = this->GetFunction2ndI2Derivative(rVariables, rVariables.Factors.Beta2);
    rVariables.Factors.Beta3  = this->GetFunction2ndI3Derivative(rVariables, rVariables.Factors.Beta3);

    // std::cout<<" Alpha["<<rVariables.Factors.Alpha1<<" "<<rVariables.Factors.Alpha2<<" "<<rVariables.Factors.Alpha3<<"]"<<std::endl;
    // std::cout<<" Beta ["<<rVariables.Factors.Beta1<<" "<<rVariables.Factors.Beta2<<" "<<rVariables.Factors.Beta3<<"]"<<std::endl;
    
    KRATOS_CATCH(" ")
  }

    
  //************// dW
    
  double& HyperElasticModel::GetFunction1stI1Derivative(HyperElasticDataType& rVariables, double& rDerivative) //dW/dI1
  {
    KRATOS_TRY
	
    KRATOS_ERROR << "calling the base class function in HyperElasticModel ... illegal operation" << std::endl;

    return rDerivative;

    KRATOS_CATCH(" ")
  }


  double& HyperElasticModel::GetFunction1stI2Derivative(HyperElasticDataType& rVariables, double& rDerivative) //dW/dI2
  {
    KRATOS_TRY
	
    KRATOS_ERROR << "calling the base class function in HyperElasticModel ... illegal operation" << std::endl;

    return rDerivative;

    KRATOS_CATCH(" ")
  }

  double& HyperElasticModel::GetFunction1stI3Derivative(HyperElasticDataType& rVariables, double& rDerivative) //dW/dI3
  {
    KRATOS_TRY
	
    KRATOS_ERROR << "calling the base class function in HyperElasticModel ... illegal operation" << std::endl;

    return rDerivative;
    
    KRATOS_CATCH(" ")
  }


  double& HyperElasticModel::GetVolumetricFunction1stJDerivative(HyperElasticDataType& rVariables, double& rDerivative) //dU/dJ
  {
    KRATOS_TRY
	
    KRATOS_ERROR << "calling the base class function in HyperElasticModel ... illegal operation" << std::endl;

    return rDerivative;

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************
  
    
  double& HyperElasticModel::GetFunction2ndI1Derivative(HyperElasticDataType& rVariables, double& rDerivative) //ddW/dI1dI1
  {
    KRATOS_TRY
	
    KRATOS_ERROR << "calling the base class function in HyperElasticModel ... illegal operation" << std::endl;

    return rDerivative;

    KRATOS_CATCH(" ")
  }

  double& HyperElasticModel::GetFunction2ndI2Derivative(HyperElasticDataType& rVariables, double& rDerivative) //ddW/dI2dI2
  {
    KRATOS_TRY
	
    KRATOS_ERROR << "calling the base class function in HyperElasticModel ... illegal operation" << std::endl;

    return rDerivative;

    KRATOS_CATCH(" ")
  }

  double& HyperElasticModel::GetFunction2ndI3Derivative(HyperElasticDataType& rVariables, double& rDerivative) //ddW/dI3dI3
  {
    KRATOS_TRY
	
    KRATOS_ERROR << "calling the base class function in HyperElasticModel ... illegal operation" << std::endl;

    return rDerivative;

    KRATOS_CATCH(" ")
  }
    
  //************************************************************************************
  //************************************************************************************

  
  double& HyperElasticModel::GetFourthOrderUnitTensor(double& rValue,
						      const double& a,
						      const double& b,
						      const double& c,
						      const double& d) //ddC/dCdC or ddb/dbdb
  {
    KRATOS_TRY

    rValue = 0.5*(msIdentityMatrix(a,c)*msIdentityMatrix(b,d)+msIdentityMatrix(a,d)*msIdentityMatrix(b,c));

    return rValue;

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************

  //isochoric volumetric slit
    
  double& HyperElasticModel::GetVolumetricFunction2ndJDerivative(HyperElasticDataType& rVariables, double& rDerivative) //ddU/dJdJ
  {
    KRATOS_TRY
	
    KRATOS_ERROR << "calling the base class function in HyperElasticModel ... illegal operation" << std::endl;

    return rDerivative;

    KRATOS_CATCH(" ")
  }
    

  HyperElasticModel::MatrixType& HyperElasticModel::GetJLeftCauchyGreenDerivative(const CauchyGreenData& rData, MatrixType& rDerivative) //dJ/db
  {
    KRATOS_TRY
	
    noalias(rDerivative)  = msIdentityMatrix;
    rDerivative *= rData.Invariants.J * 0.5;
	
    return rDerivative;

    KRATOS_CATCH(" ")
  }

    
  HyperElasticModel::MatrixType& HyperElasticModel::GetIsochoricRightCauchyGreenDerivative(const CauchyGreenData& rData, MatrixType& rDerivative) //dC'/dC
  {
    KRATOS_TRY
	
    noalias(rDerivative)  = rData.InverseCauchyGreenMatrix;
    rDerivative *= -rData.Invariants.I1/3.0;
    noalias(rDerivative) += msIdentityMatrix;

    rDerivative *= rData.Invariants.J_13 * rData.Invariants.J_13;
      
    return rDerivative;

    KRATOS_CATCH(" ")
  }


  double& HyperElasticModel::GetIsochoricRightCauchyGreenDerivative(const CauchyGreenData& rData,
								    double& rDerivative,
								    const double& a,
								    const double& b,
								    const double& c,
								    const double& d) //dC'/dC
  {
    KRATOS_TRY

    rDerivative  = GetFourthOrderUnitTensor(rDerivative,a,b,c,d);
    rDerivative -= rData.InverseCauchyGreenMatrix(a,b)*rData.CauchyGreenMatrix(c,d)/3.0;
    rDerivative *= rData.Invariants.J_13 * rData.Invariants.J_13;

    return rDerivative;

    KRATOS_CATCH(" ")
  }
    
    
  HyperElasticModel::MatrixType& HyperElasticModel::GetIsochoricLeftCauchyGreenDerivative(const CauchyGreenData& rData, MatrixType& rDerivative) //db'/db
  {
    KRATOS_TRY

    noalias(rDerivative)  = msIdentityMatrix;
    rDerivative *= -rData.Invariants.I1/3.0;
    noalias(rDerivative) += rData.CauchyGreenMatrix;

    rDerivative *= rData.Invariants.J_13 * rData.Invariants.J_13;
      
    return rDerivative;

    KRATOS_CATCH(" ")
  }

  double& HyperElasticModel::GetIsochoricLeftCauchyGreenDerivative(const CauchyGreenData& rData,
								   double& rDerivative,
								   const double& a,
								   const double& b,
								   const double& c,
								   const double& d) //db'/db
  {
    KRATOS_TRY

    rDerivative  = GetFourthOrderUnitTensor(rDerivative,a,b,c,d);
    rDerivative -= msIdentityMatrix(a,b)*rData.CauchyGreenMatrix(c,d)/3.0;
    rDerivative *= rData.Invariants.J_13 * rData.Invariants.J_13;

    return rDerivative;

    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************

  
  //************// right cauchy green: C
    
  HyperElasticModel::MatrixType& HyperElasticModel::GetI1RightCauchyGreenDerivative(const CauchyGreenData& rData, MatrixType& rDerivative) //dI1/dC
  {
    KRATOS_TRY

    noalias(rDerivative) = msIdentityMatrix;
	
    return rDerivative;

    KRATOS_CATCH(" ")
  }

  HyperElasticModel::MatrixType& HyperElasticModel::GetI2RightCauchyGreenDerivative(const CauchyGreenData& rData, MatrixType& rDerivative) //dI2/dC
  {
    KRATOS_TRY
	
    noalias(rDerivative)  = msIdentityMatrix;
    rDerivative *= rData.Invariants.I1;
    noalias(rDerivative) += rData.CauchyGreenMatrix;

    return rDerivative;

    KRATOS_CATCH(" ")
  }

  HyperElasticModel::MatrixType& HyperElasticModel::GetI3RightCauchyGreenDerivative(const CauchyGreenData& rData, MatrixType& rDerivative) //dI3/dC
  {
    KRATOS_TRY
	
    noalias(rDerivative)  = rData.InverseCauchyGreenMatrix;
    rDerivative *= rData.Invariants.I3;
	
    return rDerivative;

    KRATOS_CATCH(" ")
  }
   

  HyperElasticModel::MatrixType& HyperElasticModel::GetJRightCauchyGreenDerivative(const CauchyGreenData& rData, MatrixType& rDerivative) //dJ/dC
  {
    KRATOS_TRY
	
    noalias(rDerivative)  = rData.InverseCauchyGreenMatrix;
    rDerivative *= rData.Invariants.J * 0.5;
	
    return rDerivative;

    KRATOS_CATCH(" ")
  }


  double& HyperElasticModel::GetInverseRightCauchyGreenDerivative(const CauchyGreenData& rData,
								  double& rDerivative,
								  const double& a,
								  const double& b,
								  const double& c,
								  const double& d) //dC^-1/dC
  {
    KRATOS_TRY

    rDerivative = -0.5*(rData.InverseCauchyGreenMatrix(a,c)*rData.InverseCauchyGreenMatrix(b,d)+rData.InverseCauchyGreenMatrix(a,d)*rData.InverseCauchyGreenMatrix(b,c));

    return rDerivative;

    KRATOS_CATCH(" ")
  }

    
  //Invariants 1st derivatives by components
  double& HyperElasticModel::GetI1RightCauchyGreen1stDerivative(const CauchyGreenData& rData,
								double& rDerivative,
								const double& a,
								const double& b) //dI1/dC
  {
    KRATOS_TRY

    rDerivative = msIdentityMatrix(a,b);

    return rDerivative;

    KRATOS_CATCH(" ")
  }

  double& HyperElasticModel::GetI2RightCauchyGreen1stDerivative(const CauchyGreenData& rData,
								double& rDerivative,
								const double& a,
								const double& b) //dI2/dC
  {
    KRATOS_TRY

    rDerivative = rData.Invariants.I1 * msIdentityMatrix(a,b) + rData.CauchyGreenMatrix(a,b);
	
    return rDerivative;

    KRATOS_CATCH(" ")
  }

  double& HyperElasticModel::GetI3RightCauchyGreen1stDerivative(const CauchyGreenData& rData,
								double& rDerivative,
								const double& a,
								const double& b) //dI3/dC
  {
    KRATOS_TRY

    rDerivative  = rData.Invariants.I3 * rData.InverseCauchyGreenMatrix(a,b);

    return rDerivative;

    KRATOS_CATCH(" ")
  }


  double& HyperElasticModel::GetJRightCauchyGreen1stDerivative(const CauchyGreenData& rData,
							       double& rDerivative,
							       const double& a,
							       const double& b) ///dJ/dC
  {
    KRATOS_TRY
	
    rDerivative = 0.5 * rData.Invariants.J * rData.InverseCauchyGreenMatrix(a,b);
	
    return rDerivative;

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************
  
  //Invariants Square of the 1st derivatives by components
  double& HyperElasticModel::GetI1RightCauchyGreenSquare1stDerivative(const CauchyGreenData& rData,
								      double& rDerivative,
								      const double& a,
								      const double& b,
								      const double& c,
								      const double& d) //dI1/dC * dI1/dC
  {
    KRATOS_TRY

    rDerivative  = msIdentityMatrix(a,b);
    rDerivative *= msIdentityMatrix(c,d);
      
    return rDerivative;

    KRATOS_CATCH(" ")
  }

  double& HyperElasticModel::GetI2RightCauchyGreenSquare1stDerivative(const CauchyGreenData& rData,
								      double& rDerivative,
								      const double& a,
								      const double& b,
								      const double& c,
								      const double& d) //dI2/dC * dI2/dC
  {
    KRATOS_TRY

    rDerivative  = rData.Invariants.I1 * msIdentityMatrix(a,b) + rData.CauchyGreenMatrix(a,b);
    rDerivative *= (rData.Invariants.I1 * msIdentityMatrix(c,d) + rData.CauchyGreenMatrix(c,d));
      
    return rDerivative;

    KRATOS_CATCH(" ")
  }

  double& HyperElasticModel::GetI3RightCauchyGreenSquare1stDerivative(const CauchyGreenData& rData,
								      double& rDerivative,
								      const double& a,
								      const double& b,
								      const double& c,
								      const double& d) //dI3/dC * dI3/dC
  {
    KRATOS_TRY

    rDerivative  = rData.Invariants.I3 * rData.InverseCauchyGreenMatrix(a,b);
    rDerivative *= (rData.Invariants.I3 * rData.InverseCauchyGreenMatrix(c,d));

    return rDerivative;

    KRATOS_CATCH(" ")
  }
    
  double& HyperElasticModel::GetJRightCauchyGreenSquare1stDerivative(const CauchyGreenData& rData,
								     double& rDerivative,
								     const double& a,
								     const double& b,
								     const double& c,
								     const double& d) //dJ/dC * dJ/dC
  {
    KRATOS_TRY
	
    rDerivative  = 0.5 * rData.Invariants.J * rData.InverseCauchyGreenMatrix(a,b);
    rDerivative *= (0.5 * rData.Invariants.J * rData.InverseCauchyGreenMatrix(c,d));
	
    return rDerivative;

    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************
    
  //Invariants 2nd derivatives by components
  double& HyperElasticModel::GetI1RightCauchyGreen2ndDerivative(const CauchyGreenData& rData,
								double& rDerivative,
								const double& a,
								const double& b,
								const double& c,
								const double& d) //ddI1/dCdC
  {
    KRATOS_TRY
      
    rDerivative = 0.0;
    
    return rDerivative;

    KRATOS_CATCH(" ")
  }

  double& HyperElasticModel::GetI2RightCauchyGreen2ndDerivative(const CauchyGreenData& rData,
								double& rDerivative,
								const double& a,
								const double& b,
								const double& c,
								const double& d) //ddI2/dCdC
  {
    KRATOS_TRY

    rDerivative  = GetFourthOrderUnitTensor(rDerivative,a,b,c,d);
    rDerivative *= rData.Invariants.I1;
    rDerivative -= rData.CauchyGreenMatrix(a,b)*msIdentityMatrix(c,d);	
      
    return rDerivative;

    KRATOS_CATCH(" ")
  }

  double& HyperElasticModel::GetI3RightCauchyGreen2ndDerivative(const CauchyGreenData& rData,
								double& rDerivative,
								const double& a,
								const double& b,
								const double& c,
								const double& d) //ddI3/dCdC
  {
    KRATOS_TRY

    rDerivative  = GetInverseRightCauchyGreenDerivative(rData,rDerivative,a,b,c,d);
    rDerivative += rData.InverseCauchyGreenMatrix(a,b)*rData.InverseCauchyGreenMatrix(c,d);
    rDerivative *= rData.Invariants.I3;
    return rDerivative;

    KRATOS_CATCH(" ")
  }

  double& HyperElasticModel::GetJRightCauchyGreen2ndDerivative(const CauchyGreenData& rData,
							       double& rDerivative,
							       const double& a,
							       const double& b,
							       const double& c,
							       const double& d) //ddJ/dCdC
  {
    KRATOS_TRY

    rDerivative  = GetInverseRightCauchyGreenDerivative(rData,rDerivative,a,b,c,d);
    rDerivative += 0.5 * rData.InverseCauchyGreenMatrix(a,b)*rData.InverseCauchyGreenMatrix(c,d);
    rDerivative *= 0.5 * rData.Invariants.J;
	
    return rDerivative;

    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************

    
  //************// left cauchy green : b
    
  HyperElasticModel::MatrixType& HyperElasticModel::GetI1LeftCauchyGreenDerivative(const CauchyGreenData& rData, MatrixType& rDerivative) //dI1/db
  {
    KRATOS_TRY
	
    noalias(rDerivative) = rData.CauchyGreenMatrix;

    return rDerivative;

    KRATOS_CATCH(" ")
  }

  HyperElasticModel::MatrixType& HyperElasticModel::GetI2LeftCauchyGreenDerivative(const CauchyGreenData& rData, MatrixType& rDerivative)  //dI2/db
  {
    KRATOS_TRY
	
    noalias(rDerivative)  = rData.CauchyGreenMatrix;
    rDerivative *= rData.Invariants.I1;
    noalias(rDerivative) -= prod(rData.CauchyGreenMatrix, rData.CauchyGreenMatrix);

      
    return rDerivative;

    KRATOS_CATCH(" ")
  }

  HyperElasticModel::MatrixType& HyperElasticModel::GetI3LeftCauchyGreenDerivative(const CauchyGreenData& rData, MatrixType& rDerivative) //dI3/db
  {
    KRATOS_TRY
      
    noalias(rDerivative) = msIdentityMatrix;
    rDerivative *= rData.Invariants.I3;
    
    return rDerivative;

    KRATOS_CATCH(" ")
  }


  //Invariants 1st derivatives by components
  double& HyperElasticModel::GetI1LeftCauchyGreen1stDerivative(const CauchyGreenData& rData,
							       double& rDerivative,
							       const double& a,
							       const double& b) //dI1/db
  {
    KRATOS_TRY

    rDerivative = rData.CauchyGreenMatrix(a,b);

    return rDerivative;

    KRATOS_CATCH(" ")
  }

  double& HyperElasticModel::GetI2LeftCauchyGreen1stDerivative(const CauchyGreenData& rData,
							       double& rDerivative,
							       const double& a,
							       const double& b) //dI2/db
  {
    KRATOS_TRY

    rDerivative  = rData.Invariants.I1 * rData.CauchyGreenMatrix(a,b) + rData.CauchyGreenMatrix(a,1)*rData.CauchyGreenMatrix(1,b) + rData.CauchyGreenMatrix(a,2)*rData.CauchyGreenMatrix(2,b) + rData.CauchyGreenMatrix(a,3)*rData.CauchyGreenMatrix(3,b);
	
    return rDerivative;

    KRATOS_CATCH(" ")
  }

  double& HyperElasticModel::GetI3LeftCauchyGreen1stDerivative(const CauchyGreenData& rData,
							       double& rDerivative,
							       const double& a,
							       const double& b) //dI3/db
  {
    KRATOS_TRY

    rDerivative  = rData.Invariants.I3 * msIdentityMatrix(a,b);

    return rDerivative;

    KRATOS_CATCH(" ")
  }

  double& HyperElasticModel::GetJRightCauchyGreenSquare1stDerivative(const CauchyGreenData& rData,
								     double& rDerivative,
								     const double& a,
								     const double& b) //dJ/db
  {
    KRATOS_TRY
	
    rDerivative  = 0.5 * rData.Invariants.J * msIdentityMatrix(a,b);
	
    return rDerivative;

    KRATOS_CATCH(" ")
  }
   
  //Invariants Square of the 1st derivatives by components
  double& HyperElasticModel::GetI1LeftCauchyGreenSquare1stDerivative(const CauchyGreenData& rData,
								     double& rDerivative,
								     const double& a,
								     const double& b,
								     const double& c,
								     const double& d) //dI1/db * dI1/db
  {
    KRATOS_TRY

    rDerivative  = rData.CauchyGreenMatrix(a,b);
    rDerivative *= rData.CauchyGreenMatrix(c,d);
      
    return rDerivative;

    KRATOS_CATCH(" ")
  }

  double& HyperElasticModel::GetI2LeftCauchyGreenSquare1stDerivative(const CauchyGreenData& rData,
								     double& rDerivative,
								     const double& a,
								     const double& b,
								     const double& c,
								     const double& d) //dI2/db * dI2/db
  {
    KRATOS_TRY

    rDerivative  = rData.Invariants.I1 * rData.CauchyGreenMatrix(a,b) + rData.CauchyGreenMatrix(a,1)*rData.CauchyGreenMatrix(1,b) + rData.CauchyGreenMatrix(a,2)*rData.CauchyGreenMatrix(2,b) + rData.CauchyGreenMatrix(a,3)*rData.CauchyGreenMatrix(3,b);
    rDerivative *= (rData.Invariants.I1 * rData.CauchyGreenMatrix(c,d) + rData.CauchyGreenMatrix(c,1)*rData.CauchyGreenMatrix(1,d) + rData.CauchyGreenMatrix(c,2)*rData.CauchyGreenMatrix(2,d) + rData.CauchyGreenMatrix(c,3)*rData.CauchyGreenMatrix(3,d));
      
    return rDerivative;

    KRATOS_CATCH(" ")
  }

  double& HyperElasticModel::GetI3LeftCauchyGreenSquare1stDerivative(const CauchyGreenData& rData,
								     double& rDerivative,
								     const double& a,
								     const double& b,
								     const double& c,
								     const double& d) //dI3/db * dI3/db
  {
    KRATOS_TRY

    rDerivative  =  rData.Invariants.I3 * msIdentityMatrix(a,b);
    rDerivative *= (rData.Invariants.I3 * msIdentityMatrix(c,d));

    return rDerivative;

    KRATOS_CATCH(" ")
  }

    
  double& HyperElasticModel::GetJLeftCauchyGreenSquare1stDerivative(const CauchyGreenData& rData,
								    double& rDerivative,
								    const double& a,
								    const double& b,
								    const double& c,
								    const double& d) //dJ/db * dJ/db
  {
    KRATOS_TRY

    rDerivative  = 0.5 * rData.Invariants.J * msIdentityMatrix(a,b);
    rDerivative *= (0.5 * rData.Invariants.J * msIdentityMatrix(c,d));
	
    return rDerivative;

    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************
  
  //Invariants 2nd derivatives by components
  double& HyperElasticModel::GetI1LeftCauchyGreen2ndDerivative(const CauchyGreenData& rData,
							       double& rDerivative,
							       const double& a,
							       const double& b,
							       const double& c,
							       const double& d) //ddI1/dbdb
  {
    KRATOS_TRY
      
    rDerivative = GetFourthOrderUnitTensor(rDerivative,a,b,c,d);
    
    return rDerivative;

    KRATOS_CATCH(" ")
  }

  double& HyperElasticModel::GetI2LeftCauchyGreen2ndDerivative(const CauchyGreenData& rData,
							       double& rDerivative,
							       const double& a,
							       const double& b,
							       const double& c,
							       const double& d) //ddI2/dbdb
  {
    KRATOS_TRY
    
    rDerivative  = GetFourthOrderUnitTensor(rDerivative,a,b,c,d);
    rDerivative *= rData.Invariants.I1;
    rDerivative += rData.CauchyGreenMatrix(a,b)*rData.CauchyGreenMatrix(c,d);
    rDerivative -= (msIdentityMatrix(a,c)*msIdentityMatrix(b,d)+msIdentityMatrix(a,d)*msIdentityMatrix(b,c)) * rData.CauchyGreenMatrix(b,d);      
    return rDerivative;

    KRATOS_CATCH(" ")
  }

  double& HyperElasticModel::GetI3LeftCauchyGreen2ndDerivative(const CauchyGreenData& rData,
							       double& rDerivative,
							       const double& a,
							       const double& b,
							       const double& c,
							       const double& d) //ddI3/dbdb
  {
    KRATOS_TRY

    rDerivative  = msIdentityMatrix(a,b)*msIdentityMatrix(c,d);
    rDerivative *= rData.Invariants.I3;
      
    return rDerivative;

    KRATOS_CATCH(" ")
  }

  double& HyperElasticModel::GetJLeftCauchyGreen2ndDerivative(const CauchyGreenData& rData,
							      double& rDerivative,
							      const double& a,
							      const double& b,
							      const double& c,
							      const double& d) //ddJ/dbdb
  {
    KRATOS_TRY

    rDerivative  = msIdentityMatrix(a,b)*msIdentityMatrix(c,d);
    rDerivative *= 0.25 * rData.Invariants.J;
	
    return rDerivative;

    KRATOS_CATCH(" ")
  }


  
  //************************************************************************************
  //************************************************************************************

  int HyperElasticModel::Check(const Properties& rMaterialProperties, const ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY
      
    if(HYPERELASTIC_MODEL_PARAMETERS.Key() == 0 || rMaterialProperties[HYPERELASTIC_MODEL_PARAMETERS].size() == 0)
      KRATOS_ERROR << "HYPERELASTIC_MODEL_PARAMETERS has Key zero or invalid value" << std::endl;

    // if(BULK_MODULUS.Key() == 0 || rMaterialProperties[BULK_MODULUS] <= 0.0)
    //   KRATOS_ERROR << "BULK_MODULUS has Key zero or invalid value" << std::endl;

    // if(LAME_LAMBDA.Key() == 0 || rMaterialProperties[LAME_LAMBDA] <= 0.0)
    //   KRATOS_ERROR << "LAME_LAMBDA has Key zero or invalid value" << std::endl;

    // if(LAME_MU.Key() == 0 || rMaterialProperties[LAME_MU] <= 0.0)
    //   KRATOS_ERROR << "LAME_MU has Key zero or invalid value" << std::endl;

    return 0;

	  
    KRATOS_CATCH(" ")
  }
  


} // Namespace Kratos
