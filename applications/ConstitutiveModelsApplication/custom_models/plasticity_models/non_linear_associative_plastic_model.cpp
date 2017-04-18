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
#include "custom_models/plasticity_models/non_linear_associative_plastic_model.hpp"

namespace Kratos
{

  //*****************************DEFAULT CONSTRUCTOR************************************
  //************************************************************************************

  template<class TElasticityModel,class TYieldCriterion>
  NonLinearAssociativePlasticModel<TElasticityModel,TYieldCriterion>::NonLinearAssociativePlasticModel()
    :BaseType()
  {
    KRATOS_TRY
	    

    KRATOS_CATCH(" ")    
  }
  
  //*******************************CONSTRUCTOR******************************************
  //************************************************************************************

  template<class TElasticityModel,class TYieldCriterion>
  NonLinearAssociativePlasticModel<TElasticityModel,TYieldCriterion>::NonLinearAssociativePlasticModel(ElasticityModelTypePointer pElasticityModel, YieldCriterionTypePointer pYieldCriterion)
    :BaseType(pElasticityModel, pYieldCriterion)
  {
    KRATOS_TRY
	    

    KRATOS_CATCH(" ")    
  }

  //*******************************ASSIGMENT OPERATOR***********************************
  //************************************************************************************

  template<class TElasticityModel,class TYieldCriterion>
  NonLinearAssociativePlasticModel<TElasticityModel,TYieldCriterion>& NonLinearAssociativePlasticModel<TElasticityModel,TYieldCriterion>::operator=(NonLinearAssociativePlasticModel const& rOther)
  {
    BaseType::operator=(rOther);
    return *this;
  }

  //*******************************COPY CONSTRUCTOR*************************************
  //************************************************************************************

  template<class TElasticityModel,class TYieldCriterion>
  NonLinearAssociativePlasticModel<TElasticityModel,TYieldCriterion>::NonLinearAssociativePlasticModel(NonLinearAssociativePlasticModel const& rOther)
    :BaseType(rOther)
    ,mInternal(rOther.mInternal)
    ,mPreviousInternal(rOther.mPreviousInternal)
    ,mThermalVariables(rOther.mThermalVariables)
  {

  }

  //********************************CLONE***********************************************
  //************************************************************************************

  template<class TElasticityModel,class TYieldCriterion>
  typename PlasticityModel<TElasticityModel,TYieldCriterion>::Pointer NonLinearAssociativePlasticModel<TElasticityModel,TYieldCriterion>::Clone() const
  {
    return ( BaseType::Pointer(new NonLinearAssociativePlasticModel(*this)) );
  }

  //********************************DESTRUCTOR******************************************
  //************************************************************************************

  template<class TElasticityModel,class TYieldCriterion>
  NonLinearAssociativePlasticModel<TElasticityModel,TYieldCriterion>::~NonLinearAssociativePlasticModel()
  {
  }

  //*********************************INITIALIZE VARIABLES*******************************
  //************************************************************************************

  template<class TElasticityModel,class TYieldCriterion>
  void NonLinearAssociativePlasticModel<TElasticityModel,TYieldCriterion>::InitializeVariables(ModelDataType& rValues, PlasticDataType& rVariables)
  {
    KRATOS_TRY

    //set model data pointer
    rVariables.SetModelData(rValues);

    rValues.State.Set(ConstitutiveModelData::PLASTIC_REGION,false);

    rValues.State.Set(ConstitutiveModelData::IMPLEX_ACTIVE,false);
    if( rValues.GetProcessInfo()[IMPLEX] == 1 )
      rValues.State.Set(ConstitutiveModelData::IMPLEX_ACTIVE,true);

    rVariables.SetState(rValues.State);
    
    // RateFactor
    rVariables.RateFactor = 0;

    // EquivalentPlasticStrain    
    rVariables.Internal = mInternal;

    // DeltaGamma / DeltaPlasticStrain (asociative plasticity)
    rVariables.DeltaInternal.Variables.clear();

    // Flow Rule local variables
    rVariables.TrialStateFunction = 0;
    rVariables.StressNorm = 0;
      
    KRATOS_CATCH(" ")
  }
	
  //*********************************CALCULATE STRESS***********************************
  //************************************************************************************

  template<class TElasticityModel,class TYieldCriterion>
  void NonLinearAssociativePlasticModel<TElasticityModel,TYieldCriterion>::CalculateStressTensor(ModelDataType& rValues, MatrixType& rStressMatrix)
  {

    KRATOS_TRY

    // calculate volumetric stress
    MatrixType VolumetricStressMatrix;
    VolumetricStressMatrix.clear();
    this->mpElasticityModel->CalculateVolumetricStressTensor(rValues,VolumetricStressMatrix);

    // calculate isochoric stress
    this->CalculateIsochoricStressTensor(rValues,rStressMatrix);

    rValues.StressMatrix = rStressMatrix;  //store isochoric stress as StressMatrix
      
    rStressMatrix += VolumetricStressMatrix;
	
    
    KRATOS_CATCH(" ")
    
  }

  //****************************CALCULATE ISOCHORIC STRESS******************************
  //************************************************************************************
  
  template<class TElasticityModel,class TYieldCriterion>
  void NonLinearAssociativePlasticModel<TElasticityModel,TYieldCriterion>::CalculateIsochoricStressTensor(ModelDataType& rValues, MatrixType& rStressMatrix)
  {

    KRATOS_TRY
    
    PlasticDataType Variables;
    this->InitializeVariables(rValues,Variables);

    // calculate elastic isochoric stress
    this->mpElasticityModel->CalculateIsochoricStressTensor(rValues,rStressMatrix);

    // calculate plastic isochoric stress
    this->CalculateAndAddIsochoricStressTensor(Variables,rStressMatrix);


    if( rValues.State.Is(ConstitutiveModelData::UPDATE_INTERNAL_VARIABLES ) )
      this->UpdateInternalVariables(Variables);

    
    KRATOS_CATCH(" ")
    
  }

  //****************************CALCULATE ISOCHORIC STRESS******************************
  //************************************************************************************
  
  template<class TElasticityModel,class TYieldCriterion>
  void NonLinearAssociativePlasticModel<TElasticityModel,TYieldCriterion>::CalculateAndAddIsochoricStressTensor(PlasticDataType& rVariables, MatrixType& rStressMatrix)
  {

    KRATOS_TRY

    //1.-Isochoric stress norm   (rStressMatrix := Elastic Isochoric Stress Matrix)
    rVariables.StressNorm = ConstitutiveLawUtilities::CalculateStressNorm(rStressMatrix, rVariables.StressNorm);

    //2.-Check yield condition
    rVariables.TrialStateFunction = this->mpYieldCriterion->CalculateYieldCondition(rVariables, rVariables.TrialStateFunction);
    
    if( rVariables.State().Is(ConstitutiveModelData::IMPLEX_ACTIVE) ) 
      {
	//3.- Calculate the implex radial return
	this->CalculateImplexRadialReturn(rVariables,rStressMatrix);
      }
    else{

      if( rVariables.TrialStateFunction <= 0 )
	{
	  rVariables.State().Set(ConstitutiveModelData::PLASTIC_REGION,false);
	}
      else
	{

	  //3.- Calculate the radial return
	  bool converged = this->CalculateRadialReturn(rVariables,rStressMatrix);
	    
	  if(!converged)
	    std::cout<<" ConstitutiveLaw did not converge "<<std::endl;


	  //4.- Update back stress, plastic strain and stress
	  this->UpdateStressConfiguration(rVariables,rStressMatrix);
	    
	  //5.- Calculate thermal dissipation and delta thermal dissipation
	  this->CalculateThermalDissipation(rVariables);   

	  rVariables.State().Set(ConstitutiveModelData::PLASTIC_REGION,true);	    
	}

    }

    rVariables.State().Set(ConstitutiveModelData::COMPUTED_RETURN_MAPPING,true);
   
    KRATOS_CATCH(" ")
    
  }

  //***************************CALCULATE RADIAL RETURN MAPPING**************************
  //************************************************************************************
  
  // local newton procedure
  template<class TElasticityModel,class TYieldCriterion>
  bool NonLinearAssociativePlasticModel<TElasticityModel,TYieldCriterion>::CalculateRadialReturn(PlasticDataType& rVariables, MatrixType& rStressMatrix)
  {
    KRATOS_TRY
      
    //Set convergence parameters
    unsigned int iter    = 0;
    double Tolerance     = 1e-5;
    double MaxIterations = 50;

    //start
    double DeltaDeltaGamma    = 0;
    double DeltaStateFunction = 0;
    double DeltaPlasticStrain = 0;
    
    double& EquivalentPlasticStrainOld  = mPreviousInternal.Variables[0];
    double& EquivalentPlasticStrain     = rVariables.Internal.Variables[0];
    double& DeltaGamma                  = rVariables.DeltaInternal.Variables[0];

    EquivalentPlasticStrain = 0;
    DeltaGamma = 0;

    double StateFunction                = rVariables.TrialStateFunction;
    
    while ( fabs(StateFunction)>=Tolerance && iter<=MaxIterations)
      {
	//Calculate Delta State Function:
	DeltaStateFunction = this->mpYieldCriterion->CalculateDeltaStateFunction( DeltaStateFunction );

	//Calculate DeltaGamma:
	DeltaDeltaGamma  = StateFunction/DeltaStateFunction;
	DeltaGamma += DeltaDeltaGamma;
	       
	//Update Equivalent Plastic Strain:
	DeltaPlasticStrain      = sqrt(2.0/3.0) * DeltaGamma;
	EquivalentPlasticStrain = EquivalentPlasticStrainOld + DeltaPlasticStrain;
	       	
	//Calculate State Function:
	StateFunction = this->mpYieldCriterion->CalculateStateFunction( StateFunction );

	iter++;
      }
	   

    if(iter>MaxIterations)
      return false;


    return true;

    KRATOS_CATCH(" ")
  }

  //***************************CALCULATE IMPLEX RETURN MAPPING**************************
  //************************************************************************************

  template<class TElasticityModel,class TYieldCriterion>
  void NonLinearAssociativePlasticModel<TElasticityModel,TYieldCriterion>::CalculateImplexRadialReturn(PlasticDataType& rVariables, MatrixType& rStressMatrix)
  {
    KRATOS_TRY
      
    double& EquivalentPlasticStrainOld  = mPreviousInternal.Variables[0];
    double& EquivalentPlasticStrain     = rVariables.Internal.Variables[0];
    double& DeltaGamma                  = rVariables.DeltaInternal.Variables[0];
   
    //1.-Computation of the plastic Multiplier
    DeltaGamma = sqrt(3.0/2.0) * ( EquivalentPlasticStrain - EquivalentPlasticStrainOld );
	
    //2.- Update back stress, plastic strain and stress
    this->UpdateStressConfiguration(rVariables,rStressMatrix);

    //3.- Calculate thermal dissipation and delta thermal dissipation
    if( DeltaGamma > 0 ){
	  
      this->CalculateImplexThermalDissipation( rVariables );
      rVariables.State().Set(ConstitutiveModelData::PLASTIC_REGION,true);
      
    }
    else{
      
      mThermalVariables.PlasticDissipation = 0;
      mThermalVariables.DeltaPlasticDissipation = 0;
    } 
  
    KRATOS_CATCH(" ")
  }


  //***************************CALCULATE THERMAL DISSIPATION****************************
  //************************************************************************************


  template<class TElasticityModel,class TYieldCriterion>
  void NonLinearAssociativePlasticModel<TElasticityModel,TYieldCriterion>::CalculateThermalDissipation( PlasticDataType& rVariables )
  {
    KRATOS_TRY
      
    //1.- Thermal Dissipation:
 
    mThermalVariables.PlasticDissipation = this->mpYieldCriterion->CalculatePlasticDissipation( mThermalVariables.PlasticDissipation );
  

    //std::cout<<" PlasticDissipation "<<mThermalVariables.PlasticDissipation<<std::endl;

    //2.- Thermal Dissipation Increment:

    mThermalVariables.DeltaPlasticDissipation = this->mpYieldCriterion->CalculateDeltaPlasticDissipation( mThermalVariables.DeltaPlasticDissipation );
		    		    
    //std::cout<<" DeltaPlasticDissipation "<<mThermalVariables.DeltaPlasticDissipation<<std::endl;

    KRATOS_CATCH(" ")    
  }


  //***************************CALCULATE THERMAL DISSIPATION****************************
  //************************************************************************************

  template<class TElasticityModel,class TYieldCriterion>
  void NonLinearAssociativePlasticModel<TElasticityModel,TYieldCriterion>::CalculateImplexThermalDissipation( PlasticDataType& rVariables )
  {
    KRATOS_TRY
       
    //1.- Thermal Dissipation:
	
    mThermalVariables.PlasticDissipation = this->mpYieldCriterion->CalculateImplexPlasticDissipation( mThermalVariables.PlasticDissipation );
  
    //2.- Thermal Dissipation Increment:
      
    mThermalVariables.DeltaPlasticDissipation = this->mpYieldCriterion->CalculateImplexDeltaPlasticDissipation( mThermalVariables.DeltaPlasticDissipation );

    KRATOS_CATCH(" ")  
  
  }

  //***************************UPDATE STRESS CONFIGURATION *****************************
  //************************************************************************************

  template<class TElasticityModel,class TYieldCriterion>
  void NonLinearAssociativePlasticModel<TElasticityModel,TYieldCriterion>::UpdateStressConfiguration(PlasticDataType& rVariables, MatrixType& rStressMatrix)
  {
    KRATOS_TRY
          
    //Plastic Strain and Back Stress update
    if( rVariables.StressNorm > 0 ){

      double& DeltaGamma = rVariables.DeltaInternal.Variables[0];
      const MaterialDataType& rMaterial = rVariables.GetMaterialParameters();
             
      //Stress Update: 
      rStressMatrix -= ( rStressMatrix * ( ( 2.0 * rMaterial.GetLameMuBar() * DeltaGamma ) / rVariables.StressNorm ) );
	  
    }

    KRATOS_CATCH(" ")    
  }

  //***************************UPDATE INTERNAL VARIABLES********************************
  //************************************************************************************

  template<class TElasticityModel,class TYieldCriterion>
  void NonLinearAssociativePlasticModel<TElasticityModel,TYieldCriterion>::UpdateInternalVariables(PlasticDataType& rVariables)
  {
    KRATOS_TRY
    
    double& EquivalentPlasticStrainOld  = mPreviousInternal.Variables[0];
    double& EquivalentPlasticStrain     = rVariables.Internal.Variables[0];
    double& DeltaGamma                  = rVariables.DeltaInternal.Variables[0];

    //update mechanical variables
    EquivalentPlasticStrainOld  = EquivalentPlasticStrain;
    EquivalentPlasticStrain    += sqrt(2.0/3.0) * DeltaGamma;
	
    //update thermal variables
    //mThermalVariables = rVariables.Thermal;
    
    KRATOS_CATCH(" ")    
  }


  //**************CALCULATE SCALING FACTORS FOR THE ELASTO PLASTIC MODULI***************
  //************************************************************************************

  template<class TElasticityModel,class TYieldCriterion>
  void NonLinearAssociativePlasticModel<TElasticityModel,TYieldCriterion>::CalculateScalingFactors(PlasticDataType& rVariables, PlasticFactors& rFactors)
  {
    KRATOS_TRY

    //0.-Get needed parameters
    const ModelDataType&    rModelData = rVariables.GetModelData();
    const MaterialDataType& rMaterial  = rVariables.GetMaterialParameters();

    double&             rDeltaGamma            = rVariables.DeltaInternal.Variables[0];
    const MatrixType&   rIsochoricStressMatrix = rModelData.GetStressMatrix(); //isochoric stress stored as StressMatrix
    
    //1.-Identity build
    MatrixType Identity = identity_matrix<double> (3);
    
    //2.-Auxiliar matrices
    rFactors.Normal      = rIsochoricStressMatrix * ( 1.0 / rVariables.StressNorm );

    MatrixType Norm_Normal      = prod( rFactors.Normal, trans(rFactors.Normal) );

    double Trace_Norm_Normal    = Norm_Normal( 0, 0 ) + Norm_Normal( 1, 1 ) + Norm_Normal( 2, 2 );

    rFactors.Dev_Normal  = Norm_Normal;
    
    rFactors.Dev_Normal -= (1.0/3.0) * Trace_Norm_Normal * Identity;

    //3.-Auxiliar constants    
    if( rVariables.State().Is(ConstitutiveModelData::IMPLEX_ACTIVE) ) 
      {
	
	rFactors.Beta0 = 0;
		
	rFactors.Beta1 = 2.0 * rMaterial.GetLameMuBar() * rDeltaGamma / rVariables.StressNorm;
		
	rFactors.Beta2 = (2.0/3.0) * rVariables.StressNorm * rDeltaGamma / ( rMaterial.GetLameMuBar() );
		
	rFactors.Beta3 = ( -rFactors.Beta1 + rFactors.Beta2 );
		
	rFactors.Beta4 = ( -rFactors.Beta1 ) * rVariables.StressNorm / ( rMaterial.GetLameMuBar() );

      }
    else
      {
    
	if( rVariables.State().Is(ConstitutiveModelData::PLASTIC_RATE_REGION) )
	  rVariables.RateFactor = 1;
	else if ( rVariables.State().IsNot(ConstitutiveModelData::PLASTIC_RATE_REGION) )
	  rVariables.RateFactor = 0;
	
	double DeltaHardening = this->mpYieldCriterion->GetHardeningLaw().CalculateDeltaHardening( rVariables, DeltaHardening );

	rFactors.Beta0 = 1.0 + DeltaHardening/(3.0 * rMaterial.GetLameMuBar());
		
	rFactors.Beta1 = 2.0 * rMaterial.GetLameMuBar() * rDeltaGamma / rVariables.StressNorm;
		
	rFactors.Beta2 = ( ( 1.0 - ( 1.0 / rFactors.Beta0 ) ) * (2.0/3.0) * rVariables.StressNorm * rDeltaGamma )/( rMaterial.GetLameMuBar() );
		
	rFactors.Beta3 = ( ( 1.0 / rFactors.Beta0 ) - rFactors.Beta1 + rFactors.Beta2 );
		
	rFactors.Beta4 = ( ( 1.0 / rFactors.Beta0 ) - rFactors.Beta1 ) * rVariables.StressNorm / ( rMaterial.GetLameMuBar() );

      }
	
    //std::cout<<"FACTORS:: Beta0 "<<rFactors.Beta0<<" Beta 1 "<<rFactors.Beta1<<" Beta2 "<<rFactors.Beta2<<" Beta 3 "<<rFactors.Beta3<<" Beta4 "<<rFactors.Beta4<<std::endl;

    KRATOS_CATCH(" ")    
	  
  }

  //***********************COMPUTE ALGORITHMIC CONSTITUTIVE MATRIX**********************
  //************************************************************************************

  template<class TElasticityModel,class TYieldCriterion>
  void NonLinearAssociativePlasticModel<TElasticityModel,TYieldCriterion>::CalculateConstitutiveTensor(ModelDataType& rValues, Matrix& rConstitutiveMatrix)
  {
    KRATOS_TRY

    //Initialize ConstitutiveMatrix
    rConstitutiveMatrix.clear();
      
    PlasticDataType Variables;
    this->InitializeVariables(rValues,Variables);

    // calculate isochoric stress (radial return is needed)

    MatrixType StressMatrix = rValues.StressMatrix;
    //1.-Elastic Isochoric Stress Matrix
    this->mpElasticityModel->CalculateIsochoricStressTensor(rValues,rValues.StressMatrix);

    //2.-Calculate and Add Plastic Isochoric Stress Matrix
    this->CalculateAndAddIsochoricStressTensor(Variables,rValues.StressMatrix);
    
    //Calculate Constitutive Matrix

    // calculate elastic constitutive tensor
    this->mpElasticityModel->CalculateConstitutiveTensor(rValues,rConstitutiveMatrix);
    
    // calculate plastic constitutive tensor
    this->CalculateAndAddPlasticConstitutiveTensor(Variables,rConstitutiveMatrix);

    rValues.StressMatrix = StressMatrix;
    
    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************
  
  template<class TElasticityModel,class TYieldCriterion>
  void NonLinearAssociativePlasticModel<TElasticityModel,TYieldCriterion>::CalculateStressAndConstitutiveTensors(ModelDataType& rValues, MatrixType& rStressMatrix, Matrix& rConstitutiveMatrix)
  {
    KRATOS_TRY

    PlasticDataType Variables;
    this->InitializeVariables(rValues,Variables);

    //Calculate Stress Matrix

    // calculate volumetric stress
    MatrixType VolumetricStressMatrix;
    this->mpElasticityModel->CalculateVolumetricStressTensor(rValues,VolumetricStressMatrix);

    // calculate isochoric stress

    // calculate elastic isochoric stress
    this->mpElasticityModel->CalculateIsochoricStressTensor(rValues,rStressMatrix);
   
    // calculate plastic isochoric stress
    this->CalculateAndAddIsochoricStressTensor(Variables,rStressMatrix);

     rValues.StressMatrix = rStressMatrix;  //store isochoric stress as StressMatrix
   
    //Calculate Constitutive Matrix
    
    // calculate elastic constitutive tensor
    this->mpElasticityModel->CalculateConstitutiveTensor(rValues,rConstitutiveMatrix);
    
    // calculate plastic constitutive tensor
    this->CalculateAndAddPlasticConstitutiveTensor(Variables,rConstitutiveMatrix);

    
    rStressMatrix += VolumetricStressMatrix;
    
    
    KRATOS_CATCH(" ")
  }

  
  //************************************************************************************
  //************************************************************************************
  
  template<class TElasticityModel,class TYieldCriterion>
  void NonLinearAssociativePlasticModel<TElasticityModel,TYieldCriterion>::CalculateAndAddPlasticConstitutiveTensor(PlasticDataType& rVariables, Matrix& rConstitutiveMatrix)
  {
    KRATOS_TRY
    
    //Compute radial return
    if( rVariables.State().IsNot(ConstitutiveModelData::COMPUTED_RETURN_MAPPING) )
      KRATOS_THROW_ERROR( std::logic_error, "ReturnMapping has to be computed to perform the calculation", "" )
    
    //Algorithmic moduli factors
    PlasticFactors Factors;    
    this->CalculateScalingFactors(rVariables,Factors);
    
    //Calculate HyperElastic ConstitutiveMatrix
    const ModelDataType&  rModelData        = rVariables.GetModelData();
    const SizeType&       rVoigtSize        = rModelData.GetVoigtSize();      
    const VoigtIndexType& rIndexVoigtTensor = rModelData.GetVoigtIndexTensor();
    
    for(SizeType i=0; i<rVoigtSize; i++)
      {
	for(SizeType j=0; j<rVoigtSize; j++)
	  {
	    
	    rConstitutiveMatrix(i,j) = this->AddPlasticConstitutiveComponent(rVariables,Factors,rConstitutiveMatrix(i,j),
									     rIndexVoigtTensor[i][0],rIndexVoigtTensor[i][1],
									     rIndexVoigtTensor[j][0],rIndexVoigtTensor[j][1]);
	  }
	
      }

    rVariables.State().Set(ConstitutiveModelData::COMPUTED_CONSTITUTIVE_MATRIX);
    
    KRATOS_CATCH(" ")
  }
  


  //********************CONSTITUTIVE MATRIX PLASTIC COMPONENT***************************
  //************************************************************************************

  template<class TElasticityModel,class TYieldCriterion>
  double& NonLinearAssociativePlasticModel<TElasticityModel,TYieldCriterion>::AddPlasticConstitutiveComponent(PlasticDataType& rVariables,
													      PlasticFactors& rFactors, double& rCabcd,
													      const unsigned int& a, const unsigned int& b,
													      const unsigned int& c, const unsigned int& d)
  {
    KRATOS_TRY

    const ModelDataType&    rModelData       = rVariables.GetModelData();
    const MaterialDataType& rMaterial        = rVariables.GetMaterialParameters();

    const MatrixType& rIsochoricStressMatrix = rModelData.GetStressMatrix(); //isochoric stress stored as StressMatrix
    const MatrixType& rCauchyGreenMatrix     = rModelData.GetStrainMatrix();

    
    double Cabcd = (1.0/3.0) * ( rCauchyGreenMatrix(a,b) * rCauchyGreenMatrix(c,d) );
    
    Cabcd -= (0.5 * ( rCauchyGreenMatrix(a,c) * rCauchyGreenMatrix(b,d) + rCauchyGreenMatrix(a,d) * rCauchyGreenMatrix(b,c) ) );
    
    Cabcd *= 3.0 * rMaterial.GetLameMuBar();

    Cabcd += ( rCauchyGreenMatrix(c,d)* rIsochoricStressMatrix(a,b) + rIsochoricStressMatrix(c,d) * rCauchyGreenMatrix(a,b) );
 
    Cabcd *= (-2.0/3.0) * ( (-1) * rFactors.Beta1 );
    
    Cabcd -= rFactors.Beta3 * 2.0 * rMaterial.GetLameMuBar() * ( rFactors.Normal(a,b) * rFactors.Normal(c,d) );

    Cabcd -= rFactors.Beta4 * 2.0 * rMaterial.GetLameMuBar() * ( rFactors.Normal(a,b) * rFactors.Dev_Normal(c,d) );

    rCabcd += Cabcd;
    
    return rCabcd;     
    
    KRATOS_CATCH(" ")    
  }


}  // namespace Kratos.
