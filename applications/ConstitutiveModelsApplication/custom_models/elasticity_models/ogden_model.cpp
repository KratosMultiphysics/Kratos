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
#include "custom_models/elasticity_models/ogden_model.hpp"


namespace Kratos
{

  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

  OgdenModel::OgdenModel()
    : HyperElasticModel()
  {
  }

  //******************************COPY CONSTRUCTOR**************************************
  //************************************************************************************

  OgdenModel::OgdenModel(const OgdenModel& rOther)
    : HyperElasticModel(rOther)
  {
  }

  //********************************CLONE***********************************************
  //************************************************************************************

  ConstitutiveModel::Pointer OgdenModel::Clone() const
  {
    return ( OgdenModel::Pointer(new OgdenModel(*this)) );
  }

  //********************************ASSIGNMENT******************************************
  //************************************************************************************
  OgdenModel& OgdenModel::operator=(OgdenModel const& rOther)
  {
    HyperElasticModel::operator=(rOther);
    return *this;
  }

  //*******************************DESTRUCTOR*******************************************
  //************************************************************************************

  OgdenModel::~OgdenModel()
  {
  }

  //************************************************************************************
  //************************************************************************************

  void OgdenModel::CalculateStrainEnergy(ModelDataType& rValues, double& rDensityFunction)
  {
      KRATOS_TRY
	
      HyperElasticDataType Variables;
      this->CalculateStrainData(rValues, Variables);

      this->CalculateAndAddVolumetricStrainEnergy( Variables, rDensityFunction );

      
      const MaterialDataType& rMaterial = Variables.GetMaterialParameters();

      const std::vector<double>& rModelParameters = rMaterial.GetModelParameters(); //nu values, lambda values

      unsigned int size = (rModelParameters.size()/2.0);	

      for(unsigned int p=0; p<size; p++)
      {
	  const double& mu_p = rModelParameters[p];
	  const double& alpha_p = rModelParameters[p+size];
	  rDensityFunction += (mu_p/alpha_p) * ( std::pow(Variables.Strain.Eigen.Values[0],alpha_p) + std::pow(Variables.Strain.Eigen.Values[1],alpha_p) + std::pow(Variables.Strain.Eigen.Values[2],alpha_p) - 3.0 );
      }
      
      KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************
  
  void OgdenModel::CalculateConstitutiveTensor(ModelDataType& rValues, Matrix& rConstitutiveMatrix)
  {
    KRATOS_TRY

    //Initialize ConstitutiveMatrix  
    HyperElasticDataType Variables;
    this->CalculateStrainData(rValues,Variables);

    //Calculate Stress
    MatrixType StressMatrix;
    noalias(StressMatrix) = ZeroMatrix(3,3);
    this->CalculateAndAddStressTensor(Variables,StressMatrix);
    rValues.StressMatrix = StressMatrix; //store total stress as StressMatrix

    
    //Calculate Constitutive Matrix
    this->CalculateAndAddConstitutiveTensor(Variables,rConstitutiveMatrix);   
    	
    KRATOS_CATCH(" ")
  }


  //***********************PROTECTED OPERATIONS FROM BASE CLASS*************************
  //************************************************************************************
    
  //************// W

  //option: g(J) = (K/2)*ln(J)²  (default implemented)
  void OgdenModel::CalculateAndAddVolumetricStrainEnergy(HyperElasticDataType& rVariables, double& rVolumetricDensityFunction)
  {
      KRATOS_TRY

      const MaterialDataType& rMaterial = rVariables.GetMaterialParameters();

      //g(J) = (K/2)*ln(J)² 
      rVolumetricDensityFunction += rMaterial.GetBulkModulus() * 0.5 * std::pow(std::log(rVariables.Strain.Invariants.J),2);
	
      KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************
    
  void OgdenModel::CalculateAndAddStressTensor(HyperElasticDataType& rVariables, MatrixType& rStressMatrix)
  {
      KRATOS_TRY

      const ModelDataType&  rModelData        = rVariables.GetModelData();
      const StressMeasureType& rStressMeasure = rModelData.GetStressMeasure();

      array_1d<double,3> MainStresses;
      this->CalculateMainStresses(rVariables,MainStresses);
    
      if( rStressMeasure == ConstitutiveModelData::StressMeasure_PK2 ){ //Strain.Matrix = RightCauchyGreen (C)
	  array_1d<double,3> EigenVector;

	  for(unsigned int i=0; i<3; i++)
	  {
	      noalias(EigenVector) = matrix_row<const MatrixType>(rVariables.Strain.Eigen.Vectors,i);
	      EigenVector /= rVariables.Strain.Eigen.Values[i];
	      noalias(rStressMatrix) += MainStresses[i] * outer_prod(EigenVector,EigenVector);
	  }

      }
      else if( rStressMeasure == ConstitutiveModelData::StressMeasure_Kirchhoff ){ //Strain.Matrix = LeftCauchyGreen (b)
	  
	  array_1d<double,3> EigenVector;
	  for(unsigned int i=0; i<3; i++)
	  {
	      noalias(EigenVector) = matrix_row<const MatrixType>(rVariables.Strain.Eigen.Vectors,i);
	      noalias(rStressMatrix) += MainStresses[i] * outer_prod(EigenVector,EigenVector);
	  }

      }

      rVariables.State().Set(ConstitutiveModelData::STRESS_COMPUTED);

      KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************
  
    
  void OgdenModel::CalculateMainStresses(HyperElasticDataType& rVariables, array_1d<double,3>& rMainStresses)
  {
      KRATOS_TRY
	    
 
      const MaterialDataType& rMaterial = rVariables.GetMaterialParameters();
      const std::vector<double>& rModelParameters = rMaterial.GetModelParameters(); //nu values, lambda values
      const double& rBulkModulus = rMaterial.GetBulkModulus();
      unsigned int size = (rModelParameters.size()/2.0);
      double athird = 1.0/3.0;
	
      for(unsigned int i=0; i<3; i++)
      {
	  for(unsigned int p=0; p<size; p++)
	  {
	      const double& mu_p = rModelParameters[p];
	      const double& alpha_p = rModelParameters[p+size];
	      rMainStresses[i] += (mu_p) * std::pow(rVariables.Strain.Invariants.J, (-alpha_p*athird)) * ( std::pow(rVariables.Strain.Eigen.Values[i],alpha_p) - athird * ( std::pow(rVariables.Strain.Eigen.Values[0],alpha_p) + std::pow(rVariables.Strain.Eigen.Values[1],alpha_p) + std::pow(rVariables.Strain.Eigen.Values[2],alpha_p) ) );
	      rMainStresses[i] += rBulkModulus * std::log(rVariables.Strain.Invariants.J);
	  }	    
	    
      }
	    
      KRATOS_CATCH(" ")
   }

    
  //************************************************************************************
  //************************************************************************************
  
    
  void OgdenModel::CalculateStrainData(ModelDataType& rValues, HyperElasticDataType& rVariables)
  {
    KRATOS_TRY

    //set model data pointer
    rVariables.SetModelData(rValues);
    rVariables.SetState(rValues.State);

    //deformation gradient
    const MatrixType& rDeltaDeformationMatrix = rValues.GetDeltaDeformationMatrix();
    const MatrixType& rTotalDeformationMatrix = rValues.GetTotalDeformationMatrix();

    const StressMeasureType& rStressMeasure = rValues.GetStressMeasure();

    if( rStressMeasure == ConstitutiveModelData::StressMeasure_PK2 ){ //mStrainMatrix = RightCauchyGreen (C=FT*F)  C^-1=(FT*F)^-1=F^-1*FT^-1

	//set working strain measure
	rValues.SetStrainMeasure(ConstitutiveModelData::CauchyGreen_Right);
      
	//historical strain matrix
	rValues.StrainMatrix = ConstitutiveModelUtilities::VectorToSymmetricTensor(this->mHistoryVector,rValues.StrainMatrix);

	//current strain matrix b
	noalias(rVariables.Strain.Matrix) = prod(rValues.StrainMatrix,trans(rDeltaDeformationMatrix));
	noalias(rValues.StrainMatrix) = prod(rDeltaDeformationMatrix, rVariables.Strain.Matrix);

	//inverted total deformation gradient
	ConstitutiveModelUtilities::InvertMatrix3( rTotalDeformationMatrix, rVariables.Strain.InverseMatrix, rVariables.Strain.Invariants.I3 ); //InverseMatrix and I3 is used as wildcard here (InverseMatrix = InverseTotalDeformationGradient)

	//strain measure C
	noalias(rVariables.Strain.Matrix) = prod(rValues.StrainMatrix,trans(rVariables.Strain.InverseMatrix));
	rVariables.Strain.Matrix = prod(trans(rTotalDeformationMatrix), rVariables.Strain.Matrix);

	//inverted strain measure
	ConstitutiveModelUtilities::InvertMatrix3( rVariables.Strain.Matrix, rVariables.Strain.InverseMatrix, rVariables.Strain.Invariants.I3 );

	MatrixType StrainMatrix;
	noalias(StrainMatrix) = rVariables.Strain.Matrix;
	rVariables.Strain.Eigen.Vectors.clear();

	MathUtils<double>::EigenSystem<3> ( StrainMatrix, rVariables.Strain.Eigen.Vectors, rVariables.Strain.Matrix);

	for (unsigned int i = 0; i < 3; i++)
	    rVariables.Strain.Eigen.Values[i] = std::sqrt(rVariables.Strain.Matrix(i,i));

	rVariables.Strain.Matrix = prod( rVariables.Strain.Matrix, rVariables.Strain.Eigen.Vectors);
	rVariables.Strain.Matrix = prod( trans(rVariables.Strain.Eigen.Vectors), rVariables.Strain.Matrix);
	
	rValues.State.Set(ConstitutiveModelData::STRAIN_COMPUTED);
	
    }
    else if( rStressMeasure == ConstitutiveModelData::StressMeasure_Kirchhoff ){ //Strain.Matrix = LeftCauchyGreen (b=F*FT)

      //set working strain measure
      rValues.SetStrainMeasure(ConstitutiveModelData::CauchyGreen_Left);

      //historical strain matrix
      rValues.StrainMatrix = ConstitutiveModelUtilities::VectorToSymmetricTensor(this->mHistoryVector,rValues.StrainMatrix);

      //current strain matrix b
      noalias(rVariables.Strain.Matrix) = prod(rValues.StrainMatrix,trans(rDeltaDeformationMatrix));
      noalias(rValues.StrainMatrix) = prod(rDeltaDeformationMatrix, rVariables.Strain.Matrix);

      rVariables.Strain.Matrix.clear();
      rVariables.Strain.Eigen.Vectors.clear();

      MathUtils<double>::EigenSystem<3> ( rValues.StrainMatrix, rVariables.Strain.Eigen.Vectors, rVariables.Strain.Matrix);

      for (unsigned int i = 0; i < 3; i++)
	  rVariables.Strain.Eigen.Values[i] = std::sqrt(rVariables.Strain.Matrix(i,i));

      rVariables.Strain.Matrix = prod( rVariables.Strain.Matrix, rVariables.Strain.Eigen.Vectors);
      rVariables.Strain.Matrix = prod( trans(rVariables.Strain.Eigen.Vectors), rVariables.Strain.Matrix);

      noalias(rVariables.Strain.Matrix) = rValues.StrainMatrix;

      //inverted strain measure
      ConstitutiveModelUtilities::InvertMatrix3( rVariables.Strain.Matrix, rVariables.Strain.InverseMatrix, rVariables.Strain.Invariants.I3 );
      rValues.State.Set(ConstitutiveModelData::STRAIN_COMPUTED);

    }
    else{

      //set working strain measure
      rValues.SetStrainMeasure(ConstitutiveModelData::CauchyGreen_None);
      KRATOS_ERROR << "calling initialize OgdenModel .. StressMeasure is inconsistent"  << std::endl;

    }

    //Calculate Jacobian   
    rVariables.Strain.Invariants.J = rVariables.GetModelData().GetTotalDeformationDet();
    

    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************
       
  void OgdenModel::CalculateMainStressDerivatives(HyperElasticDataType& rVariables, MatrixType& rStressDerivatives)
  {
    KRATOS_TRY
              
    //Calculate Ogden main stress derivatives
    const ModelDataType&  rModelData        = rVariables.GetModelData();
    const StressMeasureType& rStressMeasure = rModelData.GetStressMeasure();
	
    const MaterialDataType& rMaterial = rVariables.GetMaterialParameters();

    const double& rBulkModulus = rMaterial.GetBulkModulus();
    const std::vector<double>& rModelParameters = rMaterial.GetModelParameters(); //nu values, lambda values

    unsigned int size = (rModelParameters.size()/2.0);
    double athird = 1.0/3.0;

    //Calculate Stress main streches derivatives
    if( rStressMeasure == ConstitutiveModelData::StressMeasure_PK2 ){ //Strain.Matrix = LeftCauchyGreen (C)
	for(SizeType i=0; i<3; i++)
	{
	    for(SizeType j=0; j<3; j++)
	    {
		for(unsigned int p=0; p<size; p++)
		{
		    const double& mu_p = rModelParameters[p];
		    const double& alpha_p = rModelParameters[p+size];
		    double f = athird * ( std::pow(rVariables.Strain.Eigen.Values[0],alpha_p) + std::pow(rVariables.Strain.Eigen.Values[1],alpha_p) + std::pow(rVariables.Strain.Eigen.Values[2],alpha_p) );

		    rStressDerivatives(i,j) += (mu_p * alpha_p * std::pow(rVariables.Strain.Invariants.J,(-alpha_p*athird)) * ( f - std::pow(rVariables.Strain.Eigen.Values[i],alpha_p) - std::pow(rVariables.Strain.Eigen.Values[j],alpha_p) + 3.0 * std::pow(rVariables.Strain.Eigen.Values[i],alpha_p) * this->msIdentityMatrix(i,j) ) / 6.0 ) / (rVariables.Strain.Eigen.Values[j] * rVariables.Strain.Eigen.Values[j] * rVariables.Strain.Eigen.Values[i] * rVariables.Strain.Eigen.Values[i]);
		
		    rStressDerivatives(i,j) -= (mu_p * std::pow(rVariables.Strain.Invariants.J,(-alpha_p*athird)) * ( std::pow(rVariables.Strain.Eigen.Values[i],alpha_p) - f )) * 2.0  / (rVariables.Strain.Eigen.Values[j] * rVariables.Strain.Eigen.Values[i] * rVariables.Strain.Eigen.Values[i] * rVariables.Strain.Eigen.Values[i]);
		    
		    rStressDerivatives(i,j) += 0.5 * rBulkModulus / (rVariables.Strain.Eigen.Values[j] * rVariables.Strain.Eigen.Values[j] * rVariables.Strain.Eigen.Values[i] * rVariables.Strain.Eigen.Values[i]);
		    rStressDerivatives(i,j) -= 2.0 * rBulkModulus * std::log(rVariables.Strain.Invariants.J) * this->msIdentityMatrix(i,j) / (rVariables.Strain.Eigen.Values[j] * rVariables.Strain.Eigen.Values[i] * rVariables.Strain.Eigen.Values[i] * rVariables.Strain.Eigen.Values[i]);
		}
	    }
	
	}
    }
    else if( rStressMeasure == ConstitutiveModelData::StressMeasure_Kirchhoff ){ //Strain.Matrix = LeftCauchyGreen (b)

	for(SizeType i=0; i<3; i++)
	{
	    for(SizeType j=0; j<3; j++)
	    {
		for(unsigned int p=0; p<size; p++)
		{
		    const double& mu_p = rModelParameters[p];
		    const double& alpha_p = rModelParameters[p+size];
		    double f = athird * ( std::pow(rVariables.Strain.Eigen.Values[0],alpha_p) + std::pow(rVariables.Strain.Eigen.Values[1],alpha_p) + std::pow(rVariables.Strain.Eigen.Values[2],alpha_p) );

		    rStressDerivatives(i,j) += (mu_p * alpha_p * std::pow(rVariables.Strain.Invariants.J,(-alpha_p*athird)) * ( f - std::pow(rVariables.Strain.Eigen.Values[i],alpha_p) - std::pow(rVariables.Strain.Eigen.Values[j],alpha_p) + 3.0 * std::pow(rVariables.Strain.Eigen.Values[i],alpha_p) * this->msIdentityMatrix(i,j) ) / 6.0 + 0.5 * rBulkModulus) / (rVariables.Strain.Eigen.Values[j] * rVariables.Strain.Eigen.Values[j]);
	      
		}
	    }
	
	}
    }

    KRATOS_CATCH(" ")
  }
    
  //************************************************************************************
  //************************************************************************************
   
  void OgdenModel::CalculateDerivativeFactors(array_1d<double,6>& rDerivativeFactors, MatrixType& rStressDerivatives, array_1d<double,3>& rStressEigenValues, array_1d<double,3>& rStrainEigenValues, array_1d<unsigned int,3>& rPermutation)
  {
    KRATOS_TRY

    //use strain invariants as storage or the stress eigenvalues
    double& t1 = rStressEigenValues[rPermutation[0]];
    //double& t2 = rStressEigenValues[rPermutation[1]];
    double& t3 = rStressEigenValues[rPermutation[2]];
    //get strain eigen values
    double& b1 = rStrainEigenValues[rPermutation[0]]; 
    //double& b2 = rStrainEigenValues[rPermutation[1]]; 
    double& b3 = rStrainEigenValues[rPermutation[2]]; 
	
    rDerivativeFactors[0] = (t1-t3)/((b1-b3)*(b1-b3)) + (rStressDerivatives(rPermutation[2],rPermutation[1]) - rStressDerivatives(rPermutation[2],rPermutation[2]))/(b1-b3);
    rDerivativeFactors[1] =  2*b3*(t1-t3)/((b1-b3)*(b1-b3)) + (rStressDerivatives(rPermutation[2],rPermutation[1]) - rStressDerivatives(rPermutation[2],rPermutation[2]))*(b1+b3)/(b1-b3);
    rDerivativeFactors[2] =  2*(t1-t3)/((b1-b3)*(b1-b3)*(b1-b3)) + (rStressDerivatives(rPermutation[0],rPermutation[2]) + rStressDerivatives(rPermutation[2],rPermutation[0]) - rStressDerivatives(rPermutation[0],rPermutation[0]) - rStressDerivatives(rPermutation[2],rPermutation[2]))/((b1-b3)*(b1-b3)); 

	
    rDerivativeFactors[3] = rDerivativeFactors[2] * b3 + (rStressDerivatives(rPermutation[0],rPermutation[2]) - rStressDerivatives(rPermutation[2],rPermutation[1]))/(b1-b3);    
    rDerivativeFactors[4] = rDerivativeFactors[2] * b3 + (rStressDerivatives(rPermutation[2],rPermutation[0]) - rStressDerivatives(rPermutation[2],rPermutation[1]))/(b1-b3); 
    rDerivativeFactors[5] = 2*b3*b3*(t1-t3)/((b1-b3)*(b1-b3)*(b1-b3)) + (rStressDerivatives(rPermutation[0],rPermutation[2]) - rStressDerivatives(rPermutation[2],rPermutation[0]))*b1*b3/((b1-b3)*(b1-b3)) - (rStressDerivatives(rPermutation[0],rPermutation[0]) + rStressDerivatives(rPermutation[2],rPermutation[2]))*b3*b3/((b1-b3)*(b1-b3)) - rStressDerivatives(rPermutation[2],rPermutation[1])*(b1+b3)/(b1-b3);

    // std::cout<<" Alpha["<<rVariables.Factors.Alpha1<<" "<<rVariables.Factors.Alpha2<<" "<<rVariables.Factors.Alpha3<<"]"<<std::endl;
    // std::cout<<" Beta ["<<rVariables.Factors.Beta1<<" "<<rVariables.Factors.Beta2<<" "<<rVariables.Factors.Beta3<<"]"<<std::endl;
    
    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************
  
  void OgdenModel::CalculateAndAddConstitutiveTensor(HyperElasticDataType& rVariables, Matrix& rConstitutiveMatrix)
  {
    KRATOS_TRY
              
    //Calculate Ogden ConstitutiveMatrix
    const ModelDataType&  rModelData        = rVariables.GetModelData();
    const SizeType&       rVoigtSize        = rModelData.GetVoigtSize();      
    const VoigtIndexType& rIndexVoigtTensor = rModelData.GetVoigtIndexTensor();
    
    //Calculate constitutive components OPTION 2
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
    
   double& OgdenModel::AddConstitutiveComponent(HyperElasticDataType& rVariables, double &rCabcd,
						const unsigned int& a, const unsigned int& b,
						const unsigned int& c, const unsigned int& d)
  {
    KRATOS_TRY

    const ModelDataType&  rModelData  = rVariables.GetModelData();
    const MatrixType& rStressMatrix   = rModelData.GetStressMatrix(); //stress stored as StressMatrix
    
    //Calculate Stress main streches derivatives
    MatrixType StressDerivatives;
    noalias(StressDerivatives)=ZeroMatrix(3,3);
    this->CalculateMainStressDerivatives(rVariables, StressDerivatives);

    double Cabcd = 0;
    double Cabce = 0;
   
    for(unsigned int e=0; e<3; e++)
    {
	Cabce = CalculateIsotropicTensorDerivative(rVariables,StressDerivatives,Cabce,a,b,c,e);	
	//std::cout<<" Cabce e["<<e<<"]: "<<Cabce<<std::endl;
	Cabcd += Cabce * rVariables.Strain.Matrix(e,d);
    }
	
    Cabcd *= 2.0;    	    
    Cabcd -= (1.0/rVariables.Strain.Invariants.J)*rStressMatrix(a,d)*this->msIdentityMatrix(b,c);
    //std::cout<<" Cabcd "<<Cabcd<<" stress "<<rStressMatrix(a,d)<<std::endl;
    
    rCabcd += Cabcd;
    //std::cout<<" rCabcd "<<rCabcd<<" "<<a<<" "<<b<<" "<<c<<" "<<d<<std::endl;
    return rCabcd;
	
    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************
    
  double& OgdenModel::CalculateIsotropicTensorDerivative(HyperElasticDataType& rVariables,
							 MatrixType& rStressDerivatives, double &rCabcd,
							 const unsigned int& a, const unsigned int& b,
							 const unsigned int& c, const unsigned int& d)
  {
    KRATOS_TRY

    //strain eigenvalues	
    array_1d<double,3> StrainEigenValues;
    for (unsigned int i = 0; i < 3; i++)
	StrainEigenValues[i] = rVariables.Strain.Eigen.Values[i] * rVariables.Strain.Eigen.Values[i];
    
    array_1d<double,3> StressEigenValues;
    noalias(StressEigenValues)=ZeroVector(3);
    this->CalculateMainStresses(rVariables, StressEigenValues);
        
    if( ConstitutiveModelUtilities::AreEqual(StrainEigenValues[0],StrainEigenValues[1]) && ConstitutiveModelUtilities::AreEqual(StrainEigenValues[1],StrainEigenValues[2]) ){

	//std::cout<<" option all equal "<<StrainEigenValues<<std::endl;
	double Cabcd = 0;
	
	Cabcd = ConstitutiveModelUtilities::CalculateFourthOrderUnitTensor(this->msIdentityMatrix,Cabcd,a,b,c,d);

	rCabcd += (rStressDerivatives(0,0)-rStressDerivatives(0,1)) * Cabcd;
	    
	Cabcd = ConstitutiveModelUtilities::CalculateFourthOrderTensorProduct(this->msIdentityMatrix,this->msIdentityMatrix,Cabcd,a,b,c,d);
	rCabcd += rStressDerivatives(0,1) * Cabcd;

    }
    else{
	array_1d<unsigned int,3> Permuta;
	bool set = false;
    	if( ConstitutiveModelUtilities::AreEqual(StrainEigenValues[0],StrainEigenValues[1]) ){	    
    	    Permuta[0] = 2;
    	    Permuta[1] = 0;
    	    Permuta[2] = 1;
    	}
    	else if( ConstitutiveModelUtilities::AreEqual(StrainEigenValues[2],StrainEigenValues[0]) ){
    	    Permuta[0] = 1;
    	    Permuta[1] = 2;
    	    Permuta[2] = 0;
	}
	else if( ConstitutiveModelUtilities::AreEqual(StrainEigenValues[1],StrainEigenValues[2]) ){
    	    Permuta[0] = 0;
    	    Permuta[1] = 1;
            Permuta[2] = 2;
    	}
    	else{

	    array_1d<double,3> EigenVectorA;
	    array_1d<double,3> EigenVectorB;    
	    
	    array_1d<unsigned int,5> Permutation;
	    Permutation[0] = 0;
	    Permutation[1] = 1;
	    Permutation[2] = 2;
	    Permutation[3] = 0;
	    Permutation[4] = 1;
	    
	    double Dabcd = 0;
	    double Cabcd = 0;
	    for(unsigned int i=0; i<3; i++)
	    {
		double alpha = StressEigenValues[i]/((StrainEigenValues[i]-StrainEigenValues[Permutation[i+1]])*(StrainEigenValues[i]-StrainEigenValues[Permutation[i+2]]));
		
		Dabcd = ConstitutiveModelUtilities::CalculateSquareTensorDerivative(rVariables.Strain.Matrix,this->msIdentityMatrix,Dabcd,a,b,c,d);
		Cabcd  = Dabcd;
		Dabcd = ConstitutiveModelUtilities::CalculateFourthOrderUnitTensor(this->msIdentityMatrix,Dabcd,a,b,c,d);
		Cabcd -= (StrainEigenValues[Permutation[i+1]]+StrainEigenValues[Permutation[i+2]]) * Dabcd;
		noalias(EigenVectorA) = matrix_row<const MatrixType>(rVariables.Strain.Eigen.Vectors,Permutation[i]);		
		Dabcd = ConstitutiveModelUtilities::CalculateFourthOrderTensorProduct(EigenVectorA,EigenVectorA,Dabcd,a,b,c,d);
		Cabcd -= ((StrainEigenValues[i]-StrainEigenValues[Permutation[i+1]])+(StrainEigenValues[i]-StrainEigenValues[Permutation[i+2]])) * Dabcd;
		noalias(EigenVectorA) = matrix_row<const MatrixType>(rVariables.Strain.Eigen.Vectors,Permutation[i+1]);
		Dabcd = ConstitutiveModelUtilities::CalculateFourthOrderTensorProduct(EigenVectorA,EigenVectorA,Dabcd,a,b,c,d);
		Cabcd -= (StrainEigenValues[Permutation[i+1]]-StrainEigenValues[Permutation[i+2]]) * Dabcd;
		noalias(EigenVectorA) = matrix_row<const MatrixType>(rVariables.Strain.Eigen.Vectors,Permutation[i+2]);
		Dabcd = ConstitutiveModelUtilities::CalculateFourthOrderTensorProduct(EigenVectorA,EigenVectorA,Dabcd,a,b,c,d);
		Cabcd += (StrainEigenValues[Permutation[i+1]]-StrainEigenValues[Permutation[i+2]]) * Dabcd;
		rCabcd += Cabcd * alpha;
		
	    }

	    for(unsigned int i=0; i<3; i++)
	    {
		noalias(EigenVectorA) = matrix_row<const MatrixType>(rVariables.Strain.Eigen.Vectors,i);
		for(unsigned int j=0; j<3; j++)
		{
		    noalias(EigenVectorB) = matrix_row<const MatrixType>(rVariables.Strain.Eigen.Vectors,j);	      
		    Dabcd = ConstitutiveModelUtilities::CalculateFourthOrderTensorProduct(EigenVectorA,EigenVectorB,Dabcd,a,b,c,d);

		    rCabcd += rStressDerivatives(i,j) * Dabcd;
		}
	    }
	    set = true;

	    //std::cout<<" option all different "<<StrainEigenValues<<std::endl;
    	}

	double Cabcd = 0;
	if( !set ){
	    
	    array_1d<double,6> DerivativeFactors;
	    this->CalculateDerivativeFactors(DerivativeFactors,rStressDerivatives,StrainEigenValues,StressEigenValues,Permuta);

	    //std::cout<<" option some different "<<StrainEigenValues<<std::endl;
	    
	    Cabcd = ConstitutiveModelUtilities::CalculateSquareTensorDerivative(rVariables.Strain.Matrix,this->msIdentityMatrix,Cabcd,a,b,c,d);
	    rCabcd += DerivativeFactors[0] * Cabcd;
	    Cabcd = ConstitutiveModelUtilities::CalculateFourthOrderUnitTensor(this->msIdentityMatrix,Cabcd,a,b,c,d);
	    rCabcd -= DerivativeFactors[1] * Cabcd;
	    Cabcd = ConstitutiveModelUtilities::CalculateFourthOrderTensorProduct(rVariables.Strain.Matrix,rVariables.Strain.Matrix,Cabcd,a,b,c,d);
	    rCabcd -= DerivativeFactors[2] * Cabcd;
	    Cabcd = ConstitutiveModelUtilities::CalculateFourthOrderTensorProduct(rVariables.Strain.Matrix,this->msIdentityMatrix,Cabcd,a,b,c,d);
	    rCabcd += DerivativeFactors[3] * Cabcd;
	    Cabcd = ConstitutiveModelUtilities::CalculateFourthOrderTensorProduct(this->msIdentityMatrix,rVariables.Strain.Matrix,Cabcd,a,b,c,d);
	    rCabcd += DerivativeFactors[4] * Cabcd;
	    Cabcd = ConstitutiveModelUtilities::CalculateFourthOrderTensorProduct(this->msIdentityMatrix,this->msIdentityMatrix,Cabcd,a,b,c,d);
	    rCabcd -= DerivativeFactors[5] * Cabcd;
	}
    }

    return rCabcd;
    
    KRATOS_CATCH(" ")
  }

    
  //************************************************************************************
  //************************************************************************************

  int OgdenModel::Check(const Properties& rMaterialProperties, const ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    return 0;

    KRATOS_CATCH(" ")
  }



} // Namespace Kratos
