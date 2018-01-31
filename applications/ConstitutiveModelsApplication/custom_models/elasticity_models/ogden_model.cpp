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

  void OgdenModel::CalculateAndAddStressTensor(HyperElasticDataType& rVariables, MatrixType& rStressMatrix)
  {
    KRATOS_TRY

    const ModelDataType&  rModelData        = rVariables.GetModelData();
    const StressMeasureType& rStressMeasure = rModelData.GetStressMeasure();

    if( rStressMeasure == ConstitutiveModelData::StressMeasure_PK2 ){ //Strain.Matrix = RightCauchyGreen (C)

	const MaterialDataType& rMaterial = rVariables.GetMaterialParameters();

	const double& rBulkModulus = rMaterial.GetBulkModulus();
	const std::vector<double>& rModelParameters = rMaterial.GetModelParameters(); //nu values, lambda values
	
	array_1d<double,3> MainStresses;
	array_1d<double,3> EigenVector;
	unsigned int size = (rModelParameters.size()/2.0);
	double athird = 1.0/3.0;
	for(unsigned int i=0; i<3; i++)
	{
	    MainStresses[i]=0.0;
	    for(unsigned int p=0; p<size; p++)
	    {
		const double& mu_p = rModelParameters[p];
		const double& alpha_p = rModelParameters[p+size];
		MainStresses[i] += (mu_p * std::pow(rVariables.Strain.Invariants.J,(-alpha_p*athird)) * ( std::pow(rVariables.Strain.Eigen.Values[i],alpha_p) - athird * ( std::pow(rVariables.Strain.Eigen.Values[0],alpha_p) + std::pow(rVariables.Strain.Eigen.Values[1],alpha_p) + std::pow(rVariables.Strain.Eigen.Values[2],alpha_p) ) ) + rBulkModulus * std::log(rVariables.Strain.Invariants.J) ) / (rVariables.Strain.Eigen.Values[i]*rVariables.Strain.Eigen.Values[i]);
	    }

	    noalias(EigenVector) = matrix_row<const MatrixType>(rVariables.Strain.Eigen.Vectors,i);
	    rStressMatrix += MainStresses[i] * outer_prod(EigenVector,EigenVector);
	}

    }
    else if( rStressMeasure == ConstitutiveModelData::StressMeasure_Kirchhoff ){ //Strain.Matrix = LeftCauchyGreen (b)

	const MaterialDataType& rMaterial = rVariables.GetMaterialParameters();

	const double& rBulkModulus = rMaterial.GetBulkModulus();
	const std::vector<double>& rModelParameters = rMaterial.GetModelParameters(); //nu values, lambda values
      
	array_1d<double,3> MainStresses;
	array_1d<double,3> EigenVector;
	unsigned int size = (rModelParameters.size()/2.0);
	double athird = 1.0/3.0;
	for(unsigned int i=0; i<3; i++)
	{
	    MainStresses[i]=0.0;
	    for(unsigned int p=0; p<size; p++)
	    {
		const double& mu_p = rModelParameters[p];
		const double& alpha_p = rModelParameters[p+size];
		MainStresses[i] += mu_p * std::pow(rVariables.Strain.Invariants.J,(-alpha_p*athird)) * ( std::pow(rVariables.Strain.Eigen.Values[i],alpha_p) - athird * ( std::pow(rVariables.Strain.Eigen.Values[0],alpha_p) + std::pow(rVariables.Strain.Eigen.Values[1],alpha_p) + std::pow(rVariables.Strain.Eigen.Values[2],alpha_p) ) ) + rBulkModulus * std::log(rVariables.Strain.Invariants.J);
	    }

	    noalias(EigenVector) = matrix_row<const MatrixType>(rVariables.Strain.Eigen.Vectors, i);
	    rStressMatrix += MainStresses[i] * outer_prod(EigenVector,EigenVector);
	}

    }

    rVariables.State().Set(ConstitutiveModelData::STRESS_COMPUTED);

    KRATOS_CATCH(" ")
  }


  //***********************PROTECTED OPERATIONS FROM BASE CLASS*************************
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
    
    //Calculate Invariants
    //this->CalculateInvariants(rVariables);

    //Calculate LameMuBar
    //rValues.MaterialParameters.LameMuBar = rValues.MaterialParameters.LameMu * ( rVariables.Strain.Matrix(0,0) + rVariables.Strain.Matrix(1,1) + rVariables.Strain.Matrix(2,2) ) * rVariables.Strain.Invariants.J_13 * rVariables.Strain.Invariants.J_13 * (1.0/3.0) ;

    //std::cout<<" LameMuBar "<<rValues.MaterialParameters.LameMuBar<<std::endl;

    //Algorithmic moduli factors
    //this->CalculateScalingFactors(rVariables);


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

    //Calculate Stress main streches derivatives
    MatrixType StressDerivatives;
    noalias(StressDerivatives)=ZeroMatrix(3,3);

    this->CalculateMainStressDerivatives(rVariables, StressDerivatives);
    
    //Calculate Stress left cauchy green derivatives
    array_1d<double,3> StrainEigenValues;
    for (unsigned int i = 0; i < 3; i++)
	StrainEigenValues[i] = rVariables.Strain.Eigen.Values[i] * rVariables.Strain.Eigen.Values[i];
    
    MatrixType EigenValues;
    MatrixType EigenVectors;

    MathUtils<double>::EigenSystem<3>(rModelData.StressMatrix, EigenVectors, EigenValues);

    array_1d<double,3> StressEigenValues;
    for (unsigned int i = 0; i < 3; i++)
	StressEigenValues[i] = EigenValues(i,i);
    
    //Store stress eigenvalues in strain invariants variables in a vector
    array_1d<double,6> ScalingFactors;
    //get stress eigen values
    ScalingFactors[0] = StressEigenValues[0];
    ScalingFactors[1] = StressEigenValues[1];
    ScalingFactors[2] = StressEigenValues[2];
    //get strain eigen values
    ScalingFactors[3] = StrainEigenValues[0]; 
    ScalingFactors[4] = StrainEigenValues[1]; 
    ScalingFactors[5] = StrainEigenValues[2]; 

    
    if( ConstitutiveModelUtilities::IsEqual(StrainEigenValues[0],StrainEigenValues[1]) ){
	
	if( ConstitutiveModelUtilities::IsEqual(StrainEigenValues[1],StrainEigenValues[2]) ){
	    //Calculate constitutive components OPTION 3
	    for(SizeType i=0; i<rVoigtSize; i++)
	    {
		for(SizeType j=0; j<rVoigtSize; j++)
		{
		    rConstitutiveMatrix(i,j) = this->AddConstitutiveComponentC(rVariables,StressDerivatives,rConstitutiveMatrix(i,j),
									       rIndexVoigtTensor[i][0],rIndexVoigtTensor[i][1],
									       rIndexVoigtTensor[j][0],rIndexVoigtTensor[j][1]);
		}
		
	    }
	    rVariables.State().Set(ConstitutiveModelData::CONSTITUTIVE_MATRIX_COMPUTED,true);
	}
	else{	    
	    array_1d<unsigned int,3> Permutation;
	    Permutation[0] = 2;
	    Permutation[1] = 0;
	    Permutation[2] = 1;
	    this->CalculateScalingFactors(ScalingFactors,StressDerivatives,StrainEigenValues,StressEigenValues,Permutation);
	}
    }
    else{
	if( ConstitutiveModelUtilities::IsEqual(StrainEigenValues[0],StrainEigenValues[2]) ){
	    array_1d<unsigned int,3> Permutation;
	    Permutation[0] = 1;
	    Permutation[1] = 2;
	    Permutation[2] = 0;
	    this->CalculateScalingFactors(ScalingFactors,StressDerivatives,StrainEigenValues,StressEigenValues,Permutation);
	}
	else if( ConstitutiveModelUtilities::IsEqual(StrainEigenValues[1],StrainEigenValues[2]) ){
	    array_1d<unsigned int,3> Permutation;
	    Permutation[0] = 0;
	    Permutation[1] = 1;
            Permutation[2] = 2;
	    this->CalculateScalingFactors(ScalingFactors,StressDerivatives,StrainEigenValues,StressEigenValues,Permutation);
	}
	else{
	    //Calculate constitutive components OPTION 1
	    for(SizeType i=0; i<rVoigtSize; i++)
	    {
		for(SizeType j=0; j<rVoigtSize; j++)
		{
		    rConstitutiveMatrix(i,j) = this->AddConstitutiveComponentA(rVariables,ScalingFactors,StressDerivatives,rConstitutiveMatrix(i,j),
									       rIndexVoigtTensor[i][0],rIndexVoigtTensor[i][1],
									       rIndexVoigtTensor[j][0],rIndexVoigtTensor[j][1]);
		}
		
	    }
	    rVariables.State().Set(ConstitutiveModelData::CONSTITUTIVE_MATRIX_COMPUTED,true);
	}
    }

    
    if(rVariables.State().IsNot(ConstitutiveModelData::CONSTITUTIVE_MATRIX_COMPUTED)){

	//Calculate constitutive components OPTION 2
	for(SizeType i=0; i<rVoigtSize; i++)
	{
	    for(SizeType j=0; j<rVoigtSize; j++)
	    {
		rConstitutiveMatrix(i,j) = this->AddConstitutiveComponentB(rVariables,ScalingFactors,rConstitutiveMatrix(i,j),
									   rIndexVoigtTensor[i][0],rIndexVoigtTensor[i][1],
									   rIndexVoigtTensor[j][0],rIndexVoigtTensor[j][1]);
	    }
	
	}

	rVariables.State().Set(ConstitutiveModelData::CONSTITUTIVE_MATRIX_COMPUTED,true);
    }    
    
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
    MatrixType StressDerivatives;
    noalias(StressDerivatives)=ZeroMatrix(3,3);

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

		    StressDerivatives(i,j) += (mu_p * alpha_p * std::pow(rVariables.Strain.Invariants.J,(-alpha_p*athird)) * ( f - std::pow(rVariables.Strain.Eigen.Values[i],alpha_p) - std::pow(rVariables.Strain.Eigen.Values[j],alpha_p) + 3.0 * std::pow(rVariables.Strain.Eigen.Values[i],alpha_p) * this->msIdentityMatrix(i,j) ) / 6.0 + 0.5 * rBulkModulus ) / (rVariables.Strain.Eigen.Values[j] * rVariables.Strain.Eigen.Values[j] * rVariables.Strain.Eigen.Values[i] * rVariables.Strain.Eigen.Values[i]);
		
		    StressDerivatives(i,j) -= (mu_p * std::pow(rVariables.Strain.Invariants.J,(-alpha_p*athird)) * ( std::pow(rVariables.Strain.Eigen.Values[i],alpha_p) - f ) + rBulkModulus * std::log(rVariables.Strain.Invariants.J) ) * 2.0  / (rVariables.Strain.Eigen.Values[j] * rVariables.Strain.Eigen.Values[i] * rVariables.Strain.Eigen.Values[i] * rVariables.Strain.Eigen.Values[i]);
		    
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

		    StressDerivatives(i,j) += (mu_p * alpha_p * std::pow(rVariables.Strain.Invariants.J,(-alpha_p*athird)) * ( f - std::pow(rVariables.Strain.Eigen.Values[i],alpha_p) - std::pow(rVariables.Strain.Eigen.Values[j],alpha_p) + 3.0 * std::pow(rVariables.Strain.Eigen.Values[i],alpha_p) * this->msIdentityMatrix(i,j) ) / 6.0 + 0.5 * rBulkModulus) / (rVariables.Strain.Eigen.Values[j] * rVariables.Strain.Eigen.Values[j]);
	      
		}
	    }
	
	}
    }

    KRATOS_CATCH(" ")
  }
    
  //************************************************************************************
  //************************************************************************************
   
  void OgdenModel::CalculateScalingFactors(array_1d<double,6>& rScalingFactors, MatrixType& rStressDerivatives, array_1d<double,3>& rStressEigenValues, array_1d<double,3>& rStrainEigenValues, array_1d<unsigned int,3>& rPermutation)
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
	
    rScalingFactors[0] = (t1-t3)/((b1-b3)*(b1-b3)) + (rStressDerivatives(rPermutation[2],rPermutation[1]) - rStressDerivatives(rPermutation[2],rPermutation[2]))/(b1-b3);
    rScalingFactors[1] =  2*b3*(t1-t3)/((b1-b3)*(b1-b3)) + (rStressDerivatives(rPermutation[2],rPermutation[1]) - rStressDerivatives(rPermutation[2],rPermutation[2]))*(b1+b3)/(b1-b3);
    rScalingFactors[2] =  2*(t1-t3)/((b1-b3)*(b1-b3)*(b1-b3)) + (rStressDerivatives(rPermutation[0],rPermutation[2]) + rStressDerivatives(rPermutation[2],rPermutation[0]) - rStressDerivatives(rPermutation[0],rPermutation[0]) - rStressDerivatives(rPermutation[2],rPermutation[2]))/((b1-b3)*(b1-b3)); 

	
    rScalingFactors[3] = rScalingFactors[2] * b3 + (rStressDerivatives(rPermutation[0],rPermutation[2]) - rStressDerivatives(rPermutation[2],rPermutation[1]))/(b1-b3);    
    rScalingFactors[4] = rScalingFactors[2] * b3 + (rStressDerivatives(rPermutation[2],rPermutation[0]) - rStressDerivatives(rPermutation[2],rPermutation[1]))/(b1-b3); 
    rScalingFactors[5] = 2*b3*b3*(t1-t3)/((b1-b3)*(b1-b3)*(b1-b3)) + (rStressDerivatives(rPermutation[0],rPermutation[2]) - rStressDerivatives(rPermutation[2],rPermutation[0]))*b1*b3/((b1-b3)*(b1-b3)) - (rStressDerivatives(rPermutation[0],rPermutation[0]) + rStressDerivatives(rPermutation[2],rPermutation[2]))*b3*b3/((b1-b3)*(b1-b3)) - rStressDerivatives(rPermutation[2],rPermutation[1])*(b1+b3)/(b1-b3);

    // std::cout<<" Alpha["<<rVariables.Factors.Alpha1<<" "<<rVariables.Factors.Alpha2<<" "<<rVariables.Factors.Alpha3<<"]"<<std::endl;
    // std::cout<<" Beta ["<<rVariables.Factors.Beta1<<" "<<rVariables.Factors.Beta2<<" "<<rVariables.Factors.Beta3<<"]"<<std::endl;
    
    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************
    
    double& OgdenModel::AddConstitutiveComponentA(HyperElasticDataType& rVariables, array_1d<double,6>& rScalingFactors,
						  MatrixType& rStressDerivatives, double &rCabcd,
						  const unsigned int& a, const unsigned int& b,
						  const unsigned int& c, const unsigned int& d)
  {
    KRATOS_TRY

    //this is a simplified version, more terms of the derivatives are needed to be general:

    const ModelDataType&  rModelData  = rVariables.GetModelData();
    const MatrixType& rStressMatrix   = rModelData.GetStressMatrix(); //stress stored as StressMatrix

    double Cabcd = 0;
    double Dabce = 0;
    double Cabce = 0;
    double Eabce = 0;
    

    array_1d<double,3> EigenVectorA;
    array_1d<double,3> EigenVectorB;    

    array_1d<unsigned int,5> Permutation;
    Permutation[0] = 0;
    Permutation[1] = 1;
    Permutation[2] = 2;
    Permutation[3] = 0;
    Permutation[4] = 1;
    
    for(unsigned int e=0; e<3; e++)
    {
	Cabce = 0;	    
	for(unsigned int i=0; i<3; i++)
	{
	    double alpha = rScalingFactors[Permutation[i]]/((rScalingFactors[Permutation[i]+3]-rScalingFactors[Permutation[i+1]+3])*(rScalingFactors[Permutation[i]+3]-rScalingFactors[Permutation[i+2]+3]));
	    Dabce = ConstitutiveModelUtilities::CalculateSquareTensorDerivative(rVariables.Strain.Matrix,this->msIdentityMatrix,Dabce,a,b,c,e);
	    Eabce  = Dabce;
	    Dabce = ConstitutiveModelUtilities::CalculateFourthOrderUnitTensor(this->msIdentityMatrix,Dabce,a,b,c,e);
	    Eabce -= (rScalingFactors[Permutation[i+1]+3]+rScalingFactors[Permutation[i+2]+3]) * Dabce;
	    noalias(EigenVectorA) = matrix_row<const MatrixType>(rVariables.Strain.Eigen.Vectors,Permutation[i]);		
	    Dabce = ConstitutiveModelUtilities::CalculateFourthOrderTensorProduct(EigenVectorA,EigenVectorA,Dabce,a,b,c,e);
	    Eabce -= ((rScalingFactors[Permutation[i]+3]-rScalingFactors[Permutation[i+1]+3])+(rScalingFactors[Permutation[i]+3]-rScalingFactors[Permutation[i+2]+3])) * Dabce;
	    noalias(EigenVectorA) = matrix_row<const MatrixType>(rVariables.Strain.Eigen.Vectors,Permutation[i+1]);
	    Dabce = ConstitutiveModelUtilities::CalculateFourthOrderTensorProduct(EigenVectorA,EigenVectorA,Dabce,a,b,c,e);
	    Eabce -= (rScalingFactors[Permutation[i+1]+3]-rScalingFactors[Permutation[i+2]+3]) * Dabce;
	    noalias(EigenVectorA) = matrix_row<const MatrixType>(rVariables.Strain.Eigen.Vectors,Permutation[i+2]);
	    Dabce = ConstitutiveModelUtilities::CalculateFourthOrderTensorProduct(EigenVectorA,EigenVectorA,Dabce,a,b,c,e);
	    Eabce += (rScalingFactors[Permutation[i+1]+3]-rScalingFactors[Permutation[i+2]+3]) * Dabce;
	    Cabce += Eabce * alpha;
	}

	for(unsigned int i=0; i<3; i++)
	{     
	    for(unsigned int j=0; j<3; j++)
	    {
		noalias(EigenVectorA) = matrix_row<const MatrixType>(rVariables.Strain.Eigen.Vectors,i);
		noalias(EigenVectorB) = matrix_row<const MatrixType>(rVariables.Strain.Eigen.Vectors,j);	      
		Dabce = ConstitutiveModelUtilities::CalculateFourthOrderTensorProduct(EigenVectorA,EigenVectorB,Dabce,a,b,c,e);

		Cabce += rStressDerivatives(i,j) * Dabce;
	    }
	}

	//std::cout<<" Cabce e["<<e<<"]: "<<Cabce<<std::endl;
	Cabcd += Cabce * rVariables.Strain.Matrix(e,d);
    }
	
    Cabcd *= 2.0;    	    
    Cabcd -= rStressMatrix(a,d)*this->msIdentityMatrix(b,c);
    //std::cout<<" Cabcd "<<Cabcd<<std::endl;
    
    rCabcd += Cabcd;
    //std::cout<<" rCabcd "<<rCabcd<<" "<<a<<" "<<b<<" "<<c<<" "<<d<<std::endl;
    return rCabcd;


    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************
    
  double& OgdenModel::AddConstitutiveComponentB(HyperElasticDataType& rVariables,
						array_1d<double,6>& rScalingFactors, double &rCabcd,
						const unsigned int& a, const unsigned int& b,
						const unsigned int& c, const unsigned int& d)
  {
    KRATOS_TRY

    //this is a simplified version, more terms of the derivatives are needed to be general:

    const ModelDataType&  rModelData = rVariables.GetModelData();
    const MatrixType& rStressMatrix  = rModelData.GetStressMatrix(); //stress stored as StressMatrix

    double Cabcd = 0;
    double Dabce = 0;
    double Cabce = 0;
    
    for(unsigned int e=0; e<3; e++)
    {     	    	    
	Dabce = ConstitutiveModelUtilities::CalculateSquareTensorDerivative(rVariables.Strain.Matrix,this->msIdentityMatrix,Dabce,a,b,c,e);
	Cabce = rScalingFactors[0] * Dabce;
	
	Dabce = ConstitutiveModelUtilities::CalculateFourthOrderUnitTensor(this->msIdentityMatrix,Dabce,a,b,c,e);
	Cabce -= rScalingFactors[1] * Dabce;
	
	Dabce = ConstitutiveModelUtilities::CalculateFourthOrderTensorProduct(rVariables.Strain.Eigen.Vectors,rVariables.Strain.Eigen.Vectors,Dabce,a,b,c,e);
	Cabce -= rScalingFactors[2] * Dabce;
	
	Dabce = ConstitutiveModelUtilities::CalculateFourthOrderTensorProduct(rVariables.Strain.Eigen.Vectors,this->msIdentityMatrix,Dabce,a,b,c,e);
	Cabce += rScalingFactors[3] * Dabce;
	
	Dabce = ConstitutiveModelUtilities::CalculateFourthOrderTensorProduct(this->msIdentityMatrix,rVariables.Strain.Eigen.Vectors,Dabce,a,b,c,e);
	Cabce += rScalingFactors[4] * Dabce;
	
	Dabce = ConstitutiveModelUtilities::CalculateFourthOrderTensorProduct(this->msIdentityMatrix,this->msIdentityMatrix,Dabce,a,b,c,e);

	Cabce -= rScalingFactors[5] * Dabce;
	    
	Cabcd += Cabce * rVariables.Strain.Matrix(e,d);
    }

    Cabcd *= 2.0;    	    
    Cabcd -= rStressMatrix(a,d)*this->msIdentityMatrix(b,c);
    
    rCabcd += Cabcd;
    
    return rCabcd;


    KRATOS_CATCH(" ")
  }
    
  //************************************************************************************
  //************************************************************************************
    
  double& OgdenModel::AddConstitutiveComponentC(HyperElasticDataType& rVariables,
						MatrixType& rStressDerivatives, double &rCabcd,
						const unsigned int& a, const unsigned int& b,
						const unsigned int& c, const unsigned int& d)
  {
    KRATOS_TRY

    //this is a simplified version, more terms of the derivatives are needed to be general:

    const ModelDataType&  rModelData = rVariables.GetModelData();
    const MatrixType& rStressMatrix  = rModelData.GetStressMatrix(); //stress stored as StressMatrix
    
    double Cabce = 0;
    double Dabce = 0;
    double Cabcd = 0;
    

    for(unsigned int e=0; e<3; e++)
    {     	    
	Dabce = ConstitutiveModelUtilities::CalculateFourthOrderUnitTensor(this->msIdentityMatrix,Dabce,a,b,c,e);
	Cabce = (rStressDerivatives(0,0)-rStressDerivatives(0,1)) * Dabce;
	    
	Dabce = ConstitutiveModelUtilities::CalculateFourthOrderTensorProduct(this->msIdentityMatrix,this->msIdentityMatrix,Dabce,a,b,c,e);
	Cabce -= rStressDerivatives(0,1) * Dabce;

	Cabcd += Cabce * rVariables.Strain.Matrix(e,d);
    }

    Cabcd *= 2.0;    	    
    Cabcd -= rStressMatrix(a,d)*this->msIdentityMatrix(b,c);

    rCabcd += Cabcd;
    
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
