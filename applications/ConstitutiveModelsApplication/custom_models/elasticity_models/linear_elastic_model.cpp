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
    : ElasticityModel()
  {
  }

  //******************************COPY CONSTRUCTOR**************************************
  //************************************************************************************

  LinearElasticModel::LinearElasticModel(const LinearElasticModel& rOther)
    : ElasticityModel(rOther)
  {
  }

  //********************************CLONE***********************************************
  //************************************************************************************

  ElasticityModel::Pointer LinearElasticModel::Clone() const
  {
    return ( LinearElasticModel::Pointer(new LinearElasticModel(*this)) );
  }

  //********************************ASSIGNMENT******************************************
  //************************************************************************************
  LinearElasticModel& LinearElasticModel::operator=(LinearElasticModel const& rOther)
  {
    ElasticityModel::operator=(rOther);
    return *this;
  }
  
  //*******************************DESTRUCTOR*******************************************
  //************************************************************************************

  LinearElasticModel::~LinearElasticModel()
  {
  }

  //***********************PUBLIC OPERATIONS FROM BASE CLASS****************************
  //************************************************************************************

  void LinearElasticModel::InitializeElasticData(ModelDataType& rValues, ElasticDataType& rVariables)
  {
    KRATOS_TRY

    //set model data pointer
    rVariables.SetModelData(rValues);
    rVariables.SetState(rValues.State);

       
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
    StrainVector = ConstitutiveLawUtilities::StrainTensorToVector(rValues.StrainMatrix, StrainVector);

    this->CalculateAndAddConstitutiveTensor(Variables);
    
    VectorType StressVector;
    this->CalculateAndAddStressTensor(Variables,StrainVector,StressVector);

    rStressMatrix = ConstitutiveLawUtilities::VectorToSymmetricTensor(StressVector,rStressMatrix);
    
    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************

  void LinearElasticModel::CalculateAndAddStressTensor(ElasticDataType& rVariables, VectorType& rStrainVector, VectorType& rStressVector)
  {
    KRATOS_TRY

    noalias(rStressVector) = prod(rVariables.ConstitutiveTensor,rStrainVector);
      
    rVariables.State().Set(ConstitutiveModelData::COMPUTED_STRESS);
    
    KRATOS_CATCH(" ")
  }
  
  
  //************************************************************************************
  //************************************************************************************

  void LinearElasticModel::CalculateIsochoricStressTensor(ModelDataType& rValues, MatrixType& rStressMatrix)
  {
    KRATOS_TRY

    KRATOS_ERROR << "calling the base class function in LinearElasticModel ... illegal operation" << std::endl;    
    
    KRATOS_CATCH(" ")
  }
  
  //************************************************************************************
  //************************************************************************************

  void LinearElasticModel::CalculateAndAddIsochoricStressTensor(ElasticDataType& rVariables, VectorType& rStrainVector, VectorType& rStressVector)
  {
    KRATOS_TRY
     
    KRATOS_ERROR << "calling the base class function in LinearElasticModel ... illegal operation" << std::endl;
	
    KRATOS_CATCH(" ")
  }

  
  //************************************************************************************
  //************************************************************************************

  void LinearElasticModel::CalculateVolumetricStressTensor(ModelDataType& rValues, MatrixType& rStressMatrix)
  {
    KRATOS_TRY

    KRATOS_ERROR << "calling the base class function in LinearElasticModel ... illegal operation" << std::endl;
      
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
    this->CalculateAndAddConstitutiveTensor(Variables,rConstitutiveMatrix);
    	
    KRATOS_CATCH(" ")
  }

  
  //************************************************************************************
  //************************************************************************************
  
  void LinearElasticModel::CalculateAndAddConstitutiveTensor(ElasticDataType& rVariables, Matrix& rConstitutiveMatrix)
  {
    KRATOS_TRY
              
    this->CalculateAndAddConstitutiveTensor(rVariables);
    
    rConstitutiveMatrix = ConstitutiveLawUtilities::ConstitutiveTensorToMatrix(rVariables.ConstitutiveTensor,rConstitutiveMatrix);
      
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

    
    rVariables.State().Set(ConstitutiveModelData::COMPUTED_CONSTITUTIVE_MATRIX);

    
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
    StrainVector = ConstitutiveLawUtilities::StrainTensorToVector(rValues.StrainMatrix, StrainVector);

    this->CalculateAndAddConstitutiveTensor(Variables, rConstitutiveMatrix);
    
    VectorType StressVector;
    this->CalculateAndAddStressTensor(Variables,StrainVector,StressVector);

    rStressMatrix = ConstitutiveLawUtilities::VectorToSymmetricTensor(StressVector,rStressMatrix);

    
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
    
    rConstitutiveMatrix = ConstitutiveLawUtilities::ConstitutiveTensorToMatrix(rVariables.ConstitutiveTensor, rConstitutiveMatrix);
          
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
    
    rConstitutiveMatrix = ConstitutiveLawUtilities::ConstitutiveTensorToMatrix(rVariables.ConstitutiveTensor, rConstitutiveMatrix);
      
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
      
    if(YOUNG_MODULUS.Key() == 0 || rMaterialProperties[YOUNG_MODULUS]<= 0.00)
      KRATOS_ERROR << "YOUNG_MODULUS has Key zero or invalid value" << std::endl;

    const double& nu = rMaterialProperties[POISSON_RATIO];
    const bool check = bool( (nu >0.499 && nu<0.501 ) || (nu < -0.999 && nu > -1.01 ) );

    if(POISSON_RATIO.Key() == 0 || check==true)
      KRATOS_ERROR << "POISSON_RATIO has Key zero invalid value" << std::endl;
		
    return 0;

	  
    KRATOS_CATCH(" ")
  }
  


} // Namespace Kratos
