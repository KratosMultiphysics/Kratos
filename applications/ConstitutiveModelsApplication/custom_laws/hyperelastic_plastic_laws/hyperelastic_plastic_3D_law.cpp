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
#include "custom_laws/hyperelastic_plastic_laws/hyperelastic_plastic_3D_law.hpp"


namespace Kratos
{

  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

  HyperElasticPlastic3DLaw::HyperElasticPlastic3DLaw()
    : HyperElastic3DLaw()
  {
    KRATOS_TRY

    KRATOS_CATCH(" ")	  
  }
  
  //******************************CONSTRUCTOR WITH THE MODEL****************************
  //************************************************************************************

  HyperElasticPlastic3DLaw::HyperElasticPlastic3DLaw(ModelTypePointer pModel)
    : HyperElastic3DLaw(pModel)
  {
    KRATOS_TRY
      
    KRATOS_CATCH(" ")
  }

  //******************************COPY CONSTRUCTOR**************************************
  //************************************************************************************

  HyperElasticPlastic3DLaw::HyperElasticPlastic3DLaw(const HyperElasticPlastic3DLaw& rOther)
    : HyperElastic3DLaw(rOther)   
  {    
  }

  //********************************CLONE***********************************************
  //************************************************************************************

  ConstitutiveLaw::Pointer HyperElasticPlastic3DLaw::Clone() const
  {
    return ( ConstitutiveLaw::Pointer(new HyperElasticPlastic3DLaw(*this)) );
  }

  //*******************************DESTRUCTOR*******************************************
  //************************************************************************************

  HyperElasticPlastic3DLaw::~HyperElasticPlastic3DLaw()
  {
  }
  
 

  //************* COMPUTING  METHODS
  //************************************************************************************
  //************************************************************************************

  void  HyperElasticPlastic3DLaw::InitializeModelData(Parameters& rValues,ModelDataType& rModelValues)
  {
    KRATOS_TRY
      
    rModelValues.SetOptions(rValues.GetOptions());
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

  void  HyperElasticPlastic3DLaw::FinalizeModelData(Parameters& rValues,ModelDataType& rModelValues)
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
      mCauchyGreenVector[0] = rModelValues.StressMatrix(0,0);
      mCauchyGreenVector[1] = rModelValues.StressMatrix(1,1);
      mCauchyGreenVector[2] = rModelValues.StressMatrix(2,2);
      
      mCauchyGreenVector[3] = rModelValues.StressMatrix(0,1);
      mCauchyGreenVector[4] = rModelValues.StressMatrix(1,2);
      mCauchyGreenVector[5] = rModelValues.StressMatrix(2,0);
      
      mCauchyGreenVector   *=  ( 1.0 / rModelValues.MaterialParameters.LameMu );
      
      double VolumetricPart = (rModelValues.StrainMatrix(0,0)+rModelValues.StrainMatrix(1,1)+rModelValues.StrainMatrix(2,2))/3.0;
      
      mCauchyGreenVector[0] += VolumetricPart;
      mCauchyGreenVector[1] += VolumetricPart;
      mCauchyGreenVector[2] += VolumetricPart;
    }
    
    KRATOS_CATCH(" ")      
  }


  //*****************************MATERIAL RESPONSES*************************************
  //************************************************************************************


  void  HyperElasticPlastic3DLaw::CalculateMaterialResponsePK2(Parameters& rValues)
  {
    KRATOS_TRY

    Elastic3DLaw::CalculateMaterialResponsePK2(rValues);
      
    KRATOS_CATCH(" ")      
  }

  
  //*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
  //************************************************************************************

  void HyperElasticPlastic3DLaw::GetLawFeatures(Features& rFeatures)
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



} // Namespace Kratos
