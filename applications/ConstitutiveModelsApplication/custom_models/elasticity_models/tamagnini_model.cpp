//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                MCiantia $
//   Date:                $Date:                    JULY 2018 $
//   Revision:            $Revision:                      0.0 $
//
//



// 1. it has not been tested for alpha > 0!!!!

// System includes

// External includes

// Project includes
#include "custom_models/elasticity_models/tamagnini_model.hpp"


namespace Kratos
{

   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   TamagniniModel::TamagniniModel()
      : BorjaModel()
   {
   }

   //******************************COPY CONSTRUCTOR**************************************
   //************************************************************************************

   TamagniniModel::TamagniniModel(const TamagniniModel& rOther)
      : BorjaModel( rOther )
   {
   }

   //********************************CLONE***********************************************
   //************************************************************************************

   ConstitutiveModel::Pointer TamagniniModel::Clone() const
   {
      return ( TamagniniModel::Pointer(new TamagniniModel(*this)) );
   }

   //********************************ASSIGNMENT******************************************
   //************************************************************************************
   TamagniniModel& TamagniniModel::operator=(TamagniniModel const& rOther)
   {
      BorjaModel::operator=(rOther);
      return *this;
   }

   //*******************************DESTRUCTOR*******************************************
   //************************************************************************************

   TamagniniModel::~TamagniniModel()
   {
   }




   //************************************************************************************
   //************************************************************************************
   // CalculateAndAddStressTensor
   void TamagniniModel::CalculateAndAddStressTensor( HyperElasticDataType& rVariables, MatrixType& rStressMatrix)
   {
      KRATOS_TRY

      // material parameters

      const ModelDataType & rModelData = rVariables.GetModelData();
      const Properties    & rMaterialProperties = rModelData.GetProperties();

      const double & rSwellingSlope = rMaterialProperties[SWELLING_SLOPE];
      const double & rAlphaShear = rMaterialProperties[ALPHA_SHEAR];
      const double & rReferencePressure = rMaterialProperties[PRE_CONSOLIDATION_STRESS];
      const double & rOCR = rMaterialProperties[OVER_CONSOLIDATION_RATIO];
      const double ReferencePressure = rReferencePressure / rOCR;
      const double & rConstantShearModulus = rMaterialProperties[INITIAL_SHEAR_MODULUS];

      const MatrixType& HenckyStrain = rVariables.Strain.Matrix;

      // 2.a Separate Volumetric and deviatoric part
      double VolumetricHencky;
      MatrixType DeviatoricHencky(3,3);
      noalias(DeviatoricHencky) = ZeroMatrix(3,3);
      double deviatoricNorm;

      SeparateVolumetricAndDeviatoricPart( HenckyStrain, VolumetricHencky, DeviatoricHencky, deviatoricNorm);

      // 3.a Compute Deviatoric Part
      rStressMatrix.clear();
      double Phi = 1;


      if ( -VolumetricHencky > rSwellingSlope) {
         Phi = -rSwellingSlope * ReferencePressure * exp( -VolumetricHencky / rSwellingSlope - 1.0);
      }   else {
         Phi = -ReferencePressure * VolumetricHencky - ReferencePressure * pow( VolumetricHencky - rSwellingSlope, 2) / 2.0 / rSwellingSlope;
      }
      double ShearModulus = rAlphaShear * ReferencePressure * std::exp( -VolumetricHencky / rSwellingSlope );
      rStressMatrix += 2.0 * (  rConstantShearModulus + rAlphaShear / rSwellingSlope * Phi) * DeviatoricHencky;


      // 3.b Compute Volumetric Part
      double Theta = 1; 
      if ( -VolumetricHencky > rSwellingSlope) {
         Theta = -ReferencePressure * exp( -VolumetricHencky / rSwellingSlope - 1.0);
      }   else {
         Theta = -ReferencePressure * ( -VolumetricHencky / rSwellingSlope);
      }


      double pressure = ( 1 + rAlphaShear * pow(deviatoricNorm, 2) / rSwellingSlope) * Theta;

      for (unsigned int i = 0; i< 3; i++)
         rStressMatrix(i,i) += pressure;


      KRATOS_CATCH("")
   }


   //************************************************************************************
   //************************************************************************************
   // CalculateAndAddConstitutiveTensor
   void TamagniniModel::CalculateAndAddConstitutiveTensor( HyperElasticDataType& rVariables, Matrix& rConstitutiveMatrix)
   {
      KRATOS_TRY

      // material parameters

      const ModelDataType & rModelData = rVariables.GetModelData();
      const Properties    & rMaterialProperties = rModelData.GetProperties();

      const double & rSwellingSlope = rMaterialProperties[SWELLING_SLOPE];
      const double & rAlphaShear = rMaterialProperties[ALPHA_SHEAR];
      const double & rReferencePressure = rMaterialProperties[PRE_CONSOLIDATION_STRESS];
      const double & rOCR = rMaterialProperties[OVER_CONSOLIDATION_RATIO];
      const double ReferencePressure = rReferencePressure / rOCR;
      const double & rConstantShearModulus = rMaterialProperties[INITIAL_SHEAR_MODULUS];

      // 1. Define some matrices
      Matrix FourthOrderIdentity = ZeroMatrix(6,6);
      for (unsigned int i = 0; i<3; ++i)
         FourthOrderIdentity(i,i) = 1.0;

      for (unsigned int i = 3; i<6; ++i)
         FourthOrderIdentity(i,i) = 0.50;

      Matrix IdentityCross = ZeroMatrix(6,6);
      for (unsigned int i = 0; i<3; ++i) {
         for (unsigned int j = 0; j<3; ++j) {
            IdentityCross(i,j) = 1.0;
         }
      }

      // 2. Get Elastic Left Cauchy Green tensor
      const MatrixType& HenckyStrain = rVariables.Strain.Matrix;

      // 2.a Separate Volumetric and deviatoric part
      double VolumetricHencky;
      MatrixType DeviatoricHencky(3,3);
      noalias(DeviatoricHencky) = ZeroMatrix(3,3);
      double deviatoricNorm;

      SeparateVolumetricAndDeviatoricPart( HenckyStrain, VolumetricHencky, DeviatoricHencky, deviatoricNorm);

      double Phi = 1;
      if ( -VolumetricHencky > rSwellingSlope) {
         Phi = -rSwellingSlope * ReferencePressure * exp( -VolumetricHencky / rSwellingSlope - 1.0);
      }   else {
         Phi = -ReferencePressure * VolumetricHencky - ReferencePressure * pow( VolumetricHencky - rSwellingSlope, 2) / 2.0 / rSwellingSlope;
      }
      double Theta = 1; 
      if ( -VolumetricHencky > rSwellingSlope) {
         Theta = -ReferencePressure * exp( -VolumetricHencky / rSwellingSlope - 1.0);
      }   else {
         Theta = -ReferencePressure * ( -VolumetricHencky / rSwellingSlope);
      }
      double K = 1; 
      if ( -VolumetricHencky > rSwellingSlope) {
         K = ReferencePressure / rSwellingSlope  * exp( -VolumetricHencky / rSwellingSlope - 1.0);
      }   else {
         K = ReferencePressure / rSwellingSlope;
      }

      // bulk modulus part
      double pressure = -ReferencePressure * std::exp( -VolumetricHencky / rSwellingSlope ) * ( 1 + rAlphaShear * pow(deviatoricNorm, 2) / rSwellingSlope);
      rConstitutiveMatrix = ( 1 + rAlphaShear / rSwellingSlope * pow(deviatoricNorm,2.0) ) * K * IdentityCross;

      // Shear modulus part
      rConstitutiveMatrix += 2.0*(  rConstantShearModulus + rAlphaShear / rSwellingSlope * Phi) *(FourthOrderIdentity - (1.0/3.0)*IdentityCross);

      // coupling part
      Vector StrainVector = ZeroVector(6);
      StrainVector = ConstitutiveModelUtilities::StressTensorToVector( DeviatoricHencky, StrainVector); // then I do not have to divide by 2

      double Modulus = 2.0 * rAlphaShear / rSwellingSlope * Theta; 

      for (unsigned int i = 0; i<3; ++i) {
         for (unsigned int j = 0; j<3; ++j) {
            rConstitutiveMatrix(i,j) -= Modulus * (StrainVector(i) ); 
            rConstitutiveMatrix(i,j) -= Modulus * (StrainVector(j) ); 
         }
      }

      for (unsigned int i = 0; i<3; ++i) {
         for (unsigned int j = 3; j < 6; ++j) {
            rConstitutiveMatrix(i,j) -= Modulus*(StrainVector(j));
         }
      }

      for (unsigned int i = 3; i<6; ++i) {
         for (unsigned int j = 0; j<3; ++j) {
            rConstitutiveMatrix(i,j) -= Modulus*(StrainVector(i));
         }
      }


      rVariables.State().Set(ConstitutiveModelData::CONSTITUTIVE_MATRIX_COMPUTED);

      KRATOS_CATCH("")
   }


} // Namespace Kratos
