//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_models/elasticity_models/borja_model.hpp"


namespace Kratos
{

   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   BorjaModel::BorjaModel()
      : HenckyHyperElasticModel()
   {
   }

   //******************************COPY CONSTRUCTOR**************************************
   //************************************************************************************

   BorjaModel::BorjaModel(const BorjaModel& rOther)
      : HenckyHyperElasticModel( rOther )
   {
   }

   //********************************CLONE***********************************************
   //************************************************************************************

   ConstitutiveModel::Pointer BorjaModel::Clone() const
   {
      return Kratos::make_shared<BorjaModel>(*this);
   }

   //********************************ASSIGNMENT******************************************
   //************************************************************************************
   BorjaModel& BorjaModel::operator=(BorjaModel const& rOther)
   {
      HenckyHyperElasticModel::operator=(rOther);
      return *this;
   }

   //*******************************DESTRUCTOR*******************************************
   //************************************************************************************

   BorjaModel::~BorjaModel()
   {
   }


   // ****************************************************************
   // ***************************************************************
   // not in the correct place
   void BorjaModel::SeparateVolumetricAndDeviatoricPart( const MatrixType& rA, double & rVolumetric, MatrixType& rDev, double & devNorm)
   {
      KRATOS_TRY

      rVolumetric = 0;
      for (unsigned int i = 0; i < 3; i++)
         rVolumetric += rA(i,i);

      noalias(rDev) = rA;
      for (unsigned int i = 0; i < 3; i++)
         rDev(i,i) -= rVolumetric/3.0;

      devNorm = 0;
      for (unsigned int i = 0; i < 3; i++)
         for (unsigned int j = 0; j < 3; j++)
            devNorm += pow( rDev(i,j), 2);

      devNorm = sqrt(devNorm);


      KRATOS_CATCH("")
   }



   //************************************************************************************
   //************************************************************************************
   // CalculateAndAddStressTensor
   void BorjaModel::CalculateAndAddStressTensor( HyperElasticDataType& rVariables, MatrixType& rStressMatrix)
   {
      KRATOS_TRY

      // material parameters

      const ModelDataType & rModelData = rVariables.GetModelData();
      const Properties    & rMaterialProperties = rModelData.GetMaterialProperties();

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
      double ShearModulus = rAlphaShear * ReferencePressure * std::exp( -VolumetricHencky / rSwellingSlope );
      rStressMatrix += DeviatoricHencky * 2 * ( ShearModulus  + rConstantShearModulus );

      // 3.b Compute Volumetric Part
      double pressure = -ReferencePressure * std::exp( -VolumetricHencky / rSwellingSlope ) * ( 1 + rAlphaShear * pow(deviatoricNorm, 2) / rSwellingSlope);

      for (unsigned int i = 0; i< 3; i++)
         rStressMatrix(i,i) += pressure;

      KRATOS_CATCH("")
   }


   //************************************************************************************
   //************************************************************************************
   // CalculateAndAddConstitutiveTensor
   void BorjaModel::CalculateAndAddConstitutiveTensor( HyperElasticDataType& rVariables, Matrix& rConstitutiveMatrix)
   {
      KRATOS_TRY

      // material parameters

      const ModelDataType & rModelData = rVariables.GetModelData();
      const Properties    & rMaterialProperties = rModelData.GetMaterialProperties();

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


      // bulk modulus part
      double pressure = -ReferencePressure * std::exp( -VolumetricHencky / rSwellingSlope ) * ( 1 + rAlphaShear * pow(deviatoricNorm, 2) / rSwellingSlope);
      rConstitutiveMatrix = (-pressure/rSwellingSlope) * IdentityCross;

      // Shear modulus part
      rConstitutiveMatrix += 2.0*( rAlphaShear*ReferencePressure*std::exp(-VolumetricHencky/rSwellingSlope) + rConstantShearModulus) *(FourthOrderIdentity - (1.0/3.0)*IdentityCross);

      // coupling part
      Vector StrainVector = ZeroVector(6);
      StrainVector = ConstitutiveModelUtilities::StressTensorToVector( DeviatoricHencky, StrainVector); // then I do not have to divide by 2

      double Modulus = 2.0 * ReferencePressure * exp( - VolumetricHencky/rSwellingSlope) * ( rAlphaShear / rSwellingSlope);

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
