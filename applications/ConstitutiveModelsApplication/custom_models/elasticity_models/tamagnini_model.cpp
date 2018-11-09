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

      const double & rAlphaShear = rMaterialProperties[ALPHA_SHEAR];
      const double & rYoungModulus = rMaterialProperties[YOUNG_MODULUS];
      const double & rPoissonRatio = rMaterialProperties[POISSON_RATIO];
      const double & rReferencePressure = rMaterialProperties[REFERENCE_PRESSURE];

      const MatrixType& HenckyStrain = rVariables.Strain.Matrix;



      double BulkModulus = rYoungModulus / 3.0 / ( 1.0 - 2.0 * rPoissonRatio); 
      double SwellingSlope = rReferencePressure / BulkModulus;
      double ConstantShearModulus = rYoungModulus / 2.0 / ( 1.0 + rPoissonRatio);

      
      // 2.a Separate Volumetric and deviatoric part
      double VolumetricHencky;
      MatrixType DeviatoricHencky(3,3);
      noalias(DeviatoricHencky) = ZeroMatrix(3,3);
      double deviatoricNorm;

      SeparateVolumetricAndDeviatoricPart( HenckyStrain, VolumetricHencky, DeviatoricHencky, deviatoricNorm);

      // 3.a Compute Deviatoric Part
      rStressMatrix.clear();
      double Phi = 1;


      if ( -VolumetricHencky > SwellingSlope) {
         Phi = -SwellingSlope * rReferencePressure * exp( -VolumetricHencky / SwellingSlope - 1.0);
      }   else {
         Phi = -rReferencePressure * VolumetricHencky - rReferencePressure * pow( VolumetricHencky - SwellingSlope, 2) / 2.0 / SwellingSlope;
      }
      rStressMatrix += 2.0 * (  ConstantShearModulus + rAlphaShear / SwellingSlope * Phi) * DeviatoricHencky;


      // 3.b Compute Volumetric Part
      double Theta = 1; 
      if ( -VolumetricHencky > SwellingSlope) {
         Theta = -rReferencePressure * exp( -VolumetricHencky / SwellingSlope - 1.0);
      }   else {
         Theta = -rReferencePressure * ( -VolumetricHencky / SwellingSlope);
      }


      double Pressure = ( 1 + rAlphaShear * pow(deviatoricNorm, 2) / SwellingSlope) * Theta;

      for (unsigned int i = 0; i< 3; i++)
         rStressMatrix(i,i) += Pressure;


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

      const double & rAlphaShear = rMaterialProperties[ALPHA_SHEAR];
      const double & rYoungModulus = rMaterialProperties[YOUNG_MODULUS];
      const double & rPoissonRatio = rMaterialProperties[POISSON_RATIO];
      const double & rReferencePressure = rMaterialProperties[REFERENCE_PRESSURE];

      double BulkModulus = rYoungModulus / 3.0 / ( 1.0 - 2.0 * rPoissonRatio); 
      double SwellingSlope = rReferencePressure / BulkModulus;
      double ConstantShearModulus = rYoungModulus / 2.0 / ( 1.0 + rPoissonRatio);


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
      if ( -VolumetricHencky > SwellingSlope) {
         Phi = -SwellingSlope * rReferencePressure * exp( -VolumetricHencky / SwellingSlope - 1.0);
      }   else {
         Phi = -rReferencePressure * VolumetricHencky - rReferencePressure * pow( VolumetricHencky - SwellingSlope, 2) / 2.0 / SwellingSlope;
      }
      double Theta = 1; 
      if ( -VolumetricHencky > SwellingSlope) {
         Theta = -rReferencePressure * exp( -VolumetricHencky / SwellingSlope - 1.0);
      }   else {
         Theta = -rReferencePressure * ( -VolumetricHencky / SwellingSlope);
      }
      double K = 1; 
      if ( -VolumetricHencky > SwellingSlope) {
         K = rReferencePressure / SwellingSlope  * exp( -VolumetricHencky / SwellingSlope - 1.0);
      }   else {
         K = rReferencePressure / SwellingSlope;
      }

      // bulk modulus part
      rConstitutiveMatrix = ( 1 + rAlphaShear / SwellingSlope * pow(deviatoricNorm,2.0) ) * K * IdentityCross;

      // Shear modulus part
      rConstitutiveMatrix += 2.0*(  ConstantShearModulus + rAlphaShear / SwellingSlope * Phi) *(FourthOrderIdentity - (1.0/3.0)*IdentityCross);

      // coupling part
      Vector StrainVector = ZeroVector(6);
      StrainVector = ConstitutiveModelUtilities::StressTensorToVector( DeviatoricHencky, StrainVector); // then I do not have to divide by 2

      double Modulus = 2.0 * rAlphaShear / SwellingSlope * Theta; 

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
