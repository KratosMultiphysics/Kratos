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
#include "custom_models/elasticity_models/hencky_linear_model.hpp"


namespace Kratos
{

   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   HenckyLinearModel::HenckyLinearModel()
      : HenckyHyperElasticModel()
   {
   }

   //******************************COPY CONSTRUCTOR**************************************
   //************************************************************************************

   HenckyLinearModel::HenckyLinearModel(const HenckyLinearModel& rOther)
      : HenckyHyperElasticModel( rOther )
   {
   }

   //********************************CLONE***********************************************
   //************************************************************************************

   ConstitutiveModel::Pointer HenckyLinearModel::Clone() const
   {
      return Kratos::make_shared<HenckyLinearModel>(*this);
   }

   //********************************ASSIGNMENT******************************************
   //************************************************************************************
   HenckyLinearModel& HenckyLinearModel::operator=(HenckyLinearModel const& rOther)
   {
      HenckyHyperElasticModel::operator=(rOther);
      return *this;
   }

   //*******************************DESTRUCTOR*******************************************
   //************************************************************************************

   HenckyLinearModel::~HenckyLinearModel()
   {
   }

   // ****************************************************************
   // ***************************************************************
   // not in the correct place
   void HenckyLinearModel::SeparateVolumetricAndDeviatoricPart( const MatrixType& rA, double & rVolumetric, MatrixType& rDev, double & devNorm)
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
   void HenckyLinearModel::CalculateAndAddStressTensor( HyperElasticDataType& rVariables, MatrixType& rStressMatrix)
   {
      KRATOS_TRY

      // material parameters

      const ModelDataType & rModelData = rVariables.GetModelData();
      const Properties    & rMaterialProperties = rModelData.GetProperties();

      const double & rYoungModulus = rMaterialProperties[YOUNG_MODULUS];
      const double & rPoissonRatio = rMaterialProperties[POISSON_RATIO];

      double BulkModulus = rYoungModulus / 3.0 / ( 1.0 - 2.0 * rPoissonRatio);
      double ShearModulus = rYoungModulus / 2.0 / (1.0 + rPoissonRatio);


      MatrixType HenckyStrain(3,3);
      HenckyStrain = rVariables.Strain.Matrix;

      if ( this->mSetStressState )
      {
         this->mSetStressState = false;
         SetStressState( HenckyStrain, rYoungModulus, rPoissonRatio);
      }

      // 2.a Separate Volumetric and deviatoric part
      double VolumetricHencky;
      MatrixType DeviatoricHencky(3,3);
      noalias(DeviatoricHencky) = ZeroMatrix(3,3);
      double deviatoricNorm;

      SeparateVolumetricAndDeviatoricPart( HenckyStrain, VolumetricHencky, DeviatoricHencky, deviatoricNorm);

      // 3.a Compute Deviatoric Part
      rStressMatrix.clear();
      rStressMatrix += DeviatoricHencky * 2.0 * ShearModulus;

      // 3.b Compute Volumetric Part
      double pressure = VolumetricHencky * BulkModulus;

      for (unsigned int i = 0; i< 3; i++)
         rStressMatrix(i,i) += pressure;

      KRATOS_CATCH("")
   }


   //************************************************************************************
   //************************************************************************************
   // CalculateAndAddConstitutiveTensor
   void HenckyLinearModel::CalculateAndAddConstitutiveTensor( HyperElasticDataType& rVariables, Matrix& rConstitutiveMatrix)
   {
      KRATOS_TRY

      // material parameters
      const ModelDataType & rModelData = rVariables.GetModelData();
      const Properties    & rMaterialProperties = rModelData.GetProperties();

      const double & rYoung = rMaterialProperties[YOUNG_MODULUS];
      const double & rNu = rMaterialProperties[POISSON_RATIO];

      double diagonal =   rYoung/(1.0+rNu)/(1.0-2.0*rNu) * (1.0-rNu);
      double nodiagonal = rYoung/(1.0+rNu)/(1.0-2.0*rNu) * ( rNu);
      double corte      = rYoung/(1.0+rNu)/2.0;

      rConstitutiveMatrix.clear();

      for (unsigned int i = 0; i<3; ++i) {
         for (unsigned int j = 0; j<3; ++j) {
            if (i == j) {
               rConstitutiveMatrix(i,i) = diagonal;
            }
            else {
               rConstitutiveMatrix(i,j) = nodiagonal;
            }
         }
      }

      for (unsigned int j = 3; j<6; ++j)
         rConstitutiveMatrix(j,j) = corte;


      rVariables.State().Set(ConstitutiveModelData::CONSTITUTIVE_MATRIX_COMPUTED);

      KRATOS_CATCH("")
   }


   // *********************************************************************************
   // Set Stress state
   void HenckyLinearModel::SetStressState( MatrixType & rHenckyStrain, const double & rE, const double & rNu)
   {

      KRATOS_TRY
      Matrix StressMat(3,3);
      noalias( StressMat) = ZeroMatrix(3,3);

      for (int i = 0; i < 3; ++i)
         StressMat(i,i) = this->mInitialStressState(i);


      Vector EigenStress(3);
      MatrixType EigenStressM;
      MatrixType  EigenV;
      noalias( EigenStressM ) = ZeroMatrix(3,3);
      noalias( EigenV) = ZeroMatrix(3,3);
      MatrixType OriginalHencky;
      noalias( OriginalHencky) = ZeroMatrix(3,3);
      OriginalHencky = rHenckyStrain;

      //SolidMechanicsMathUtilities<double>::EigenVectors( StressMat, EigenV, EigenStress);
      MathUtils<double>::GaussSeidelEigenSystem( StressMat, EigenV, EigenStressM);
      for (unsigned int i = 0; i < 3; i++)
         EigenStress(i) = EigenStressM(i,i);


      Vector ElasticHenckyStrain(3);
      Matrix InverseElastic(3,3);
      const double & YoungModulus = rE;
      const double & PoissonCoef  = rNu;


      for (unsigned int i = 0; i < 3; ++i) {
         for (unsigned int j = 0; j < 3; ++j) {
            if (i == j ) {
               InverseElastic(i,i) = 1.0/YoungModulus;
            }
            else {
               InverseElastic(i,j) = - PoissonCoef / YoungModulus;
            }
         }
      }

      noalias( ElasticHenckyStrain )  = prod( InverseElastic, EigenStress);


      MatrixType ElasticLeftCauchy(3,3);
      noalias( ElasticLeftCauchy ) = ZeroMatrix(3,3);
      noalias( rHenckyStrain ) = ZeroMatrix(3,3);

      for (unsigned int i = 0; i < 3; ++i) {
         ElasticLeftCauchy(i,i) = std::exp(2.0*ElasticHenckyStrain(i));
         rHenckyStrain(i,i) = ElasticHenckyStrain(i);
      }

      noalias( ElasticLeftCauchy ) = prod( EigenV, ElasticLeftCauchy);
      noalias( ElasticLeftCauchy )  = prod( ElasticLeftCauchy, trans(EigenV) );

      noalias( rHenckyStrain ) = prod( EigenV, rHenckyStrain);
      noalias( rHenckyStrain ) = prod( rHenckyStrain, trans(EigenV) );

      rHenckyStrain += OriginalHencky;

      Vector NewHistoryVector(6);
      NewHistoryVector = ConstitutiveModelUtilities::StrainTensorToVector( ElasticLeftCauchy, NewHistoryVector);
      this->mHistoryVector = NewHistoryVector;


      KRATOS_CATCH("")
   }


} // Namespace Kratos
