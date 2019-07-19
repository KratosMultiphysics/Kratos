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
      const Properties    & rProperties = rModelData.GetProperties();

      const double & rSwellingSlope = rProperties[SWELLING_SLOPE];
      const double & rAlphaShear = rProperties[ALPHA_SHEAR];
      const double & rReferencePressure = rProperties[PRE_CONSOLIDATION_STRESS];
      const double & rOCR = rProperties[OVER_CONSOLIDATION_RATIO];
      const double ReferencePressure = rReferencePressure / rOCR;
      const double & rConstantShearModulus = rProperties[INITIAL_SHEAR_MODULUS];

      MatrixType HenckyStrain(3,3);
      HenckyStrain = rVariables.Strain.Matrix;

      if ( this->mSetStressState )
      {
         this->mSetStressState = false;
         SetStressState( HenckyStrain, ReferencePressure, rSwellingSlope, rAlphaShear, rConstantShearModulus);
      }

      // 2.a Separate Volumetric and deviatoric part
      double VolumetricHencky;
      MatrixType DeviatoricHencky(3,3);
      noalias(DeviatoricHencky) = ZeroMatrix(3,3);
      double deviatoricNorm;

      SeparateVolumetricAndDeviatoricPart( HenckyStrain, VolumetricHencky, DeviatoricHencky, deviatoricNorm);

      // 3.a Compute Deviatoric Part
      rStressMatrix.clear();
      double ShearModulus = rAlphaShear * ReferencePressure * std::exp( -VolumetricHencky / rSwellingSlope );
      rStressMatrix += DeviatoricHencky * 2.0 * ( ShearModulus  + rConstantShearModulus );

      // 3.b Compute Volumetric Part
      double pressure = -ReferencePressure * std::exp( -VolumetricHencky / rSwellingSlope ) * ( 1.0 + rAlphaShear * pow(deviatoricNorm, 2) / rSwellingSlope);

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
      const Properties    & rProperties = rModelData.GetProperties();

      const double & rSwellingSlope = rProperties[SWELLING_SLOPE];
      const double & rAlphaShear = rProperties[ALPHA_SHEAR];
      const double & rReferencePressure = rProperties[PRE_CONSOLIDATION_STRESS];
      const double & rOCR = rProperties[OVER_CONSOLIDATION_RATIO];
      const double ReferencePressure = rReferencePressure / rOCR;
      const double & rConstantShearModulus = rProperties[INITIAL_SHEAR_MODULUS];

      // 1. Define some matrices
      Matrix FourthOrderIdentity(6,6);
      noalias( FourthOrderIdentity ) = ZeroMatrix(6,6);
      for (unsigned int i = 0; i<3; ++i)
         FourthOrderIdentity(i,i) = 1.0;

      for (unsigned int i = 3; i<6; ++i)
         FourthOrderIdentity(i,i) = 0.50;

      Matrix IdentityCross(6,6);
      noalias( IdentityCross )  = ZeroMatrix(6,6);
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
      double pressure = -ReferencePressure * std::exp( -VolumetricHencky / rSwellingSlope ) * ( 1.0 + rAlphaShear * pow(deviatoricNorm, 2) / rSwellingSlope);
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


   // *********************************************************************************
   // Set Stress state
   void BorjaModel::SetStressState( MatrixType & rHenckyStrain, const double & rReferencePressure, const double & rSwellingSlope, const double & rAlphaShear, const double & rConstantShearModulus)
   {

      KRATOS_TRY

      MatrixType OriginalHencky;
      noalias( OriginalHencky) = ZeroMatrix(3,3);
      OriginalHencky = rHenckyStrain;

      Vector Objective(3);
      Objective(0) = this->mInitialStressState(0) + this->mInitialStressState(1) + this->mInitialStressState(2);
      Objective(0) /= 3.0;
      Objective(1) = this->mInitialStressState(1) - Objective(0);
      Objective(2) = this->mInitialStressState(2) - Objective(0);

      Vector Guess = ZeroVector(3);

      bool NotConverged = true;

      Vector Y = ZeroVector(3);

      Matrix TangentMatrix = ZeroMatrix(3,3);
      Matrix InverseTangent = ZeroMatrix(3,3);
      Vector Residual = ZeroVector(3);
      Vector dGuess = Residual;
      double DeviatoricNorm2;
      double error, detI, ShearModulus; 
      unsigned nIter = 0;

      while (NotConverged) {

         //1 COMPUTE SOME ERROR
         DeviatoricNorm2 = pow( Guess(1), 2) + 2.0*pow( Guess(2), 2);
         ShearModulus = rAlphaShear*rReferencePressure*std::exp( -Guess(0) / rSwellingSlope) + rConstantShearModulus;
         Y(0) = -rReferencePressure *std::exp( - Guess(0) / rSwellingSlope) * ( 1.0 + rAlphaShear/rSwellingSlope * DeviatoricNorm2) ;
         Y(1) = 2.0*ShearModulus*Guess(1);
         Y(2) = 2.0*ShearModulus*Guess(2);


         Residual = Y - Objective;
         error = 0.0;
         for (unsigned int i = 0; i < 3; ++i)
            error += pow( Residual(i), 2);

         if (error < 1.0e-12) {
            NotConverged = false;
         }

         //1.1 Compute the Tangent Matrix
         TangentMatrix(0,0) =  rReferencePressure*std::exp( -Guess(0)/rSwellingSlope) * (1.0 + rAlphaShear/rSwellingSlope*DeviatoricNorm2) / rSwellingSlope;
         TangentMatrix(0,1) = -rReferencePressure*std::exp( -Guess(0)/rSwellingSlope) * ( rAlphaShear/rSwellingSlope) * 2.0 *  Guess(1) ;
         TangentMatrix(0,2) = -rReferencePressure*std::exp( -Guess(0)/rSwellingSlope) * ( rAlphaShear/rSwellingSlope) * 2.0 *  Guess(2) * 2.0 ;


         TangentMatrix(1,0) =  2.0* rAlphaShear * rReferencePressure*std::exp(-Guess(0) / rSwellingSlope)*(-1.0/rSwellingSlope)*Guess(1);
         TangentMatrix(1,1) = 2.0*ShearModulus;
         TangentMatrix(1,2) = 0.0;

         TangentMatrix(2,0) =  2.0* rAlphaShear * rReferencePressure*std::exp(-Guess(0) / rSwellingSlope)*(-1.0/rSwellingSlope)*Guess(2);
         TangentMatrix(2,1) = 0.0;
         TangentMatrix(2,2) = 2.0*ShearModulus;

         MathUtils<double>::InvertMatrix( TangentMatrix, InverseTangent, detI);

         dGuess = prod( InverseTangent, Residual);

         // (Try to solve some convergence problems)
         double dGNorm = 0.0;
         for (unsigned int i = 0; i <3 ; ++i)
            dGNorm += pow( dGuess(i), 2);
         dGNorm = sqrt( dGNorm);
         if ( dGNorm > 0.01)
            dGuess *= 0.1;


         Guess -= dGuess;
         nIter += 1;

         if (nIter > 1000) {
            std::cout << " NONCONVERGING::Initial Stress State of BORJA " << std::endl;
            return;
         }

      }

      // 2. Set the ElasticLeftCauchy
      double Hencky1 = Guess(1) + Guess(0)/3.0;
      double Hencky2 = Guess(2) + Guess(0)/3.0; 

      this->mHistoryVector = ZeroVector(6);

      this->mHistoryVector(1) = std::exp( 2.0*Hencky1);

      this->mHistoryVector(0) = std::exp( 2.0*Hencky2);
      this->mHistoryVector(2) = std::exp( 2.0*Hencky2);

      noalias( rHenckyStrain ) = ZeroMatrix(3,3);

      rHenckyStrain(1,1) = Hencky1;

      rHenckyStrain(0,0) = Hencky2;
      rHenckyStrain(2,2) = Hencky2;

      rHenckyStrain += OriginalHencky;

      KRATOS_CATCH("")
   }
} // Namespace Kratos
