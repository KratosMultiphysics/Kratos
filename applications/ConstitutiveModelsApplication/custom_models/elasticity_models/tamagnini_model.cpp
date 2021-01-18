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
      mSetStressState = false;
   }

   //******************************COPY CONSTRUCTOR**************************************
   //************************************************************************************

   TamagniniModel::TamagniniModel(const TamagniniModel& rOther)
      : BorjaModel( rOther ), mSetStressState(rOther.mSetStressState)
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
      mSetStressState = rOther.mSetStressState;
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

      MatrixType HenckyStrain(3,3);
      HenckyStrain  = rVariables.Strain.Matrix;

      if ( mSetStressState )
      {
         mSetStressState = false;
         SetStressState( HenckyStrain, rYoungModulus, rPoissonRatio);
      }



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

   // ***************************************************************************
   // Set Value (Vector)
   void TamagniniModel::SetValue( const Variable<Vector> & rThisVariable, const Vector & rValue, const ProcessInfo& rCurrentProcessInfo)
   {
      KRATOS_TRY

      if ( rThisVariable == ELASTIC_LEFT_CAUCHY_FROM_KIRCHHOFF_STRESS)
      {
         mInitialStressState.resize(6);
         noalias( mInitialStressState) = ZeroVector(6);
         mInitialStressState = rValue;
         mSetStressState = true;

      }

      KRATOS_CATCH("")
   }

   // *********************************************************************************
   // Set Stress state
   void TamagniniModel::SetStressState( MatrixType & rHenckyStrain, const double & rE, const double & rNu)
   {

      KRATOS_TRY
      Matrix StressMat(3,3);
      noalias( StressMat) = ZeroMatrix(3,3);

      for (int i = 0; i < 3; ++i)
         StressMat(i,i) = mInitialStressState(i);


      Vector EigenStress(3);
      MatrixType EigenStressM;
      MatrixType  EigenV;
      noalias( EigenStressM ) = ZeroMatrix(3,3);
      noalias( EigenV) = ZeroMatrix(3,3);
      MatrixType OriginalHencky;
      noalias( OriginalHencky) = ZeroMatrix(3,3);
      OriginalHencky = rHenckyStrain;

   MathUtils<double>::GaussSeidelEigenSystem<MatrixType, MatrixType>(StressMat, EigenV, EigenStressM);
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

   noalias(ElasticLeftCauchy) = prod(EigenV, ElasticLeftCauchy);
   noalias(ElasticLeftCauchy) = prod(ElasticLeftCauchy, trans(EigenV));

   noalias(rHenckyStrain) = prod(EigenV, rHenckyStrain);
   noalias(rHenckyStrain) = prod(rHenckyStrain, trans(EigenV));

      rHenckyStrain += OriginalHencky;

      Vector NewHistoryVector(6);
      NewHistoryVector = ConstitutiveModelUtilities::StrainTensorToVector( ElasticLeftCauchy, NewHistoryVector);
      this->mHistoryVector = NewHistoryVector;




      KRATOS_CATCH("")
   }


} // Namespace Kratos
