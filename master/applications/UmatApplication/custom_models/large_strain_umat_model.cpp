//
//   Project Name:        KratosUmatApplication        $
//   Created by:          $Author:           LMonforte $
//   Last modified by:    $Co-Author:                  $
//   Date:                $Date:          October 2017 $
//   Revision:            $Revision:               0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_models/large_strain_umat_model.hpp"



namespace Kratos
{

   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   LargeStrainUmatModel::LargeStrainUmatModel()
      : SmallStrainUmatModel()
   {
      mInitializedModel = false;
   }

   //******************************COPY CONSTRUCTOR**************************************
   //************************************************************************************

   LargeStrainUmatModel::LargeStrainUmatModel(const LargeStrainUmatModel& rOther)
      : SmallStrainUmatModel(rOther)
   {
   }

   //********************************CLONE***********************************************
   //************************************************************************************

   ConstitutiveModel::Pointer LargeStrainUmatModel::Clone() const
   {
     return Kratos::make_shared<LargeStrainUmatModel>(*this);
   }

   //********************************ASSIGNMENT******************************************
   //************************************************************************************
   LargeStrainUmatModel& LargeStrainUmatModel::operator=(LargeStrainUmatModel const& rOther)
   {
      SmallStrainUmatModel::operator=(rOther);
      return *this;
   }

   //*******************************DESTRUCTOR*******************************************
   //************************************************************************************

   LargeStrainUmatModel::~LargeStrainUmatModel()
   {
   }


   // ************************* CREATE STRAIN VECTORS **********************************
   //***********************************************************************************
   void LargeStrainUmatModel::CreateStrainsVectors( UmatDataType & rVariables,  double* & rpStrain, double* & rpIncrementalStrain)
   {
      KRATOS_TRY

      MatrixType StrainMatrix;
      double det;
      ConstitutiveModelUtilities::InvertMatrix3( rVariables.IncrementalDeformation, StrainMatrix, det);

      StrainMatrix = prod( trans(StrainMatrix), StrainMatrix);
      for (unsigned int i = 0; i < 3; i++)
         StrainMatrix(i,i) -= 1.0;
      StrainMatrix *= (-0.5);

      VectorType StrainVector;
      ConstitutiveModelUtilities::StrainTensorToVector(StrainMatrix, StrainVector);

      rpIncrementalStrain = new double[6];
      rpStrain = new double[6];
      for (unsigned int i = 0; i < 6; i++) {
         rpIncrementalStrain[i] = StrainVector(i);
         rpStrain[i] = 0.0;
      }


      KRATOS_CATCH("")
   }


   // ****************************** TIME n STRESS STATE *******************************
   // **********************************************************************************
   void LargeStrainUmatModel::CreateStressAtInitialState( UmatDataType & rVariables, double* & rpStressVector)
   {
      KRATOS_TRY


      MatrixType PreviousStressTensor;
      PreviousStressTensor = ConstitutiveModelUtilities::StressVectorToTensor( mStressVectorFinalized, PreviousStressTensor);
      PreviousStressTensor = prod( rVariables.IncrementalDeformation, Matrix( prod( PreviousStressTensor, trans(rVariables.IncrementalDeformation) ) ) );

      Vector PreviousStressVector(6);
      PreviousStressVector = ConstitutiveModelUtilities::StressTensorToVector(PreviousStressTensor, PreviousStressVector);

      rpStressVector = new double[6];
      for (unsigned int i = 0; i < 6; i++)
         rpStressVector[i] = PreviousStressVector(i);


      KRATOS_CATCH("")
   }


   //***********************
   //************************************************************************************
   //************************************************************************************
   void LargeStrainUmatModel::SetConstitutiveMatrix( Matrix & rC, const Matrix & rpCBig, const MatrixType & rStressMatrix)
   {
      KRATOS_TRY

      Matrix ExtraMatrix = ZeroMatrix(6);
      ExtraMatrix = CalculateExtraMatrix( rStressMatrix, ExtraMatrix);
      for (unsigned int i = 0; i < 6; i++){
         for (unsigned int j = 0; j < 6; j++){
            ExtraMatrix(i,j) += rpCBig(i,j);
         }
      }

      SmallStrainUmatModel::SetConstitutiveMatrix( rC, ExtraMatrix, rStressMatrix);

      KRATOS_CATCH("")
   }
   // ************************** COMPUTE MORE CONSTITUTIVE TERMS *********************
   // ********************************************************************************
   // ****** SOME SORT OF TERMS
   Matrix& LargeStrainUmatModel::CalculateExtraMatrix( const MatrixType& rStressMatrix, Matrix& rExtraMatrix)
   {
      KRATOS_TRY

      // NO SE M'ACUT UNA MANERA MÉS GUARRA DE FER AIXÒ.
      Matrix ExtraMatrix = ZeroMatrix(6,6);
      Matrix Identity = identity_matrix<double>(3);

      unsigned indexi, indexj;
      for ( unsigned int i = 0; i < 3; i++) {
         for (unsigned int j = 0; j < 3; j++) {
            for (unsigned int k = 0; k < 3; k++) {
               for (unsigned int l = 0; l < 3; l++) {
                  unsigned int auxi, auxj;
                  if ( i < j) {
                     auxi = i;
                     auxj = j;
                  } else {
                     auxi = j;
                     auxj = i;
                  }
                  unsigned int auxk, auxl;
                  if ( k < l) {
                     auxk = k;
                     auxl = l;
                  } else {
                     auxk = l;
                     auxl = k;
                  }

                  // now I have to decide where i put and what I put
                  if ( auxi == auxj) {
                     indexi = auxi;
                  }
                  else if ( auxi ==0) {
                     if ( auxj == 1) {
                        indexi = 3;
                     }
                     else {
                        indexi = 4;
                     }
                  }
                  else {
                     indexi = 5;
                  }
                  if (auxk == auxl) {
                     indexj = auxk;
                  }
                  else if ( auxk ==0) {
                     if ( auxl == 1) {
                        indexj = 3;
                     }
                     else {
                        indexj = 4;
                     }
                  }
                  else {
                     indexj = 5;
                  }
                  double voigtnumber = 1.0;
                  if ( indexj > 2)
                     voigtnumber *= 0.50;
                  if ( indexi > 2)
                     voigtnumber *=0.5;
                  ExtraMatrix(indexi, indexj) -= voigtnumber * (Identity(i,k) * rStressMatrix(j,l) +  Identity(j,k)*rStressMatrix(i,l) ); // vale, molt bé, però no sé on queda la "voigt notation"

               }
            }
         }
      }


      noalias(rExtraMatrix) = ExtraMatrix;
      return rExtraMatrix;

      KRATOS_CATCH("")
   }

} // Namespace Kratos
