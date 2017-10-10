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
#include "custom_models/large_strain_Umat_model.hpp"



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
      return ( LargeStrainUmatModel::Pointer(new LargeStrainUmatModel(*this)) );
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
      StrainVector = ConstitutiveModelUtilities::StrainTensorToVector(StrainMatrix, StrainVector);

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

} // Namespace Kratos
