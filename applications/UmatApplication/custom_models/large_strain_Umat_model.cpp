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

      Matrix AuxMatrix = rVariables.IncrementalDeformation;
      Matrix StrainMatrix;
      double det;

      MathUtils<double>::InvertMatrix( AuxMatrix, StrainMatrix, det);

      StrainMatrix = prod( trans(StrainMatrix), StrainMatrix);
      for (unsigned int i = 0; i < 3; i++)
         StrainMatrix(i,i) -= 1.0;
      StrainMatrix *= (-0.5);

      Vector StrainVector = ZeroVector(6);
      StrainVector = MathUtils<double>::StrainTensorToVector( StrainMatrix, 6);

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

      Vector PreviousStressVector = ZeroVector(6);
      for (unsigned int i = 0; i < 6; i++)
         PreviousStressVector(i) = mpStressVectorFinalized[i];

      Matrix PreviousStressTensor;
      PreviousStressTensor = MathUtils<double>::StressVectorToTensor( PreviousStressVector);
      PreviousStressTensor = prod( rVariables.IncrementalDeformation, Matrix( prod( PreviousStressTensor, trans(rVariables.IncrementalDeformation) ) ) );

      PreviousStressVector = MathUtils<double>::StressTensorToVector( PreviousStressTensor, 6);

      rpStressVector = new double[6];
      for (unsigned int i = 0; i < 6; i++)
         rpStressVector[i] = PreviousStressVector(i);


      KRATOS_CATCH("")
   }

} // Namespace Kratos
