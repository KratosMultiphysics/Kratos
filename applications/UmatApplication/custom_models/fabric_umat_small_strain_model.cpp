//
//   Project Name:        KratosUmatApplication        $
//   Created by:          $Author:           LMonforte $
//   Last modified by:    $Co-Author:                  $
//   Date:                $Date:            April 2018 $
//   Revision:            $Revision:               0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_models/fabric_umat_small_strain_model.hpp"



namespace Kratos
{

   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   FabricSmallStrainUmatModel::FabricSmallStrainUmatModel()
      : SmallStrainUmatModel()
   {
   }

   //******************************COPY CONSTRUCTOR**************************************
   //************************************************************************************

   FabricSmallStrainUmatModel::FabricSmallStrainUmatModel(const FabricSmallStrainUmatModel& rOther)
      : SmallStrainUmatModel(rOther)
   {
   }

   //********************************CLONE***********************************************
   //************************************************************************************

   ConstitutiveModel::Pointer FabricSmallStrainUmatModel::Clone() const
   {
      return ( FabricSmallStrainUmatModel::Pointer(new FabricSmallStrainUmatModel(*this)) );
   }

   //********************************ASSIGNMENT******************************************
   //************************************************************************************
   FabricSmallStrainUmatModel& FabricSmallStrainUmatModel::operator=(FabricSmallStrainUmatModel const& rOther)
   {
      SmallStrainUmatModel::operator=(rOther);
      return *this;
   }

   //*******************************DESTRUCTOR*******************************************
   //************************************************************************************

   FabricSmallStrainUmatModel::~FabricSmallStrainUmatModel()
   {
   }



   // ***********************************************************************************
   // more functions. Constitutive parameters (temporary solution)
   void FabricSmallStrainUmatModel::CreateConstitutiveParametersVector(double* & pVector, int & rNumberParameters, const Properties & rMaterialProperties) 
   {
      KRATOS_TRY

      rNumberParameters = 12;

      pVector = new double[rNumberParameters];
      pVector[0] = 42.0E3; // E 
      pVector[1] = 0.17; // nu
      pVector[2] = 1.0; // alpha
      pVector[3] = 1.0; // beta
      pVector[4] = 2.2; // Mf
      pVector[5] = 0.8; // cc
      pVector[6] = -0.5; // mm
      pVector[7] = 5.0; // rhos
      pVector[8] = -1.0; // ksis
      pVector[9] = 10.0; // rhom
      pVector[10] = -2.0; // ksim
      pVector[11] = 500.0; // pc0

      KRATOS_CATCH("")
   }

   // **********************************************************************************
   // set initial stress state (very temporary solution)
   void FabricSmallStrainUmatModel::SetInitialStressTEMPORARY( VectorType & rStressVector) 
   {
      KRATOS_TRY

      rStressVector(0) = -20.0;
      rStressVector(1) = -20.0;
      rStressVector(2) = -20.0;

      KRATOS_CATCH("")
   }

   void FabricSmallStrainUmatModel::InitializeStateVariables( Vector & rStateVariables)
   {
      KRATOS_TRY

      rStateVariables(0) = 6.92; // void ratio
      rStateVariables(7) = -1000.0; // PS
      rStateVariables(8) = -1000.0; // PM
      rStateVariables(9) = -2000.0; // plastic multiplier

      KRATOS_CATCH("")
   }

} // Namespace Kratos
