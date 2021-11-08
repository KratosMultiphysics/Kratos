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
      pVector[0] = rMaterialProperties[YOUNG_MODULUS]; // E 
      pVector[1] = rMaterialProperties[POISSON_RATIO]; // nu
      pVector[2] = rMaterialProperties[ALPHA]; // alpha
      pVector[3] = rMaterialProperties[BETA]; // beta
      pVector[4] = rMaterialProperties[MF]; // Mf
      pVector[5] = rMaterialProperties[CC]; // cc
      pVector[6] = rMaterialProperties[MM]; // mm
      pVector[7] = rMaterialProperties[RHOS]; // rhos
      pVector[8] = rMaterialProperties[KSIS]; // ksis
      pVector[9] = rMaterialProperties[RHOM]; // rhom
      pVector[10] =rMaterialProperties[KSIM]; // ksim
      pVector[11] =rMaterialProperties[PC0]; // pc0

      KRATOS_CATCH("")
   }

   void FabricSmallStrainUmatModel::InitializeStateVariables( Vector & rStateVariables, const Properties & rMaterialProperties)
   {
      KRATOS_TRY

      rStateVariables(0) = 6.92; // void ratio
      rStateVariables(7) = -1000.0; // PS
      rStateVariables(8) = -1000.0; // PM
      rStateVariables(9) = -2000.0; // PS+PM == PC

      rStateVariables(0) = rMaterialProperties[VOID_RATIO]; // void ratio
      rStateVariables(7) = rMaterialProperties[PS]; // PS
      rStateVariables(8) = rMaterialProperties[PM]; // PM
      rStateVariables(9) = rMaterialProperties[PLASTIC_MULTIPLIER]; // plastic multiplier

      KRATOS_CATCH("")
   }

   double& FabricSmallStrainUmatModel::GetValue(const Variable<double> & rVariable, double& rValue)
   {
      KRATOS_TRY
     
      if ( rVariable == PS)
      {
         rValue = mStateVariablesFinalized[7];
      } else if ( rVariable == PM ) {
         rValue = mStateVariablesFinalized[8];
      } else if ( rVariable == VOID_RATIO) {
         rValue = mStateVariablesFinalized[9];
      } else if ( rVariable == PLASTIC_MULTIPLIER) {
         rValue = mStateVariablesFinalized[0];
      } else {
         rValue = SmallStrainUmatModel::GetValue( rVariable, rValue);
      }

      return rValue;
      
      KRATOS_CATCH("")
   }
} // Namespace Kratos
