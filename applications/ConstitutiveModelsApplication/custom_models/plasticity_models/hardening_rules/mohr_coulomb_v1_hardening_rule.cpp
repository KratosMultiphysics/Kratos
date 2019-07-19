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
#include "custom_models/plasticity_models/hardening_rules/mohr_coulomb_v1_hardening_rule.hpp"

namespace Kratos
{

   //*******************************CONSTRUCTOR******************************************
   //************************************************************************************

   MohrCoulombV1HardeningRule::MohrCoulombV1HardeningRule()
      :HardeningRule()
   {
   }


   //*******************************ASSIGMENT OPERATOR***********************************
   //************************************************************************************

   MohrCoulombV1HardeningRule& MohrCoulombV1HardeningRule::operator=(MohrCoulombV1HardeningRule const& rOther)
   {
      HardeningRule::operator=(rOther);
      return *this;
   }

   //*******************************COPY CONSTRUCTOR*************************************
   //************************************************************************************

   MohrCoulombV1HardeningRule::MohrCoulombV1HardeningRule(MohrCoulombV1HardeningRule const& rOther)
      :HardeningRule(rOther)
   {

   }


   //********************************CLONE***********************************************
   //************************************************************************************

   HardeningRule::Pointer MohrCoulombV1HardeningRule::Clone() const
   {
      return Kratos::make_shared<MohrCoulombV1HardeningRule>(*this);
   }


   //********************************DESTRUCTOR******************************************
   //************************************************************************************

   MohrCoulombV1HardeningRule::~MohrCoulombV1HardeningRule()
   {
   }

   /// Operations.

   //*******************************CALCULATE TOTAL HARDENING****************************
   //************************************************************************************

   double& MohrCoulombV1HardeningRule::CalculateHardening(const PlasticDataType& rVariables, double& rHardening)
   {
      KRATOS_TRY


      rHardening = 0;
      return rHardening;

      KRATOS_CATCH(" ")
   }



   //*******************************CALCULATE HARDENING DERIVATIVE***********************
   //************************************************************************************

   double& MohrCoulombV1HardeningRule::CalculateDeltaHardening(const PlasticDataType& rVariables, double& rDeltaHardening)
   {
      KRATOS_TRY

      rDeltaHardening = 0;
      return rDeltaHardening;

      KRATOS_CATCH(" ")
   }

   //*******************************CALCULATE HARDENING DERIVATIVE***********************
   //************************************************************************************
   double& MohrCoulombV1HardeningRule::CalculateDeltaHardening(const PlasticDataType& rVariables, double& rDeltaHardening, const MatrixType & rPlasticPotentialDerivative)
   {
      KRATOS_TRY

      rDeltaHardening = 0;
      return rDeltaHardening;
      KRATOS_CATCH(" ")
   }
}  // namespace Kratos.
