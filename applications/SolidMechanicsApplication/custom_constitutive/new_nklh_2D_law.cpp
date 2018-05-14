//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:               DAbadias $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2018 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes
#include <iostream>

// External includes
#include<cmath>

// Project includes
#include "includes/properties.h"
#include "custom_constitutive/new_nklh_2D_law.hpp"

#include "solid_mechanics_application_variables.h"

namespace Kratos
{

   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   NewNKLH2DLaw::NewNKLH2DLaw()
      : NewNKLH3DLaw()
   {
   }

   //******************************COPY CONSTRUCTOR**************************************
   //************************************************************************************

   NewNKLH2DLaw::NewNKLH2DLaw(const NewNKLH2DLaw& rOther)
      : NewNKLH3DLaw(rOther)
   {
   }

   //********************************CLONE***********************************************
   //************************************************************************************

   ConstitutiveLaw::Pointer NewNKLH2DLaw::Clone() const
   {
      NewNKLH2DLaw::Pointer p_clone(new NewNKLH2DLaw(*this));
      return p_clone;
   }

   //*******************************DESTRUCTOR*******************************************
   //************************************************************************************

   NewNKLH2DLaw::~NewNKLH2DLaw()
   {
   }
   
   void NewNKLH2DLaw::GetLawFeatures(Features& rFeatures)
   {
      //Set the type of law
      rFeatures.mOptions.Set( AXISYMMETRIC_LAW );
      rFeatures.mOptions.Set( INFINITESIMAL_STRAINS );
      rFeatures.mOptions.Set( ISOTROPIC );

      //Set strain measure required by the consitutive law
      rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
      rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

      //Set the strain size
      rFeatures.mStrainSize = GetStrainSize();

      //Set the spacedimension
      rFeatures.mSpaceDimension = WorkingSpaceDimension();

   }

} // Namespace Kratos
