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
#include "custom_models/von_mises_umat_small_strain_model.hpp"



namespace Kratos
{

   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   VonMisesSmallStrainUmatModel::VonMisesSmallStrainUmatModel()
      : SmallStrainUmatModel()
   {
   }

   //******************************COPY CONSTRUCTOR**************************************
   //************************************************************************************

   VonMisesSmallStrainUmatModel::VonMisesSmallStrainUmatModel(const VonMisesSmallStrainUmatModel& rOther)
      : SmallStrainUmatModel(rOther)
   {
   }

   //********************************CLONE***********************************************
   //************************************************************************************

   ConstitutiveModel::Pointer VonMisesSmallStrainUmatModel::Clone() const
   {
     return Kratos::make_shared<VonMisesSmallStrainUmatModel>(*this);
   }

   //********************************ASSIGNMENT******************************************
   //************************************************************************************
   VonMisesSmallStrainUmatModel& VonMisesSmallStrainUmatModel::operator=(VonMisesSmallStrainUmatModel const& rOther)
   {
      SmallStrainUmatModel::operator=(rOther);
      return *this;
   }

   //*******************************DESTRUCTOR*******************************************
   //************************************************************************************

   VonMisesSmallStrainUmatModel::~VonMisesSmallStrainUmatModel()
   {
   }


} // Namespace Kratos
