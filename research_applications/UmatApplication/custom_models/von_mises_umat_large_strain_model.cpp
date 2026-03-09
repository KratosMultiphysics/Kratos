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
#include "custom_models/von_mises_umat_large_strain_model.hpp"



namespace Kratos
{

   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   VonMisesLargeStrainUmatModel::VonMisesLargeStrainUmatModel()
      : LargeStrainUmatModel()
   {
   }

   //******************************COPY CONSTRUCTOR**************************************
   //************************************************************************************

   VonMisesLargeStrainUmatModel::VonMisesLargeStrainUmatModel(const VonMisesLargeStrainUmatModel& rOther)
      : LargeStrainUmatModel(rOther)
   {
   }

   //********************************CLONE***********************************************
   //************************************************************************************

   ConstitutiveModel::Pointer VonMisesLargeStrainUmatModel::Clone() const
   {
     return Kratos::make_shared<VonMisesLargeStrainUmatModel>(*this);
   }

   //********************************ASSIGNMENT******************************************
   //************************************************************************************
   VonMisesLargeStrainUmatModel& VonMisesLargeStrainUmatModel::operator=(VonMisesLargeStrainUmatModel const& rOther)
   {
      LargeStrainUmatModel::operator=(rOther);
      return *this;
   }

   //*******************************DESTRUCTOR*******************************************
   //************************************************************************************

   VonMisesLargeStrainUmatModel::~VonMisesLargeStrainUmatModel()
   {
   }


} // Namespace Kratos
