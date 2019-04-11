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
#include "custom_models/sanisand_small_strain_model.hpp"



namespace Kratos
{

   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   SanisandSmallStrainUmatModel::SanisandSmallStrainUmatModel()
      : SmallStrainUmatModel()
   {
   }

   //******************************COPY CONSTRUCTOR**************************************
   //************************************************************************************

   SanisandSmallStrainUmatModel::SanisandSmallStrainUmatModel(const SanisandSmallStrainUmatModel& rOther)
      : SmallStrainUmatModel(rOther)
   {
   }

   //********************************CLONE***********************************************
   //************************************************************************************

   ConstitutiveModel::Pointer SanisandSmallStrainUmatModel::Clone() const
   {
     return Kratos::make_shared<SanisandSmallStrainUmatModel>(*this);
   }

   //********************************ASSIGNMENT******************************************
   //************************************************************************************
   SanisandSmallStrainUmatModel& SanisandSmallStrainUmatModel::operator=(SanisandSmallStrainUmatModel const& rOther)
   {
      SmallStrainUmatModel::operator=(rOther);
      return *this;
   }

   //*******************************DESTRUCTOR*******************************************
   //************************************************************************************

   SanisandSmallStrainUmatModel::~SanisandSmallStrainUmatModel()
   {
   }


} // Namespace Kratos
