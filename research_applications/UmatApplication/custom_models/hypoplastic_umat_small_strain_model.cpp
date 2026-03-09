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
#include "custom_models/hypoplastic_umat_small_strain_model.hpp"



namespace Kratos
{

   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   HypoplasticSmallStrainUmatModel::HypoplasticSmallStrainUmatModel()
      : SmallStrainUmatModel()
   {

   }

   //******************************COPY CONSTRUCTOR**************************************
   //************************************************************************************

   HypoplasticSmallStrainUmatModel::HypoplasticSmallStrainUmatModel(const HypoplasticSmallStrainUmatModel& rOther)
      : SmallStrainUmatModel(rOther)
   {
   }

   //********************************CLONE***********************************************
   //************************************************************************************

   ConstitutiveModel::Pointer HypoplasticSmallStrainUmatModel::Clone() const
   {
     return Kratos::make_shared<HypoplasticSmallStrainUmatModel>(*this);
   }

   //********************************ASSIGNMENT******************************************
   //************************************************************************************
   HypoplasticSmallStrainUmatModel& HypoplasticSmallStrainUmatModel::operator=(HypoplasticSmallStrainUmatModel const& rOther)
   {
      SmallStrainUmatModel::operator=(rOther);
      return *this;
   }

   //*******************************DESTRUCTOR*******************************************
   //************************************************************************************

   HypoplasticSmallStrainUmatModel::~HypoplasticSmallStrainUmatModel()
   {
   }


} // Namespace Kratos
