//
//   Project Name:        KratosFemToDemApplication $
//   Created by:          $Author:Alejandro Cornejo $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                Sept 2016 $
//   Revision:            $Revision:                  0.0 $
//

// System includes
#include <iostream>

// External includes
#include <cmath>

// Project includes
#include "includes/properties.h"
#include "custom_constitutive/zarate_law.hpp"

#include "fem_to_dem_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

ZarateLaw::ZarateLaw()
	: LinearElasticPlaneStrain2DLaw()
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

ZarateLaw::ZarateLaw(const ZarateLaw &rOther)
	: LinearElasticPlaneStrain2DLaw(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer ZarateLaw::Clone() const
{
	ZarateLaw::Pointer p_clone(new ZarateLaw(*this));
	return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

ZarateLaw ::~ZarateLaw()
{
}

} // namespace Kratos