//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:             BSD License
//                               Kratos default license:
//kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
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

ConstitutiveLaw::Pointer ZarateLaw::Clone() 
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