//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//

/* Project includes */
#include "custom_saturation/brooksandcorey_law.hpp"

namespace Kratos
{

int BrooksAndCoreyLaw::Check(const Properties& rMaterialProperties,
                           const GeometryType& rElementGeometry,
                           const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    int ierr = SaturationLaw::Check(rCurrentProcessInfo);
    if(ierr != 0)
        return ierr;

    return ierr;

    KRATOS_CATCH("");
}

//------------------------------------------------------------------------------------------------

void BrooksAndCoreyLaw::InitializeMaterial(const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues)
{

}

//------------------------------------------------------------------------------------------------

void BrooksAndCoreyLaw::CalculateMaterialResponse (Parameters& rValues)
{
    
}

//------------------------------------------------------------------------------------------------

void BrooksAndCoreyLaw::InitializeSaturationLawVariables (SaturationLawVariables& rVariables, Parameters& rValues)
{
    
}

}