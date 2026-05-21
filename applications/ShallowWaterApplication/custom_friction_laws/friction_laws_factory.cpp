//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

// System includes


// External includes


// Project includes
#include "friction_laws_factory.h"
#include "nodal_manning_law.h"
#include "manning_law.h"
#include "chezy_law.h"
#include "wind_water_friction.h"


namespace Kratos
{

FrictionLaw::Pointer FrictionLawsFactory::CreateBottomFrictionLaw(
    const GeometryType& rGeometry,
    const Properties& rProperty,
    const ProcessInfo& rProcessInfo)
{
    if (rProperty.Has(MANNING)) {
        return Kratos::make_shared<ManningLaw>(rGeometry, rProperty, rProcessInfo);
    }
    else if (rProperty.Has(CHEZY)) {
        return Kratos::make_shared<ChezyLaw>(rGeometry, rProperty, rProcessInfo);
    }
    else if (rGeometry[0].SolutionStepsDataHas(MANNING)) {
        return Kratos::make_shared<NodalManningLaw>(rGeometry, rProperty, rProcessInfo);
    }
    else {
        return Kratos::make_shared<FrictionLaw>();
    }
}

FrictionLaw::Pointer FrictionLawsFactory::CreateSurfaceFrictionLaw(
    const GeometryType& rGeometry,
    const Properties& rProperty,
    const ProcessInfo& rProcessInfo)
{
    if (rProcessInfo.Has(DENSITY_AIR) && rGeometry[0].SolutionStepsDataHas(WIND)) {
        return Kratos::make_shared<WindWaterFriction>(rGeometry, rProperty, rProcessInfo);
    }
    else {
        return Kratos::make_shared<FrictionLaw>();
    }
}

}  // namespace Kratos
