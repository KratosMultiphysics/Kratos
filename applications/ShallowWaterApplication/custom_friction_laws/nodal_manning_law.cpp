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
#include "nodal_manning_law.h"
#include "shallow_water_application_variables.h"
#include "custom_utilities/shallow_water_utilities.h"


namespace Kratos
{

NodalManningLaw::NodalManningLaw(
    const GeometryType& rGeometry,
    const Properties& rProperty,
    const ProcessInfo& rProcessInfo)
{
    this->Initialize(rGeometry, rProperty, rProcessInfo);
}

void NodalManningLaw::Initialize(
    const GeometryType& rGeometry,
    const Properties& rProperty,
    const ProcessInfo& rProcessInfo)
{
    double manning = 0.0;
    for (auto& r_node : rGeometry) {
        manning += r_node.FastGetSolutionStepValue(MANNING);
    }
    manning /= rGeometry.size();
    mManning2 = std::pow(manning, 2);

    mEpsilon = rGeometry.Length() * rProcessInfo[RELATIVE_DRY_HEIGHT];
}

}  // namespace Kratos
