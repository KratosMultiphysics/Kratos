//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . <Kratos> utility for beam section properties
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Tobias Siemer
//

#pragma once

// System includes
#include <string>
#include <vector>

// External includes

// Project includes
#include "includes/define.h"

namespace Kratos
{

class KRATOS_API(KRATOS_CORE) BeamSectionPropertiesUtility
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(BeamSectionPropertiesUtility);

    struct SectionProperties
    {
        double Area = 0.0;
        double I22 = 0.0;
        double I33 = 0.0;
        double TorsionalInertia = 0.0;
        double ShearFactorY = 0.0;
        double ShearFactorZ = 0.0;
    };

    BeamSectionPropertiesUtility() = delete;

    static SectionProperties CalculateProperties(
        const std::string& rSectionType,
        const std::vector<double>& rDimensions);
};

} // namespace Kratos
