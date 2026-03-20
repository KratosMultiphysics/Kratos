// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Tobias Siemer
//

// System includes
#include <array>
#include <vector>

// External includes

// Project includes
#include "custom_utilities/cross_section_properties_utility.h"
#include "structural_mechanics_application_variables.h"
#include "utilities/beam_section_properties_utility.h"

namespace Kratos
{
namespace
{

const std::array<const Variable<double>*, 10> msDimensionVariables = {{
    &DIM1, &DIM2, &DIM3, &DIM4, &DIM5, &DIM6, &DIM7, &DIM8, &DIM9, &DIM10
}};

std::vector<double> CollectDimensions(const Properties& rProperties)
{
    std::vector<double> dimensions;
    dimensions.reserve(msDimensionVariables.size());

    for (const auto p_variable : msDimensionVariables) {
        if (!rProperties.Has(*p_variable)) {
            break;
        }
        dimensions.push_back(rProperties.GetValue(*p_variable));
    }

    return dimensions;
}

bool HasCompleteBeamSectionProperties(const Properties& rProperties)
{
    return rProperties.Has(CROSS_AREA)
        && rProperties.Has(I22)
        && rProperties.Has(I33)
        && rProperties.Has(TORSIONAL_INERTIA)
        && rProperties.Has(AREA_EFFECTIVE_Y)
        && rProperties.Has(AREA_EFFECTIVE_Z);
}

} // namespace

void CrossSectionPropertiesUtility::CalculateBeamProperties(ModelPart& rModelPart)
{
    for (auto& r_properties : rModelPart.rProperties()) {
        if (!r_properties.Has(CROSS_SECTION_TYPE) || HasCompleteBeamSectionProperties(r_properties)) {
            continue;
        }

        const auto section_properties = BeamSectionPropertiesUtility::CalculateProperties(
            r_properties.GetValue(CROSS_SECTION_TYPE),
            CollectDimensions(r_properties));

        r_properties.SetValue(CROSS_AREA, section_properties.Area);
        r_properties.SetValue(I22, section_properties.I22);
        r_properties.SetValue(I33, section_properties.I33);
        r_properties.SetValue(TORSIONAL_INERTIA, section_properties.TorsionalInertia);
        r_properties.SetValue(AREA_EFFECTIVE_Y, section_properties.ShearFactorY * section_properties.Area);
        r_properties.SetValue(AREA_EFFECTIVE_Z, section_properties.ShearFactorZ * section_properties.Area);
    }
}

} // namespace Kratos
