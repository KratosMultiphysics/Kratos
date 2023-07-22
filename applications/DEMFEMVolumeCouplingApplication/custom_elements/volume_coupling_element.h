// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"
// #include "includes/element.h"
// #include "utilities/integration_utilities.h"
// #include "structural_mechanics_application_variables.h"
// #include "utilities/geometrical_sensitivity_utility.h"
// #include "custom_utilities/structural_mechanics_element_utilities.h"
#include "custom_elements/small_displacement.h"

namespace Kratos
{

class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) VolumeCouplingElement
    : public SmallDisplacement
{


public:

   
    virtual double GetIntegrationWeight(
        const GeometryType::IntegrationPointsArrayType& rThisIntegrationPoints,
        const IndexType PointNumber,
        const double detJ
        ) const override;


};
} // namespace Kratos.
