// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Vicente Mataix Ferrandiz
//                   Alejandro Cornejo Velazquez
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "utilities/integration_utilities.h"
#include "structural_mechanics_application_variables.h"
#include "utilities/geometrical_sensitivity_utility.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    /// The definition of the index type
    typedef std::size_t IndexType;

    /// The definition of the sizetype
    typedef std::size_t SizeType;

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class BaseSolidElement
 * @ingroup StructuralMechanicsApplication
 * @brief This is base class used to define the solid elements
 * @details The elements derived from this class are the small displacement element, the total lagrangian element and the updated lagrangian element
 * @author Riccardo Rossi
 * @author Vicente Mataix Ferrandiz
 */
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
