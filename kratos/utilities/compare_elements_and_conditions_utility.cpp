//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Carlos Roig
//                   Vicente Mataix Ferrandiz
//


// System includes
#include <typeinfo>

// External includes

// Project includes
#include "includes/kratos_components.h"
#include "utilities/compare_elements_and_conditions_utility.h"

namespace Kratos
{
void CompareElementsAndConditionsUtility::GetRegisteredName(
    const Geometry<Node>& rGeometry,
    std::string& rName
    )
{
    KRATOS_TRY;

    for(auto const& component: KratosComponents<Geometry<Node>>::GetComponents()) {
        if (Geometry<Node>::IsSame(*(component.second), rGeometry)) {
            rName = component.first;
            return;
        }
    }

    KRATOS_ERROR << "Geometry \"" << typeid(rGeometry).name() << "\" not found!" << std::endl;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void CompareElementsAndConditionsUtility::GetRegisteredName(
    const Element& rElement,
    std::string& rName
    )
{
    KRATOS_TRY;

    for(auto const& component: KratosComponents<Element>::GetComponents()) {
        if (GeometricalObject::IsSame(*(component.second), rElement)) {
            rName = component.first;
            return;
        }
    }

    KRATOS_ERROR << "Element \"" << typeid(rElement).name() << "\" not found!" << std::endl;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void CompareElementsAndConditionsUtility::GetRegisteredName(
    const Condition& rCondition,
    std::string& rName
    )
{
    KRATOS_TRY;

    for(auto const& component: KratosComponents<Condition>::GetComponents()) {
        if (GeometricalObject::IsSame(*(component.second), rCondition)) {
            rName = component.first;
            return;
        }
    }

    KRATOS_ERROR << "Condition \"" << typeid(rCondition).name() << "\" not found!" << std::endl;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void CompareElementsAndConditionsUtility::GetRegisteredName(
    const MasterSlaveConstraint& rMasterSlaveConstraint,
    std::string& rName
    )
{
    KRATOS_TRY;

    for (const auto& r_component : KratosComponents<MasterSlaveConstraint>::GetComponents()) {
        if (typeid(*(r_component.second)) == typeid(rMasterSlaveConstraint)) {
            rName = r_component.first;
            return;
        }
    }

    KRATOS_ERROR << "MasterSlaveConstraint \"" << typeid(rMasterSlaveConstraint).name() << "\" not found!" << std::endl;

    KRATOS_CATCH("");
}


} // namespace Kratos.

