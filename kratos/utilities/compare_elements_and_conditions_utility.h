//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main author:    
//

#if !defined(KRATOS_COMPARE_ELEMENTS_AND_CONDITIONS_UTILITY_H_INCLUDED)
#define KRATOS_COMPARE_ELEMENTS_AND_CONDITIONS_UTILITY_H_INCLUDED

// System includes
#include <typeinfo>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/kratos_components.h"
#include "includes/geometrical_object.h"

namespace Kratos
{
///@name Kratos Classes
///@{
// A utility for comparison of dynamic types of elements and conditions.
class CompareElementsAndConditionsUtility
{
public:
    ///@name Type Definitions
    ///@{
    ///@}
    ///@name Operations
    ///@{

    static void GetRegisteredName(const Element& rElement, std::string& rName) {
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

    static void GetRegisteredName(const Element* pElement, std::string& rName) {
        CompareElementsAndConditionsUtility::GetRegisteredName(*pElement, rName);
    }

    static void GetRegisteredName(const Condition& rCondition, std::string& rName) {
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

    static void GetRegisteredName(const Condition* pCondition, std::string& rName) {
        CompareElementsAndConditionsUtility::GetRegisteredName(*pCondition, rName);
    }

    ///@}

private:
    ///@name Private Operations
    ///@{
    ///@}
};
///@} // Kratos Classes
} // namespace Kratos.

#endif // KRATOS_COMPARE_ELEMENTS_AND_CONDITIONS_UTILITY_H_INCLUDED defined
