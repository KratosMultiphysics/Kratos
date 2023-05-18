//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main author: Carlos Roig
//               Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_COMPARE_ELEMENTS_AND_CONDITIONS_UTILITY_H_INCLUDED)
#define KRATOS_COMPARE_ELEMENTS_AND_CONDITIONS_UTILITY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/condition.h"
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

    static void KRATOS_API(KRATOS_CORE) GetRegisteredName(
        const Geometry<Node>& rGeometry,
        std::string& rName);

    static void GetRegisteredName(
        const Geometry<Node>* pGeometry,
        std::string& rName)
    {
        CompareElementsAndConditionsUtility::GetRegisteredName(*pGeometry, rName);
    }

    static void KRATOS_API(KRATOS_CORE) GetRegisteredName(
        const Element& rElement,
        std::string& rName);

    static void GetRegisteredName(const Element* pElement, std::string& rName) {
        CompareElementsAndConditionsUtility::GetRegisteredName(*pElement, rName);
    }

    static void KRATOS_API(KRATOS_CORE) GetRegisteredName(
        const Condition& rCondition,
        std::string& rName);

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
