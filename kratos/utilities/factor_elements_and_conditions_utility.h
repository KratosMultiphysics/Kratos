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

#if !defined(KRATOS_FACTOR_ELEMENTS_AND_CONDITIONS_UTILITY_H_INCLUDED)
#define KRATOS_FACTOR_ELEMENTS_AND_CONDITIONS_UTILITY_H_INCLUDED

// System includes
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/model_part.h"
#include "includes/kratos_components.h"
#include "utilities/compare_elements_and_conditions_utility.h"

namespace Kratos
{
///@name Kratos Classes
///@{

class FactorElementsAndConditionsUtility
{
public:
    ///@name Type Definitions
    ///@{

    typedef ModelPart::ElementsContainerType ElementsContainerType;

    typedef ModelPart::ConditionsContainerType ConditionsContainerType;

    ///@}
    ///@name Operations
    ///@{

    static void Execute(const ElementsContainerType& rMixedElements, std::vector<ElementsContainerType>& rFactoredElements)
    {
        KRATOS_TRY;

        if (rMixedElements.size() == 0)
        {
            rFactoredElements.resize(0);
            return;
        }

        rFactoredElements.reserve(5);
        rFactoredElements.resize(1);
        unsigned pos = 0;
        rFactoredElements[pos].push_back(*rMixedElements.ptr_begin());
        for(auto it = rMixedElements.ptr_begin() + 1; it != rMixedElements.ptr_end(); ++it)
        {
            if (CompareElementAndConditionUtility::IsSame(**it, rFactoredElements[pos].front()) == false)
            {
                // Check if the new type already has a container allocated.
                bool found = false;
                for (unsigned k=0; k < rFactoredElements.size(); ++k)
                {
                    if (CompareElementAndConditionUtility::IsSame(**it, rFactoredElements[k].front()))
                    {
                        pos = k;
                        found = true;
                        break;
                    }
                }

                if (!found)
                {
                    // Allocate a new container for the new type.
                    pos = rFactoredElements.size();
                    rFactoredElements.resize(pos + 1);
                }
            }

            rFactoredElements[pos].push_back(*it);
        }

        KRATOS_CATCH("");
    }

    static void Execute(const ConditionsContainerType& rMixedConditions, std::vector<ConditionsContainerType>& rFactoredConditions)
    {
        KRATOS_TRY;

        if (rMixedConditions.size() == 0)
        {
            rFactoredConditions.resize(0);
            return;
        }

        rFactoredConditions.reserve(5);
        rFactoredConditions.resize(1);
        unsigned pos = 0;
        rFactoredConditions[pos].push_back(*rMixedConditions.ptr_begin());
        for(auto it = rMixedConditions.ptr_begin() + 1; it != rMixedConditions.ptr_end(); ++it)
        {
            if (CompareElementAndConditionUtility::IsSame(**it, rFactoredConditions[pos].front()) == false)
            {
                // Check if the new type already has a container allocated.
                bool found = false;
                for (unsigned k=0; k < rFactoredConditions.size(); ++k)
                {
                    if (CompareElementAndConditionUtility::IsSame(**it, rFactoredConditions[k].front()))
                    {
                        pos = k;
                        found = true;
                        break;
                    }
                }

                if (!found)
                {
                    // Allocate a new container for the new type.
                    pos = rFactoredConditions.size();
                    rFactoredConditions.resize(pos + 1);
                }
            }

            rFactoredConditions[pos].push_back(*it);
        }

        KRATOS_CATCH("");
    }
    
    ///@}

private:
    ///@name Private Operations
    ///@{
    ///@}
};
///@} // Kratos Classes
} // namespace Kratos.

#endif // KRATOS_FACTOR_ELEMENTS_AND_CONDITIONS_UTILITY_H_INCLUDED defined
