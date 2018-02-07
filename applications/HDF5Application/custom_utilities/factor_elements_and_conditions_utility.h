//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: HDF5Application/license.txt
//
//  Main author:    Michael Andre, https://github.com/msandre
//

#if !defined(KRATOS_FACTOR_ELEMENTS_AND_CONDITIONS_UTILITY_H_INCLUDED)
#define KRATOS_FACTOR_ELEMENTS_AND_CONDITIONS_UTILITY_H_INCLUDED

// System includes
#include <vector>
#include <string>
#include <utility>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/model_part.h"

namespace Kratos
{
///@name Kratos Classes
///@{

class FactorElementsUtility
{
public:
    ///@name Type Definitions
    ///@{

    typedef ModelPart::ElementsContainerType ElementsContainerType;

    typedef std::vector<std::pair<std::string, ElementsContainerType>>::iterator iterator;

    typedef std::vector<std::pair<std::string, ElementsContainerType>>::const_iterator const_iterator;

    ///@}
    ///@name Life Cycle
    ///@{

    explicit FactorElementsUtility(const ElementsContainerType& rElements);

    ///@}
    ///@name Operations
    ///@{

    iterator begin()
    {
        return iterator(mFactoredElements.begin());
    }

    const_iterator begin() const
    {
        return const_iterator(mFactoredElements.begin());
    }

    iterator end()
    {
        return iterator(mFactoredElements.end());
    }

    const_iterator end() const
    {
        return const_iterator(mFactoredElements.end());
    }
    
    ///@}

private:
    ///@name Member Variables
    ///@{

    std::vector<std::pair<std::string, ElementsContainerType>> mFactoredElements;

    ///@}
    ///@name Private Operations
    ///@{
    ///@}
};

class FactorConditionsUtility
{
public:
    ///@name Type Definitions
    ///@{

    typedef ModelPart::ConditionsContainerType ConditionsContainerType;

    typedef std::vector<std::pair<std::string, ConditionsContainerType>>::iterator iterator;

    typedef std::vector<std::pair<std::string, ConditionsContainerType>>::const_iterator const_iterator;

    ///@}
    ///@name Life Cycle
    ///@{

    explicit FactorConditionsUtility(const ConditionsContainerType& rConditions);

    ///@}
    ///@name Operations
    ///@{

    iterator begin()
    {
        return iterator(mFactoredConditions.begin());
    }

    const_iterator begin() const
    {
        return const_iterator(mFactoredConditions.begin());
    }

    iterator end()
    {
        return iterator(mFactoredConditions.end());
    }

    const_iterator end() const
    {
        return const_iterator(mFactoredConditions.end());
    }
    
    ///@}

private:
    ///@name Member Variables
    ///@{

    std::vector<std::pair<std::string, ConditionsContainerType>> mFactoredConditions;

    ///@}
    ///@name Private Operations
    ///@{
    ///@}
};

///@} // Kratos Classes
} // namespace Kratos.

#endif // KRATOS_FACTOR_ELEMENTS_AND_CONDITIONS_UTILITY_H_INCLUDED defined
