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

// External includes

// Project includes
#include "includes/element.h"
#include "includes/condition.h"
#include "utilities/indexed_object.h"
#include "containers/pointer_vector_set.h"

namespace Kratos
{

typedef PointerVectorSet<Element, IndexedObject> ElementsContainerType;

typedef PointerVectorSet<Condition, IndexedObject> ConditionsContainerType;

/// Factor a collection of elements into uniform containers.
/**
 * Each container consists of all elements of a single type. In MPI,
 * the sequence of element containers is global (i.e., the ith container on
 * each process corresponds to the same element type), but the contents of
 * each container on a process are local. If a process has no local elements 
 * corresponding to the ith container, it is empty.
 */
std::vector<ElementsContainerType> FactorElements(ElementsContainerType const& rElements);

std::vector<ConditionsContainerType> FactorConditions(ConditionsContainerType const& rConditions);

} // namespace Kratos.

#endif // KRATOS_FACTOR_ELEMENTS_AND_CONDITIONS_UTILITY_H_INCLUDED defined
