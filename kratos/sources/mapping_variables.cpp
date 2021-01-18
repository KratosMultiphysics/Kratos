//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    msandre
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/mapping_variables.h"
#include "includes/kernel.h"

namespace Kratos
{
KRATOS_CREATE_VARIABLE( IndexSet::Pointer, INDEX_SET )         // An unordened set of which contains the indexes with the paired
KRATOS_CREATE_VARIABLE( IndexMap::Pointer, INDEX_MAP )         // An unordened map of which contains the indexes with the paired
KRATOS_CREATE_VARIABLE( double, TANGENT_FACTOR )               // The factor between the tangent and normal behaviour

void KratosApplication::RegisterMappingVariables()
{
    KRATOS_REGISTER_VARIABLE( INDEX_SET )                      // An unordened set of which contains the indexes with the paired
    KRATOS_REGISTER_VARIABLE( INDEX_MAP )                      // An unordened map of which contains the indexes with the paired
    KRATOS_REGISTER_VARIABLE( TANGENT_FACTOR )                 // The factor between the tangent and normal behaviour
}

}  // namespace Kratos.
