//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "utilities/cpp_tests_utilities.h"

namespace Kratos
{
namespace CppTestsUtilities
{
void Create2DGeometry(
    ModelPart& rModelPart, 
    const std::string& rElementName
    )
{
    Properties::Pointer p_elem_prop = rModelPart.CreateNewProperties(0);

    // First we create the nodes
    rModelPart.CreateNewNode(1, 0.0 , 0.0 , 0.0);
    rModelPart.CreateNewNode(2, 1.0 , 0.0 , 0.0);
    rModelPart.CreateNewNode(3, 1.0 , 1.0 , 0.0);
    rModelPart.CreateNewNode(4, 0.0 , 1.0 , 0.0);
    rModelPart.CreateNewNode(5, 2.0 , 0.0 , 0.0);
    rModelPart.CreateNewNode(6, 2.0 , 1.0 , 0.0);

    rModelPart.CreateNewElement(rElementName, 1, {{1,2,3}}, p_elem_prop);
    rModelPart.CreateNewElement(rElementName, 2, {{1,3,4}}, p_elem_prop);
    rModelPart.CreateNewElement(rElementName, 3, {{2,5,3}}, p_elem_prop);
    rModelPart.CreateNewElement(rElementName, 4, {{5,6,3}}, p_elem_prop);
}

/***********************************************************************************/
/***********************************************************************************/

void Create3DGeometry(
    ModelPart& rModelPart, 
    const std::string& rElementName
    )
{
    Properties::Pointer p_elem_prop = rModelPart.CreateNewProperties(0);

    // First we create the nodes
    rModelPart.CreateNewNode(1 , 0.0 , 1.0 , 1.0);
    rModelPart.CreateNewNode(2 , 0.0 , 1.0 , 0.0);
    rModelPart.CreateNewNode(3 , 0.0 , 0.0 , 1.0);
    rModelPart.CreateNewNode(4 , 1.0 , 1.0 , 1.0);
    rModelPart.CreateNewNode(5 , 0.0 , 0.0 , 0.0);
    rModelPart.CreateNewNode(6 , 1.0 , 1.0 , 0.0);

    rModelPart.CreateNewNode(7 , 1.0 , 0.0 , 1.0);
    rModelPart.CreateNewNode(8 , 1.0 , 0.0 , 0.0);
    rModelPart.CreateNewNode(9 , 2.0 , 1.0 , 1.0);
    rModelPart.CreateNewNode(10 , 2.0 , 1.0 , 0.0);
    rModelPart.CreateNewNode(11 , 2.0 , 0.0 , 1.0);
    rModelPart.CreateNewNode(12 , 2.0 , 0.0 , 0.0);

    rModelPart.CreateNewElement(rElementName, 1, {{12,10,8,9}}, p_elem_prop);
    rModelPart.CreateNewElement(rElementName, 2, {{4,6,9,7}}, p_elem_prop);
    rModelPart.CreateNewElement(rElementName, 3, {{11,7,9,8}}, p_elem_prop);
    rModelPart.CreateNewElement(rElementName, 4, {{5,3,8,6}}, p_elem_prop);
    rModelPart.CreateNewElement(rElementName, 5, {{4,6,7,3}}, p_elem_prop);
    rModelPart.CreateNewElement(rElementName, 6, {{2,3,5,6}}, p_elem_prop);
    rModelPart.CreateNewElement(rElementName, 7, {{10,9,6,8}}, p_elem_prop);
    rModelPart.CreateNewElement(rElementName, 8, {{7,8,3,6}}, p_elem_prop);
    rModelPart.CreateNewElement(rElementName, 9, {{7,8,6,9}}, p_elem_prop);
    rModelPart.CreateNewElement(rElementName, 10, {{4,1,6,3}}, p_elem_prop);
    rModelPart.CreateNewElement(rElementName, 11, {{9,12,11,8}}, p_elem_prop);
    rModelPart.CreateNewElement(rElementName, 12, {{3,2,1,6}}, p_elem_prop);
}

} // namespace ConstraintUtilities
} // namespace Kratos
