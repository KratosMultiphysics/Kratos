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
    const std::string& rEntityName,
    const bool Initialize,
    const bool Elements
    )
{
    Properties::Pointer p_prop = rModelPart.HasProperties(0) ? rModelPart.pGetProperties(0) : rModelPart.CreateNewProperties(0);

    // First we create the nodes
    rModelPart.CreateNewNode(1, 0.0 , 0.0 , 0.0);
    rModelPart.CreateNewNode(2, 1.0 , 0.0 , 0.0);
    rModelPart.CreateNewNode(3, 1.0 , 1.0 , 0.0);
    rModelPart.CreateNewNode(4, 0.0 , 1.0 , 0.0);
    rModelPart.CreateNewNode(5, 2.0 , 0.0 , 0.0);
    rModelPart.CreateNewNode(6, 2.0 , 1.0 , 0.0);

    if (Elements) {
        rModelPart.CreateNewElement(rEntityName, 1, {{1,2,3}}, p_prop);
        rModelPart.CreateNewElement(rEntityName, 2, {{1,3,4}}, p_prop);
        rModelPart.CreateNewElement(rEntityName, 3, {{2,5,3}}, p_prop);
        rModelPart.CreateNewElement(rEntityName, 4, {{5,6,3}}, p_prop);

        // Initialize Elements
        if (Initialize) {
            const auto& r_process_info = rModelPart.GetProcessInfo();
            for (auto& r_elem : rModelPart.Elements())
                r_elem.Initialize(r_process_info);
        }
    } else {
        rModelPart.CreateNewCondition(rEntityName, 1, {{1,2,3}}, p_prop);
        rModelPart.CreateNewCondition(rEntityName, 2, {{1,3,4}}, p_prop);
        rModelPart.CreateNewCondition(rEntityName, 3, {{2,5,3}}, p_prop);
        rModelPart.CreateNewCondition(rEntityName, 4, {{5,6,3}}, p_prop);

        // Initialize Elements
        if (Initialize) {
            const auto& r_process_info = rModelPart.GetProcessInfo();
            for (auto& r_cond : rModelPart.Conditions())
                r_cond.Initialize(r_process_info);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void Create2DQuadrilateralsGeometry(
    ModelPart& rModelPart, 
    const std::string& rEntityName,
    const bool Initialize,
    const bool Elements
    )
{
    Properties::Pointer p_prop = rModelPart.HasProperties(0) ? rModelPart.pGetProperties(0) : rModelPart.CreateNewProperties(0);

    // First we create the nodes
    rModelPart.CreateNewNode(1, 0.0 , 0.0 , 0.0);
    rModelPart.CreateNewNode(2, 1.0 , 0.0 , 0.0);
    rModelPart.CreateNewNode(3, 1.0 , 1.0 , 0.0);
    rModelPart.CreateNewNode(4, 0.0 , 1.0 , 0.0);
    rModelPart.CreateNewNode(5, 2.0 , 0.0 , 0.0);
    rModelPart.CreateNewNode(6, 2.0 , 1.0 , 0.0);

    if (Elements) {
        rModelPart.CreateNewElement(rEntityName, 1, {{1,2,3,4}}, p_prop);
        rModelPart.CreateNewElement(rEntityName, 2, {{2,5,6,3}}, p_prop);

        // Initialize Elements
        if (Initialize) {
            const auto& r_process_info = rModelPart.GetProcessInfo();
            for (auto& r_elem : rModelPart.Elements())
                r_elem.Initialize(r_process_info);
        }
    } else {
        rModelPart.CreateNewCondition(rEntityName, 1, {{1,2,3,4}}, p_prop);
        rModelPart.CreateNewCondition(rEntityName, 2, {{2,5,6,3}}, p_prop);

        // Initialize Elements
        if (Initialize) {
            const auto& r_process_info = rModelPart.GetProcessInfo();
            for (auto& r_cond : rModelPart.Conditions())
                r_cond.Initialize(r_process_info);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void Create3DGeometry(
    ModelPart& rModelPart,
    const std::string& rElementName,
    const bool Initialize
    )
{
    Properties::Pointer p_prop = rModelPart.HasProperties(0) ? rModelPart.pGetProperties(0) : rModelPart.CreateNewProperties(0);

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

    rModelPart.CreateNewElement(rElementName, 1, {{12,10,8,9}}, p_prop);
    rModelPart.CreateNewElement(rElementName, 2, {{4,6,9,7}}, p_prop);
    rModelPart.CreateNewElement(rElementName, 3, {{11,7,9,8}}, p_prop);
    rModelPart.CreateNewElement(rElementName, 4, {{5,3,8,6}}, p_prop);
    rModelPart.CreateNewElement(rElementName, 5, {{4,6,7,3}}, p_prop);
    rModelPart.CreateNewElement(rElementName, 6, {{2,3,5,6}}, p_prop);
    rModelPart.CreateNewElement(rElementName, 7, {{10,9,6,8}}, p_prop);
    rModelPart.CreateNewElement(rElementName, 8, {{7,8,3,6}}, p_prop);
    rModelPart.CreateNewElement(rElementName, 9, {{7,8,6,9}}, p_prop);
    rModelPart.CreateNewElement(rElementName, 10, {{4,1,6,3}}, p_prop);
    rModelPart.CreateNewElement(rElementName, 11, {{9,12,11,8}}, p_prop);
    rModelPart.CreateNewElement(rElementName, 12, {{3,2,1,6}}, p_prop);

    // Initialize Elements
    if (Initialize) {
        const auto& r_process_info = rModelPart.GetProcessInfo();
        for (auto& r_elem : rModelPart.Elements())
            r_elem.Initialize(r_process_info);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void Create3DHexahedraGeometry(
    ModelPart& rModelPart,
    const std::string& rElementName,
    const bool Initialize
    )
{
    Properties::Pointer p_prop = rModelPart.HasProperties(0) ? rModelPart.pGetProperties(0) : rModelPart.CreateNewProperties(0);

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

    rModelPart.CreateNewElement(rElementName, 1, {{5,8,6,2,3,7,4,1}}, p_prop);
    rModelPart.CreateNewElement(rElementName, 2, {{8,12,10,6,7,11,9,4}}, p_prop);

    // Initialize Elements
    if (Initialize) {
        const auto& r_process_info = rModelPart.GetProcessInfo();
        for (auto& r_elem : rModelPart.Elements())
            r_elem.Initialize(r_process_info);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void Create3DQuadraticGeometry(
    ModelPart& rModelPart, 
    const std::string& rElementName,
    const bool Initialize
    )
{
    Properties::Pointer p_prop = rModelPart.HasProperties(0) ? rModelPart.pGetProperties(0) : rModelPart.CreateNewProperties(0);

    // First we create the nodes
    rModelPart.CreateNewNode(1, 1.0000000000, 0.0000000000, 1.0000000000);
    rModelPart.CreateNewNode(2, 0.5000000000, 0.0000000000, 1.0000000000);
    rModelPart.CreateNewNode(3, 1.0000000000, 0.5000000000, 1.0000000000);
    rModelPart.CreateNewNode(4, 1.0000000000, 0.0000000000, 0.5000000000);
    rModelPart.CreateNewNode(5, 0.5000000000, 0.0000000000, 0.5000000000);
    rModelPart.CreateNewNode(6, 1.0000000000, 0.5000000000, 0.5000000000);
    rModelPart.CreateNewNode(7, 0.5000000000, 0.5000000000, 1.0000000000);
    rModelPart.CreateNewNode(8, 0.5000000000, 0.5000000000, 0.5000000000);
    rModelPart.CreateNewNode(9, 0.2019246055, 0.3959160307, 0.6930668948);
    rModelPart.CreateNewNode(10, 0.7019246055, 0.8959160307, 0.6930668948);
    rModelPart.CreateNewNode(11, 1.0000000000, 0.0000000000, 0.0000000000);
    rModelPart.CreateNewNode(12, 0.0000000000, 0.0000000000, 1.0000000000);
    rModelPart.CreateNewNode(13, 1.0000000000, 1.0000000000, 1.0000000000);
    rModelPart.CreateNewNode(14, 0.5000000000, 0.0000000000, 0.0000000000);
    rModelPart.CreateNewNode(15, 1.0000000000, 0.5000000000, 0.0000000000);
    rModelPart.CreateNewNode(16, 0.5000000000, 1.0000000000, 1.0000000000);
    rModelPart.CreateNewNode(17, 0.0000000000, 0.5000000000, 1.0000000000);
    rModelPart.CreateNewNode(18, 0.0000000000, 0.0000000000, 0.5000000000);
    rModelPart.CreateNewNode(19, 1.0000000000, 1.0000000000, 0.5000000000);
    rModelPart.CreateNewNode(20, 0.4038492111, 0.7918320615, 0.3861337896);
    rModelPart.CreateNewNode(21, 0.2019246055, 0.3959160307, 0.1930668948);
    rModelPart.CreateNewNode(22, 0.5000000000, 0.5000000000, 0.0000000000);
    rModelPart.CreateNewNode(23, 0.5000000000, 1.0000000000, 0.5000000000);
    rModelPart.CreateNewNode(24, 0.0000000000, 0.5000000000, 0.5000000000);
    rModelPart.CreateNewNode(25, 0.2019246055, 0.8959160307, 0.6930668948);
    rModelPart.CreateNewNode(26, 0.7019246055, 0.8959160307, 0.1930668948);
    rModelPart.CreateNewNode(27, 0.0000000000, 0.0000000000, 0.0000000000);
    rModelPart.CreateNewNode(28, 1.0000000000, 1.0000000000, 0.0000000000);
    rModelPart.CreateNewNode(29, 0.0000000000, 1.0000000000, 1.0000000000);
    rModelPart.CreateNewNode(30, 0.2019246055, 0.8959160307, 0.1930668948);
    rModelPart.CreateNewNode(31, 0.5000000000, 1.0000000000, 0.0000000000);
    rModelPart.CreateNewNode(32, 0.0000000000, 0.5000000000, 0.0000000000);
    rModelPart.CreateNewNode(33, 0.0000000000, 1.0000000000, 0.5000000000);
    rModelPart.CreateNewNode(34, 0.0000000000, 1.0000000000, 0.0000000000);

    // Now we create the elements
    rModelPart.CreateNewElement(rElementName, 1,  {{27, 11, 28, 13, 14, 15, 22,  8,  6, 19}}, p_prop);
    rModelPart.CreateNewElement(rElementName, 2,  {{34, 12, 27, 20, 24, 18, 32, 30,  9, 21}}, p_prop);
    rModelPart.CreateNewElement(rElementName, 3,  {{27, 34, 20, 28, 32, 30, 21, 22, 31, 26}}, p_prop);
    rModelPart.CreateNewElement(rElementName, 4,  {{34, 20, 28, 29, 30, 26, 31, 33, 25, 23}}, p_prop);
    rModelPart.CreateNewElement(rElementName, 5,  {{34, 20, 29, 12, 30, 25, 33, 24,  9, 17}}, p_prop);
    rModelPart.CreateNewElement(rElementName, 6,  {{20, 27, 28, 13, 21, 22, 26, 10,  8, 19}}, p_prop);
    rModelPart.CreateNewElement(rElementName, 7,  {{28, 20, 13, 29, 26, 10, 19, 23, 25, 16}}, p_prop);
    rModelPart.CreateNewElement(rElementName, 8,  {{20, 13, 29, 12, 10, 16, 25,  9,  7, 17}}, p_prop);
    rModelPart.CreateNewElement(rElementName, 9,  {{27,  1, 11, 13,  5,  4, 14,  8,  3,  6}}, p_prop);
    rModelPart.CreateNewElement(rElementName, 10, {{12, 13,  1, 27,  7,  3,  2, 18,  8,  5}}, p_prop);
    rModelPart.CreateNewElement(rElementName, 11, {{12, 27, 20, 13, 18, 21,  9,  7,  8, 10}}, p_prop);

    // Initialize Elements
    if (Initialize) {
        const auto& r_process_info = rModelPart.GetProcessInfo();
        for (auto& r_elem : rModelPart.Elements())
            r_elem.Initialize(r_process_info);
    }
}

} // namespace ConstraintUtilities
} // namespace Kratos
