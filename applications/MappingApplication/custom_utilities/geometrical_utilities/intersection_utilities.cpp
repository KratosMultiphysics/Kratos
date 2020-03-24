//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

// System includes

// External includes

// Project includes
#include "intersection_utilities.h"

namespace Kratos
{

void IntersectionUtilities::FindIntersection1DGeometries2D(
    ModelPart& rModelPartDomainA,
    ModelPart& rModelPartDomainB,
    ModelPart& rModelPartResult,
    double Tolerance)
{
    KRATOS_ERROR_IF(rModelPartDomainA.ConditionsBegin()->GetGeometry().Dimension() != 1)
        << "Can compare only line segments with other line segments." << std::endl;

    for (auto condition_a_itr = rModelPartDomainA.ConditionsBegin();
        condition_a_itr != rModelPartDomainA.ConditionsEnd();
        ++condition_a_itr)
    {
        for (auto condition_b_itr = rModelPartDomainB.ConditionsBegin();
            condition_b_itr != rModelPartDomainB.ConditionsEnd();
            ++condition_b_itr)
        {
            // ERROR HasIntersection does not properly show parallel curves, which are the significant subset
            if (condition_a_itr->GetGeometry().HasIntersection(condition_b_itr->GetGeometry()))
            {
                rModelPartResult.AddGeometry(Kratos::make_shared<CouplingGeometry<NodeType>>(
                    condition_a_itr->pGetGeometry(), condition_b_itr->pGetGeometry()));
            }
        }
    }
}

} // namespace Kratos.
