//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Marc Nunez, Inigo Lopez based on R.Rossi and V.Mataix work
//
//
//

/* System includes */

/* External includes */

/* Project includes */
#include "find_cut_skin_entities_process.h"
#include "includes/cfd_variables.h"

namespace Kratos
{

// Default constructor
FindCutSkinEntitiesProcess::FindCutSkinEntitiesProcess(
    ModelPart& rModelPart,
    ModelPart& rSectionModelPart,
    const array_1d<double,3>& rVersor,
    const array_1d<double,3>& rOrigin)
    :mrModelPart(rModelPart),
    mrSectionModelPart(rSectionModelPart),
    mrVersor(rVersor),
    mrOrigin(rOrigin)
{
    KRATOS_TRY

    for (auto& r_node : rModelPart.Nodes()) {
        auto direction = r_node.Coordinates() - mrOrigin;
        auto distance = inner_prod(direction, rVersor);
        if (std::abs(distance) < 1e-9) {
            r_node.SetValue(DISTANCE, 1e-9);
        } else {
            r_node.SetValue(DISTANCE, distance);
        }
    }
    std::size_t node_index = 0;

    for (auto& r_cond : rModelPart.Conditions()) {
        auto& r_geometry = r_cond.GetGeometry();
        std::size_t n_pos = 0;
        std::size_t n_neg = 0;

        for (auto& r_node : r_geometry) {
            auto distance = r_node.GetValue(DISTANCE);
            if (distance > 0.0) {
                n_pos++;
            }
            else {
                n_neg++;
            }

        }

        if (n_neg > 0 && n_pos > 0) {
            node_index++;
            auto p_node = rSectionModelPart.CreateNewNode(node_index, r_geometry.Center().X(), r_geometry.Center().Y(), r_geometry.Center().Z());
            p_node->SetValue(PRESSURE_COEFFICIENT, r_cond.GetValue(PRESSURE_COEFFICIENT));
        }
    }



    KRATOS_CATCH("")
}

void FindCutSkinEntitiesProcess::Execute()
{
    KRATOS_TRY;

    KRATOS_CATCH("")
}

} /* namespace Kratos.*/
