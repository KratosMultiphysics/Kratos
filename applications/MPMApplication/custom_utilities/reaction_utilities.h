//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Nicolo Crescenzio

#pragma once

// System includes

// External includes

// Project includes

namespace Kratos
{

// Auxiliary utility to compute the reaction over a line or surface
// The line or surface can be a part of the boundary of the background grid
// (conforming boundary conditions) or can be defined through material
// point conditions (non-conforming boundary conditions)
class KRATOS_API(MPM_APPLICATION) ReactionUtilities
{
public:

    /**
    * @brief Computes the sum of the REACTION variable for all nodes in a given model part
    * @param rModelPart reference to the model part where the force is to be computed
    * @return An array containing the force value
    */
    static array_1d<double, 3> CalculateGridConformingReaction(ModelPart &rModelPart);

    /**
    * @brief Computes the sum of the MPC_CONTACT_FORCE variable for all material point
    * conditions in a given model part
    * @param rModelPart reference to the model part where the force is to be computed
    * @return An array containing the force value
    */
    static array_1d<double, 3> CalculateNonConformingReaction(ModelPart &rModelPart);

}; // Class ReactionUtilities

} // namespace Kratos
