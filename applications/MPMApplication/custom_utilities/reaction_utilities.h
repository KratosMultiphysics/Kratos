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
  ///@addtogroup MPMApplication
  ///@{

  // Auxiliary utility to compute the reaction over a line or surface
  // The line or surface can be a part of the boundary of the background grid
  // (conforming boundary conditions) or can be defined through material
  // point conditions (non-conforming boundary conditions)
  class KRATOS_API(MPM_APPLICATION) ReactionUtilities
  {
  public:

    /// Pointer definition of ReactionUtilities
    KRATOS_CLASS_POINTER_DEFINITION(ReactionUtilities);

    ///@name Operations
    ///@{

    /**
    * @brief Computes the sum of the reaction over the nodes of the given modelpart
    * @param rModelPart reference to the model part where the force is to be computed
    * @return An array containing the force value
    */
    static array_1d<double, 3> CalculateGridConformingReaction(ModelPart &rModelPart);

    /**
    * @brief Computes the sum of the reaction over the material point conditions of the given modelpart
    * @param rModelPart reference to the model part where the force is to be computed
    * @return An array containing the force value
    */
    static array_1d<double, 3> CalculateNonConformingReaction(ModelPart &rModelPart);

}; // Class ReactionUtilities

} // namespace Kratos
