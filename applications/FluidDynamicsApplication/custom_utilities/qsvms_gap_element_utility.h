//     |  /           |
//     ' /   __| _` | __|  _ \   __|
//     . \  |   (   | |   (   |\__ `
//    _|\_\_|  \__,_|\__|\___/ ____/
//                    Multi-Physics
//
//   Fluid Dynamics Application
//
//   License:         BSD License
//                    Kratos default license: kratos/license.txt
//
//   Main author:     Nicolò Antonelli
//

#pragma once

#include "includes/element.h"

namespace Kratos
{
/** @brief Create a QSVMS gap integration element using the parent Q1 nodes.
 * @details The quadrature geometry must already contain the parent nodes and
 * its custom integration rule. Element data needed by QSVMS is retained.
 * The caller owns insertion and initialization in the destination ModelPart.
 */
class KRATOS_API(FLUID_DYNAMICS_APPLICATION) QSVMSGapElementUtility
{
public:
    static Element::Pointer Create(
        Element::IndexType NewId,
        Element& rParent,
        Geometry<Node>::Pointer pQuadratureGeometry);
};
} // namespace Kratos
