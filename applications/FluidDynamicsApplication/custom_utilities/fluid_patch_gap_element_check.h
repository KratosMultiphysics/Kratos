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
/** @brief Check the element types supported by fluid patch gap coupling. */
class KRATOS_API(FLUID_DYNAMICS_APPLICATION) FluidPatchGapElementCheck
{
public:
    static void Check(const Element& rElement, bool UseQSVMS, bool IsMovingInnerPatch);
};
} // namespace Kratos
