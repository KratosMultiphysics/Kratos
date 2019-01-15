//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//

// System includes

// External includes

// Project includes
#include "structural_mechanics_element_utilities.h"
#include "includes/variables.h"

namespace Kratos {
namespace StructuralMechanicsElementUtilities {

bool ComputeLumpedMassMatrix(
    const Properties& rProperites,
    const ProcessInfo& rCurrentProcessInfo)
{
    // giving the globally defined setting (through ProcessInfo) priority
    // over the locally defined one (through Properties)
    // this is needed for the explicit solver, which requires the lumped
    // mass matrix and specifies it's computation through the ProcessInfo
    if (rCurrentProcessInfo.Has(COMPUTE_LUMPED_MASS_MATRIX)) {
        return rCurrentProcessInfo[COMPUTE_LUMPED_MASS_MATRIX];
    }
    else if (rProperites.Has(COMPUTE_LUMPED_MASS_MATRIX)) {
        return rProperites[COMPUTE_LUMPED_MASS_MATRIX];
    }

    // the default for all elements in StructuralMechanics is
    // to use the consistent-mass-matrix, hence returning false here
    return false;
}

} // namespace StructuralMechanicsElementUtilities.
}  // namespace Kratos.


