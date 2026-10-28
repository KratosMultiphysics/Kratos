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

#include "custom_utilities/fluid_patch_gap_element_check.h"

#include "custom_elements/symbolic_stokes.h"
#include "custom_elements/symbolic_stokes_ale.h"
#include "custom_elements/data_containers/stokes/stokes_data.h"
#include "custom_elements/qs_vms.h"
#include "custom_elements/data_containers/time_integrated_qs_vms/time_integrated_qs_vms_data.h"

namespace Kratos
{
void FluidPatchGapElementCheck::Check(
    const Element& rElement, const bool UseQSVMS, const bool IsMovingInnerPatch)
{
    using StokesElementType = SymbolicStokes<SymbolicStokesData<2, 4>>;
    using QSVMSQ1ElementType = QSVMS<TimeIntegratedQSVMSData<2, 4>>;
    if (UseQSVMS) {
        KRATOS_ERROR_IF(dynamic_cast<const QSVMSQ1ElementType*>(&rElement) == nullptr)
            << "QSVMS gap coupling requires TimeIntegratedQSVMS2D4N on both patches." << std::endl;
    } else {
        KRATOS_ERROR_IF(dynamic_cast<const StokesElementType*>(&rElement) == nullptr)
            << "Gap-SBM currently requires SymbolicStokes2D4N bulk elements." << std::endl;
        KRATOS_ERROR_IF(IsMovingInnerPatch && dynamic_cast<const SymbolicStokesALE2D4N*>(&rElement) == nullptr)
            << "A moving patch requires SymbolicStokesALE2D4N." << std::endl;
    }
}
} // namespace Kratos
