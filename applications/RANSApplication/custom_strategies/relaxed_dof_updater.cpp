//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// Project includes
#include "spaces/ublas_space.h"
#include "utilities/parallel_utilities.h"

// Include base h
#include "relaxed_dof_updater.h"

namespace Kratos
{
template <>
void RelaxedDofUpdater<UblasSpace<double, CompressedMatrix, Vector>>::Initialize(
    const DofsArrayType& rDofSet,
    const SystemVectorType& rDx)
{
}

template <>
void RelaxedDofUpdater<UblasSpace<double, CompressedMatrix, Vector>>::Clear()
{
}

template <>
void RelaxedDofUpdater<UblasSpace<double, CompressedMatrix, Vector>>::UpdateDofs(
    DofsArrayType& rDofSet,
    const SystemVectorType& rDx)
{
    block_for_each(rDofSet, [&](DofType& rDof) {
        if (rDof.IsFree()) {
            rDof.GetSolutionStepValue() += rDx[rDof.EquationId()] * mRelaxationFactor;
        }
    });
}

template <>
std::string RelaxedDofUpdater<UblasSpace<double, CompressedMatrix, Vector>>::Info() const
{
    std::stringstream buffer;
    buffer << "RelaxedDofUpdater - UblasSpace";
    return buffer.str();
}

///@}

//class template instantiations
template class RelaxedDofUpdater<UblasSpace<double, CompressedMatrix, Vector>>;
}; // namespace Kratos

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block