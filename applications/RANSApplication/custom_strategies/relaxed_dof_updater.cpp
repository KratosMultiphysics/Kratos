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
#include "relaxed_dof_updater.h"
#include "spaces/ublas_space.h"

namespace Kratos
{
template <>
void RelaxedDofUpdater<UblasSpace<double, CompressedMatrix, Vector>>::Initialize(
    const RelaxedDofUpdater<UblasSpace<double, CompressedMatrix, Vector>>::DofsArrayType& rDofSet,
    const RelaxedDofUpdater<UblasSpace<double, CompressedMatrix, Vector>>::SystemVectorType& rDx)
{
}

template <>
void RelaxedDofUpdater<UblasSpace<double, CompressedMatrix, Vector>>::Clear()
{
}

template <>
void RelaxedDofUpdater<UblasSpace<double, CompressedMatrix, Vector>>::UpdateDofs(
    RelaxedDofUpdater<UblasSpace<double, CompressedMatrix, Vector>>::DofsArrayType& rDofSet,
    const RelaxedDofUpdater<UblasSpace<double, CompressedMatrix, Vector>>::SystemVectorType& rDx,
    const double RelaxationFactor)
{
    const int num_dof = static_cast<int>(rDofSet.size());

#pragma omp parallel for
    for (int i = 0; i < num_dof; ++i)
    {
        auto it_dof = rDofSet.begin() + i;

        if (it_dof->IsFree())
            it_dof->GetSolutionStepValue() += rDx[it_dof->EquationId()] * RelaxationFactor;
    }
}

template <>
void RelaxedDofUpdater<UblasSpace<double, CompressedMatrix, Vector>>::AssignDofs(
    RelaxedDofUpdater<UblasSpace<double, CompressedMatrix, Vector>>::DofsArrayType& rDofSet,
    const RelaxedDofUpdater<UblasSpace<double, CompressedMatrix, Vector>>::SystemVectorType& rX)
{
    const int num_dof = static_cast<int>(rDofSet.size());

#pragma omp parallel for
    for (int i = 0; i < num_dof; ++i)
    {
        auto it_dof = rDofSet.begin() + i;
        if (it_dof->IsFree())
            it_dof->GetSolutionStepValue() = rX[it_dof->EquationId()];
    }
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

