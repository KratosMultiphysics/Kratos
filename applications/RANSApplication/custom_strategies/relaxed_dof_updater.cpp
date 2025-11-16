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
    mOldRelaxationFactor = 0.0;
}

template <>
void RelaxedDofUpdater<UblasSpace<double, CompressedMatrix, Vector>>::Clear()
{
}

template <>
std::string RelaxedDofUpdater<UblasSpace<double, CompressedMatrix, Vector>>::Info() const
{
    std::stringstream buffer;
    buffer << "RelaxedDofUpdater - UblasSpace";
    return buffer.str();
}

template <>
void RelaxedDofUpdater<UblasSpace<double, CompressedMatrix, Vector>>::UpdateDofs(
    DofsArrayType& rDofSet,
    const SystemVectorType& rDx)
{

    double relaxation_factor = mRelaxationFactor;
    if (mRelaxationFactor <= 0.0) {
        if (mOldDx.size() != rDx.size() || mOldRelaxationFactor == 0.0) {
            mOldDx.resize(rDx.size());
            mDiffDx.resize(rDx.size());

            relaxation_factor = mInitialRelaxationFactor;
        } else {
            noalias(mDiffDx) = rDx - mOldDx;
            const double numerator = inner_prod(mOldDx, mDiffDx);
            const double denominator = inner_prod(mDiffDx, mDiffDx);
            relaxation_factor = std::max(
                std::min(-mOldRelaxationFactor * numerator / denominator, mMaxRelaxationFactor),
                mMinRelaxationFactor);
        }

        noalias(mOldDx) = rDx;
        mOldRelaxationFactor = relaxation_factor;
        KRATOS_INFO(this->Info()) << "Using relaxation factor " << relaxation_factor << ".\n";
    }

    block_for_each(rDofSet, [&](DofType& rDof) {
        if (rDof.IsFree()) {
            rDof.GetSolutionStepValue() += rDx[rDof.EquationId()] * relaxation_factor;
        }
    });
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