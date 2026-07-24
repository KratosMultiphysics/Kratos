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

using SparseSpace = UblasSpace<double, CompressedMatrix, Vector>;

template<class TSparseSpace>
void RelaxedDofUpdater<TSparseSpace>::UpdateDofs(
    DofsArrayType& rDofSet,
    const SystemVectorType& rDx)
{
    if constexpr (std::is_same_v<TSparseSpace, SparseSpace>) {
        block_for_each(rDofSet, [&](DofType& rDof) {
            if (rDof.IsFree()) {
                rDof.GetSolutionStepValue() += rDx[rDof.EquationId()] * mRelaxationFactor;
            }
        });
    }
}

template<class TSparseSpace>
std::string RelaxedDofUpdater<TSparseSpace>::Info() const
{
    if constexpr (std::is_same_v<TSparseSpace, SparseSpace>) {
        std::stringstream buffer;
        buffer << "RelaxedDofUpdater - UblasSpace";
        return buffer.str();
    }
}

///@}

//class template instantiations
template class DofUpdater<SparseSpace>;
template class RelaxedDofUpdater<SparseSpace>;
}; // namespace Kratos

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block