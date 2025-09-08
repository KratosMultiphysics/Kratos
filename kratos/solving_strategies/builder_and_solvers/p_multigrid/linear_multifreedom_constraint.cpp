//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Máté Kelemen
//

// Project includes
#include "solving_strategies/builder_and_solvers/p_multigrid/linear_multifreedom_constraint.hpp" // MultifreedomConstraint
#include "includes/variables.h" // GEOMETRIC_STIFFNESS_MATRIX, INTERNAL_FORCES_VECTOR


namespace Kratos {


LinearMultifreedomConstraint::LinearMultifreedomConstraint(const IndexType Id,
                                                           DofPointerVectorType&& rDofs,
                                                           const std::vector<std::size_t>& rConstraintLabels,
                                                           const MatrixType& rRelationMatrix,
                                                           const VectorType& rConstraintGaps)
    : MultifreedomConstraint(Id, std::move(rDofs), rConstraintLabels)
{
    // Sanity checks.
    KRATOS_ERROR_IF_NOT(rRelationMatrix.size1() == rConstraintGaps.size())
        << "relation matrix (" << rRelationMatrix.size1() << "x" << rRelationMatrix.size2() << ") "
        << "is incompatible with its constraint gap vector (" << rConstraintGaps.size() << ")";

    KRATOS_ERROR_IF_NOT(this->GetDofs().size() == rRelationMatrix.size2())
        << "provided number of DoFs (" << this->GetDofs().size() << ") "
        << "is incompatible with the relation matrix "
        << "(" << rRelationMatrix.size1() << "x" << rRelationMatrix.size2() << ")";

    KRATOS_ERROR_IF_NOT(rConstraintLabels.size() == rRelationMatrix.size1())
        << "constraint label vector (" << rConstraintLabels.size() << ") "
        << "in incompatible with the relation matrix "
        << "(" << rRelationMatrix.size1() << "x" << rRelationMatrix.size2() << ")";

    // Store input data.
    this->SetValue(GEOMETRIC_STIFFNESS_MATRIX, rRelationMatrix);
    this->SetValue(INTERNAL_FORCES_VECTOR, rConstraintGaps);
}


void LinearMultifreedomConstraint::CalculateLocalSystem(MatrixType& rRelationMatrix,
                                                        VectorType& rConstraintGaps,
                                                        const ProcessInfo&) const
{
    rRelationMatrix = this->GetData().GetValue(GEOMETRIC_STIFFNESS_MATRIX);
    rConstraintGaps = this->GetData().GetValue(INTERNAL_FORCES_VECTOR);
}


} // namespace Kratos
