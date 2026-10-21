// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Gennady Markelov
//

#pragma once

#include "solving_strategies/convergencecriterias/residual_criteria.h"

#include <utility>

namespace Kratos
{

template <class TSparseSpace, class TDenseSpace>
class GeoMechanicsResidualCriteria : public ResidualCriteria<TSparseSpace, TDenseSpace>
{
public:
    using BaseType                = ResidualCriteria<TSparseSpace, TDenseSpace>;
    using ConvergenceCriteriaType = typename BaseType::BaseType;
    using DofsArrayType           = typename BaseType::DofsArrayType;
    using SizeType                = typename BaseType::SizeType;
    using TDataType               = typename BaseType::TDataType;
    using TSystemVectorType       = typename BaseType::TSystemVectorType;

    enum class ResidualRepresentation { Raw, ConstraintReduced };

    explicit GeoMechanicsResidualCriteria(ResidualRepresentation Representation = ResidualRepresentation::Raw)
        : BaseType(), mResidualRepresentation(Representation)
    {
    }

    explicit GeoMechanicsResidualCriteria(Parameters ThisParameters,
                                          ResidualRepresentation Representation = ResidualRepresentation::Raw)
        : BaseType(std::move(ThisParameters)), mResidualRepresentation(Representation)
    {
    }

    typename ConvergenceCriteriaType::Pointer Create(Parameters ThisParameters) const override
    {
        return Kratos::make_shared<GeoMechanicsResidualCriteria>(std::move(ThisParameters), mResidualRepresentation);
    }

    ResidualRepresentation GetResidualRepresentation() const { return mResidualRepresentation; }

protected:
    void ComputeActiveDofs(ModelPart& rModelPart, const DofsArrayType& rDofSet) override
    {
        this->mActiveDofs.assign(rDofSet.size(), 1);

        for (const auto& r_dof : rDofSet) {
            if (r_dof.IsFixed()) {
                this->mActiveDofs[r_dof.EquationId()] = 0;
            }
        }

        // Raw residuals include both sides of an MPC; reduced residuals retain the independent master equations.
        for (const auto& r_constraint : rModelPart.MasterSlaveConstraints()) {
            if (r_constraint.IsActive()) {
                if (mResidualRepresentation == ResidualRepresentation::Raw) {
                    for (const auto& p_master_dof : r_constraint.GetMasterDofsVector()) {
                        this->mActiveDofs[p_master_dof->EquationId()] = 0;
                    }
                }
                for (const auto& p_slave_dof : r_constraint.GetSlaveDofsVector()) {
                    this->mActiveDofs[p_slave_dof->EquationId()] = 0;
                }
            }
        }
    }

    void CalculateResidualNorm(ModelPart&               rModelPart,
                               TDataType&               rResidualSolutionNorm,
                               SizeType&                rDofNum,
                               DofsArrayType&           rDofSet,
                               const TSystemVectorType& rb) override
    {
        BaseType::CalculateResidualNorm(rModelPart, rResidualSolutionNorm, rDofNum, rDofSet, rb);
        if (rDofNum == 0) {
            // An empty active set has a zero norm; use a nonzero divisor in ResidualCriteria::PostCriteria.
            rDofNum = 1;
        }
    }

private:
    ResidualRepresentation mResidualRepresentation;
};

} // namespace Kratos
