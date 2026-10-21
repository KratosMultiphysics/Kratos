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

namespace Kratos
{

template <class TSparseSpace, class TDenseSpace>
class GeoMechanicsResidualCriteria : public ResidualCriteria<TSparseSpace, TDenseSpace>
{
public:
    using BaseType = ResidualCriteria<TSparseSpace, TDenseSpace>;
    using BaseType::BaseType;
    using DofsArrayType = typename BaseType::DofsArrayType;

protected:
    void ComputeActiveDofs(ModelPart& rModelPart, const DofsArrayType& rDofSet) override
    {
        this->mActiveDofs.assign(rDofSet.size(), 1);

        for (const auto& r_dof : rDofSet) {
            if (r_dof.IsFixed()) {
                this->mActiveDofs[r_dof.EquationId()] = 0;
            }
        }

        // Master DoFs remain independent; only active slave DoFs are eliminated.
        for (const auto& r_constraint : rModelPart.MasterSlaveConstraints()) {
            if (r_constraint.IsActive()) {
                for (const auto& p_slave_dof : r_constraint.GetSlaveDofsVector()) {
                    this->mActiveDofs[p_slave_dof->EquationId()] = 0;
                }
            }
        }
    }
};

} // namespace Kratos
