// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//

#pragma once

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"

// Application includes
#include "custom_strategies/schemes/newmark_quasistatic_U_Pw_scheme.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos {

template <class TSparseSpace, class TDenseSpace>
class NewmarkQuasistaticPwScheme
    : public NewmarkQuasistaticUPwScheme<TSparseSpace, TDenseSpace> {
public:
    KRATOS_CLASS_POINTER_DEFINITION(NewmarkQuasistaticPwScheme);

    using BaseType = Scheme<TSparseSpace, TDenseSpace>;
    using TSystemMatrixType = typename BaseType::TSystemMatrixType;
    using TSystemVectorType = typename BaseType::TSystemVectorType;
    using MotherType = NewmarkQuasistaticUPwScheme<TSparseSpace, TDenseSpace>;
    // mBeta and mGamma are not really used
    using MotherType::mBeta;
    using MotherType::mGamma;
    using MotherType::mTheta;

    explicit NewmarkQuasistaticPwScheme(double theta)
        : NewmarkQuasistaticUPwScheme<TSparseSpace, TDenseSpace>(0.25, 0.5, theta)
    {
    }

    int Check(const ModelPart& rModelPart) const override
    {
        KRATOS_TRY

        GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>::Check(rModelPart);

        // Check beta, gamma and theta
        KRATOS_ERROR_IF(mBeta <= 0.0 || mGamma <= 0.0 || mTheta <= 0.0)
            << "Some of the scheme variables: beta, gamma or theta has an "
               "invalid value "
            << std::endl;

        return 0;

        KRATOS_CATCH("")
    }

protected:
    inline void UpdateVariablesDerivatives(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        // Update DtPressure
        block_for_each(rModelPart.Nodes(), [&](Node& rNode) {
            this->UpdateScalarTimeDerivative(rNode, WATER_PRESSURE, DT_WATER_PRESSURE);
        });

        KRATOS_CATCH("")
    }

    void CheckAllocatedVariables(const ModelPart& rModelPart) const override
    {
        // check that variables are correctly allocated
        for (const auto& rNode : rModelPart.Nodes()) {
            this->CheckSolutionStepsData(rNode, WATER_PRESSURE);
            this->CheckSolutionStepsData(rNode, DT_WATER_PRESSURE);
            this->CheckDof(rNode, WATER_PRESSURE);
        }
    }

}; // Class NewmarkQuasistaticPwScheme

} // namespace Kratos