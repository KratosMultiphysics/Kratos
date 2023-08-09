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
#include "geo_mechanics_application_variables.h"
#include "custom_strategies/schemes/newmark_quasistatic_U_Pw_scheme.hpp"

namespace Kratos
{

template<class TSparseSpace, class TDenseSpace>
class NewmarkQuasistaticPwScheme : public NewmarkQuasistaticUPwScheme<TSparseSpace,TDenseSpace>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION( NewmarkQuasistaticPwScheme );

    using BaseType              = Scheme<TSparseSpace,TDenseSpace>;
    using TSystemMatrixType     = typename BaseType::TSystemMatrixType;
    using TSystemVectorType     = typename BaseType::TSystemVectorType;
    using MotherType = NewmarkQuasistaticUPwScheme<TSparseSpace,TDenseSpace>;
    using MotherType::mDeltaTime;
// mBeta and mGamma are not really used
    using MotherType::mBeta;
    using MotherType::mGamma;
    using MotherType::mTheta;

    explicit NewmarkQuasistaticPwScheme(double theta) :
        NewmarkQuasistaticUPwScheme<TSparseSpace,TDenseSpace>(0.25, 0.5, theta)
    {}

    int Check(const ModelPart& rModelPart) const override
    {
        KRATOS_TRY

        BaseType::Check(rModelPart);

        //check that variables are correctly allocated
        for (const auto& rNode : rModelPart.Nodes())
        {
            KRATOS_ERROR_IF_NOT(rNode.SolutionStepsDataHas(WATER_PRESSURE))
                << "WATER_PRESSURE variable is not allocated for node "
                << rNode.Id()
                << std::endl;

            KRATOS_ERROR_IF_NOT(rNode.SolutionStepsDataHas(DT_WATER_PRESSURE))
                << "DT_WATER_PRESSURE variable is not allocated for node "
                << rNode.Id()
                << std::endl;

            KRATOS_ERROR_IF_NOT(rNode.HasDofFor(WATER_PRESSURE))
                << "missing WATER_PRESSURE dof on node "
                << rNode.Id()
                << std::endl;
        }

        //check for minimum value of the buffer index.
        KRATOS_ERROR_IF(rModelPart.GetBufferSize() < 2)
            << "insufficient buffer size. Buffer size should be greater than 2. Current size is "
            << rModelPart.GetBufferSize()
            << std::endl;

        // Check beta, gamma and theta
        KRATOS_ERROR_IF(mBeta <= 0.0 || mGamma<= 0.0 || mTheta <= 0.0)
            << "Some of the scheme variables: beta, gamma or theta has an invalid value "
            << std::endl;

        return 0;

        KRATOS_CATCH( "" )
    }

    void FinalizeSolutionStep( ModelPart& rModelPart,
                               TSystemMatrixType& A,
                               TSystemVectorType& Dx,
                               TSystemVectorType& b) override
    {
        KRATOS_TRY

        MotherType::FinalizeSolutionStepActiveEntities(rModelPart,A,Dx,b);

        KRATOS_CATCH("")
    }

protected:
    inline void UpdateVariablesDerivatives(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        //Update DtPressure
        block_for_each(rModelPart.Nodes(), [&](Node& rNode){
            const double DeltaPressure =  rNode.FastGetSolutionStepValue(WATER_PRESSURE)
                                        - rNode.FastGetSolutionStepValue(WATER_PRESSURE, 1);
            const auto &PreviousDtPressure = rNode.FastGetSolutionStepValue(DT_WATER_PRESSURE, 1);

            rNode.FastGetSolutionStepValue(DT_WATER_PRESSURE) =  (1.0/(mTheta*mDeltaTime))*(DeltaPressure - (1.0-mTheta)*mDeltaTime*PreviousDtPressure);
        });

        KRATOS_CATCH( "" )
    }
}; // Class NewmarkQuasistaticPwScheme

}