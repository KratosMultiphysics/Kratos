// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Mohamed Nabi
//                   John van Esch
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

class NewmarkQuasistaticTScheme : public NewmarkQuasistaticUPwScheme<TSparseSpace,TDenseSpace>
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( NewmarkQuasistaticTScheme );

    using BaseType = Scheme<TSparseSpace, TDenseSpace>;
    using DofsArrayType = typename BaseType::DofsArrayType;
    using TSystemMatrixType = typename BaseType::TSystemMatrixType;
    using TSystemVectorType = typename BaseType::TSystemVectorType;
    using LocalSystemVectorType = typename BaseType::LocalSystemVectorType;
    using LocalSystemMatrixType = typename BaseType::LocalSystemMatrixType;
    using MotherType = NewmarkQuasistaticUPwScheme<TSparseSpace, TDenseSpace>;
    using MotherType::mDeltaTime;
    using MotherType::mBeta;
    using MotherType::mGamma;
    using MotherType::mTheta;

    ///Constructor
    // ============================================================================================
    // ============================================================================================
    explicit NewmarkQuasistaticTScheme(double theta) : 
        NewmarkQuasistaticUPwScheme<TSparseSpace,TDenseSpace>(0.25, 0.5, theta)
    { }

    ///Destructor
    // ============================================================================================
    // ============================================================================================
    ~NewmarkQuasistaticTScheme() override = default;

    // ============================================================================================
    // ============================================================================================
    int Check(const ModelPart& rModelPart) const override
    {
        KRATOS_TRY

        BaseType::Check(rModelPart);

        //check that variables are correctly allocated
        for (const auto& rNode : rModelPart.Nodes())
        {
            if (rNode.SolutionStepsDataHas(TEMPERATURE) == false)
                KRATOS_ERROR << "TEMPERATURE variable is not allocated for node "
                             << rNode.Id()
                             << std::endl;

            if (rNode.SolutionStepsDataHas(DT_TEMPERATURE) == false)
                KRATOS_ERROR << "DT_TEMPERATURE variable is not allocated for node "
                             << rNode.Id()
                             << std::endl;

            if (rNode.HasDofFor(TEMPERATURE) == false)
                KRATOS_ERROR << "missing TEMPERATURE dof on node "
                             << rNode.Id()
                             << std::endl;
        }

        //check for minimum value of the buffer index.
        if (rModelPart.GetBufferSize() < 2)
            KRATOS_ERROR << "insufficient buffer size. Buffer size should be greater or equal to 2. Current size is "
                         << rModelPart.GetBufferSize()
                         << std::endl;

        // Check beta, gamma and theta
        if (mBeta <= 0.0 || mGamma<= 0.0 || mTheta <= 0.0)
            KRATOS_ERROR << "Some of the scheme variables: beta, gamma or theta has an invalid value "
                         << std::endl;

        return 0;

        KRATOS_CATCH( "" )
    }


    // ============================================================================================
    // ============================================================================================
    void FinalizeSolutionStep(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        MotherType::FinalizeSolutionStepActiveEntities(rModelPart,A,Dx,b);

        KRATOS_CATCH("")
    }


protected:

    /// Member Variables
    // ============================================================================================
    // ============================================================================================
    inline void UpdateVariablesDerivatives(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        //Update DtPressure
        block_for_each(rModelPart.Nodes(), [this](Node& rNode){
            const double DeltaTemperature =  rNode.FastGetSolutionStepValue(TEMPERATURE)
                                        - rNode.FastGetSolutionStepValue(TEMPERATURE, 1);
            const auto &PreviousDtTemperature = rNode.FastGetSolutionStepValue(DT_TEMPERATURE, 1);

            rNode.FastGetSolutionStepValue(DT_TEMPERATURE) =  1.0/(mTheta*mDeltaTime)*(DeltaTemperature - (1.0-mTheta)*mDeltaTime*PreviousDtTemperature);
        });

        KRATOS_CATCH( "" )
    }

    // ============================================================================================
    // ============================================================================================
    inline void SetTimeFactors(ModelPart& rModelPart) override
    {
        KRATOS_TRY

    	mDeltaTime = rModelPart.GetProcessInfo()[DELTA_TIME];
        rModelPart.GetProcessInfo()[DT_TEMPERATURE_COEFFICIENT] = 1.0 / (mTheta * mDeltaTime);

        KRATOS_CATCH("")
    }

}; // Class NewmarkQuasistaticTScheme
}  // namespace Kratos