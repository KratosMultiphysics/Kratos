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
//

#if !defined(KRATOS_NEWMARK_QUASISTATIC_T_SCHEME )
#define  KRATOS_NEWMARK_QUASISTATIC_T_SCHEME

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

    typedef Scheme<TSparseSpace,TDenseSpace>                      BaseType;
    typedef typename BaseType::DofsArrayType                 DofsArrayType;
    typedef typename BaseType::TSystemMatrixType         TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType         TSystemVectorType;
    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;
    typedef NewmarkQuasistaticUPwScheme<TSparseSpace,TDenseSpace> MotherType;
    using MotherType::mDeltaTime;
    using MotherType::mBeta;
    using MotherType::mGamma;
    using MotherType::mTheta;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ///Constructor
    NewmarkQuasistaticTScheme(double theta) : 
        NewmarkQuasistaticUPwScheme<TSparseSpace,TDenseSpace>(0.25, 0.5, theta)
    { }
	
    ///Constructor
    NewmarkQuasistaticTScheme(double beta, double gamma, double theta) : Scheme<TSparseSpace,TDenseSpace>()
    {
        mBeta = beta;
        mGamma = gamma;
        mTheta = theta;
    }

    //------------------------------------------------------------------------------------

    ///Destructor
    ~NewmarkQuasistaticTScheme() override {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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

            //if (rNode.SolutionStepsDataHas(DT_WATER_PRESSURE) == false)
            //    KRATOS_ERROR << "DT_WATER_PRESSURE variable is not allocated for node "
            //                 << rNode.Id()
            //                 << std::endl;

            if (rNode.HasDofFor(TEMPERATURE) == false)
                KRATOS_ERROR << "missing TEMPERATURE dof on node "
                             << rNode.Id()
                             << std::endl;
        }

        //check for minimum value of the buffer index.
        if (rModelPart.GetBufferSize() < 2)
            KRATOS_ERROR << "insufficient buffer size. Buffer size should be greater than 2. Current size is "
                         << rModelPart.GetBufferSize()
                         << std::endl;

        // Check beta, gamma and theta
        if (mBeta <= 0.0 || mGamma<= 0.0 || mTheta <= 0.0)
            KRATOS_ERROR << "Some of the scheme variables: beta, gamma or theta has an invalid value "
                         << std::endl;

        return 0;

        KRATOS_CATCH( "" )
    }


//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
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


//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    inline void UpdateVariablesDerivatives(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        //Update DtPressure
        block_for_each(rModelPart.Nodes(), [&](Node<3>& rNode){
            const double DeltaTemperature =  rNode.FastGetSolutionStepValue(TEMPERATURE)
                                        - rNode.FastGetSolutionStepValue(TEMPERATURE, 1);
            //const auto &PreviousDtPressure = rNode.FastGetSolutionStepValue(DT_WATER_PRESSURE, 1);

            //rNode.FastGetSolutionStepValue(DT_WATER_PRESSURE) =  1.0/(mTheta*mDeltaTime)*(DeltaPressure - (1.0-mTheta)*mDeltaTime*PreviousDtPressure);
        });

        KRATOS_CATCH( "" )
    }

}; // Class NewmarkQuasistaticTScheme
}  // namespace Kratos

#endif // KRATOS_NEWMARK_QUASISTATIC_T_SCHEME defined
