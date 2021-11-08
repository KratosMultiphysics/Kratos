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

#if !defined(KRATOS_NEWMARK_QUASISTATIC_PW_SCHEME )
#define  KRATOS_NEWMARK_QUASISTATIC_PW_SCHEME

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"

// Application includes
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

template<class TSparseSpace, class TDenseSpace>

class NewmarkQuasistaticPwScheme : public NewmarkQuasistaticUPwScheme<TSparseSpace,TDenseSpace>
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( NewmarkQuasistaticPwScheme );

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
    NewmarkQuasistaticPwScheme(double theta) : 
        NewmarkQuasistaticUPwScheme<TSparseSpace,TDenseSpace>(0.25, 0.5, theta)
    { }

    //------------------------------------------------------------------------------------

    ///Destructor
    ~NewmarkQuasistaticPwScheme() override {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    int Check(const ModelPart& rModelPart) const override
    {
        KRATOS_TRY

        BaseType::Check(rModelPart);

        //check that variables are correctly allocated
        for (const auto& rNode : rModelPart.Nodes())
        {
            if (rNode.SolutionStepsDataHas(WATER_PRESSURE) == false)
                KRATOS_ERROR << "WATER_PRESSURE variable is not allocated for node "
                             << rNode.Id()
                             << std::endl;

            if (rNode.SolutionStepsDataHas(DT_WATER_PRESSURE) == false)
                KRATOS_ERROR << "DT_WATER_PRESSURE variable is not allocated for node "
                             << rNode.Id()
                             << std::endl;

            if (rNode.HasDofFor(WATER_PRESSURE) == false)
                KRATOS_ERROR << "missing WATER_PRESSURE dof on node "
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

        //Update Acceleration, Velocity and DtPressure

        double DeltaPressure;

        const int NNodes = static_cast<int>(rModelPart.Nodes().size());
        ModelPart::NodesContainerType::iterator node_begin = rModelPart.NodesBegin();

        #pragma omp parallel for private(DeltaPressure)
        for(int i = 0; i < NNodes; i++)
        {
            ModelPart::NodesContainerType::iterator itNode = node_begin + i;

            double& CurrentDtPressure = itNode->FastGetSolutionStepValue(DT_WATER_PRESSURE);
            DeltaPressure = itNode->FastGetSolutionStepValue(WATER_PRESSURE) - itNode->FastGetSolutionStepValue(WATER_PRESSURE, 1);
            const double& PreviousDtPressure = itNode->FastGetSolutionStepValue(DT_WATER_PRESSURE, 1);

            CurrentDtPressure = 1.0/(mTheta*mDeltaTime)*(DeltaPressure - (1.0-mTheta)*mDeltaTime*PreviousDtPressure);
        }

        KRATOS_CATCH( "" )
    }

}; // Class NewmarkQuasistaticPwScheme
}  // namespace Kratos

#endif // KRATOS_NEWMARK_QUASISTATIC_PW_SCHEME defined
