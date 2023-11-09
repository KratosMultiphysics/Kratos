// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Ignasi de Pouplana,
//                   Vahid Galavi
//

#pragma once

#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"

#include "geo_mechanics_application_variables.h"
#include "geo_mechanics_scheme.hpp"

namespace Kratos
{

template<class TSparseSpace, class TDenseSpace>
class NewmarkQuasistaticUPwScheme : public GeoMechanicsScheme<TSparseSpace,TDenseSpace>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION( NewmarkQuasistaticUPwScheme );

    using BaseType              = Scheme<TSparseSpace,TDenseSpace>;
    using DofsArrayType         = typename BaseType::DofsArrayType;
    using TSystemMatrixType     = typename BaseType::TSystemMatrixType;
    using TSystemVectorType     = typename BaseType::TSystemVectorType;
    using LocalSystemVectorType = typename BaseType::LocalSystemVectorType;
    using LocalSystemMatrixType = typename BaseType::LocalSystemMatrixType;

    NewmarkQuasistaticUPwScheme(double beta, double gamma, double theta)
        : GeoMechanicsScheme<TSparseSpace, TDenseSpace>()
        , mBeta(beta)
        , mGamma(gamma)
        , mTheta(theta)
    {}

    int Check(const ModelPart& rModelPart) const override
    {
        KRATOS_TRY

        BaseType::Check(rModelPart);

        //check that variables are correctly allocated
        for (const auto& rNode : rModelPart.Nodes()) {
            KRATOS_ERROR_IF_NOT(rNode.SolutionStepsDataHas(DISPLACEMENT))
                << "DISPLACEMENT variable is not allocated for node "
                << rNode.Id()
                << std::endl;

            KRATOS_ERROR_IF_NOT(rNode.SolutionStepsDataHas(VELOCITY))
                << "VELOCITY variable is not allocated for node "
                << rNode.Id()
                << std::endl;

            KRATOS_ERROR_IF_NOT(rNode.SolutionStepsDataHas(ACCELERATION))
                << "ACCELERATION variable is not allocated for node "
                << rNode.Id()
                << std::endl;

            KRATOS_ERROR_IF_NOT(rNode.SolutionStepsDataHas(WATER_PRESSURE))
                << "WATER_PRESSURE variable is not allocated for node "
                << rNode.Id()
                << std::endl;

            KRATOS_ERROR_IF_NOT(rNode.SolutionStepsDataHas(DT_WATER_PRESSURE))
                << "DT_WATER_PRESSURE variable is not allocated for node "
                << rNode.Id()
                << std::endl;

            KRATOS_ERROR_IF_NOT(rNode.HasDofFor(DISPLACEMENT_X))
                << "missing DISPLACEMENT_X dof on node "
                << rNode.Id()
                << std::endl;

            KRATOS_ERROR_IF_NOT(rNode.HasDofFor(DISPLACEMENT_Y))
                << "missing DISPLACEMENT_Y dof on node "
                << rNode.Id()
                << std::endl;

            KRATOS_ERROR_IF_NOT(rNode.HasDofFor(DISPLACEMENT_Z))
                << "missing DISPLACEMENT_Z dof on node "
                << rNode.Id()
                << std::endl;

            KRATOS_ERROR_IF_NOT(rNode.HasDofFor(WATER_PRESSURE))
                << "missing WATER_PRESSURE dof on node "
                << rNode.Id()
                << std::endl;
        }

        CheckBufferSize(rModelPart);

        // Check beta, gamma and theta
        KRATOS_ERROR_IF(mBeta <= 0.0 || mGamma<= 0.0 || mTheta <= 0.0)
            << "Some of the scheme variables: beta, gamma or theta has an invalid value"
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

        if (rModelPart.GetProcessInfo()[NODAL_SMOOTHING]) {
            unsigned int Dim = rModelPart.GetProcessInfo()[DOMAIN_SIZE];

            SizeType StressTensorSize = STRESS_TENSOR_SIZE_2D;
            if (Dim == N_DIM_3D) StressTensorSize = STRESS_TENSOR_SIZE_3D; 

            // Clear nodal variables
            block_for_each(rModelPart.Nodes(), [&StressTensorSize](Node& rNode) {
                rNode.FastGetSolutionStepValue(NODAL_AREA) = 0.0;
                Matrix& rNodalStress = rNode.FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR);
                if (rNodalStress.size1() != StressTensorSize)
                    rNodalStress.resize(StressTensorSize,StressTensorSize,false);
                noalias(rNodalStress) = ZeroMatrix(StressTensorSize, StressTensorSize);
                rNode.FastGetSolutionStepValue(NODAL_DAMAGE_VARIABLE) = 0.0;
                rNode.FastGetSolutionStepValue(NODAL_JOINT_AREA)      = 0.0;
                rNode.FastGetSolutionStepValue(NODAL_JOINT_WIDTH)     = 0.0;
                rNode.FastGetSolutionStepValue(NODAL_JOINT_DAMAGE)    = 0.0;
            });

            FinalizeSolutionStepActiveEntities(rModelPart,A,Dx,b);

            // Compute smoothed nodal variables
            block_for_each(rModelPart.Nodes(), [&](Node& rNode) {
                const double& NodalArea = rNode.FastGetSolutionStepValue(NODAL_AREA);
                if (NodalArea > 1.0e-20) {
                    const double InvNodalArea = 1.0/NodalArea;
                    Matrix& rNodalStress = rNode.FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR);
                    for(unsigned int i = 0; i < rNodalStress.size1(); ++i) {
                        for(unsigned int j = 0; j < rNodalStress.size2(); ++j) {
                            rNodalStress(i,j) *= InvNodalArea;
                        }
                    }
                    rNode.FastGetSolutionStepValue(NODAL_DAMAGE_VARIABLE) *= InvNodalArea;
                }

                const double& NodalJointArea = rNode.FastGetSolutionStepValue(NODAL_JOINT_AREA);
                if (NodalJointArea > 1.0e-20) {
                    const double InvNodalJointArea = 1.0/NodalJointArea;
                    rNode.FastGetSolutionStepValue(NODAL_JOINT_WIDTH)  *= InvNodalJointArea;
                    rNode.FastGetSolutionStepValue(NODAL_JOINT_DAMAGE) *= InvNodalJointArea;
                }
            });

        } else {
            FinalizeSolutionStepActiveEntities(rModelPart,A,Dx,b);
        }

        KRATOS_CATCH("")
    }

protected:
    double mBeta  = 0.25;
    double mGamma = 0.5;
    double mTheta = 0.5;

    virtual inline void SetTimeFactors(ModelPart& rModelPart)
    {
        KRATOS_TRY

        SetDeltaTime(rModelPart.GetProcessInfo()[DELTA_TIME]);
        rModelPart.GetProcessInfo()[VELOCITY_COEFFICIENT]    = mGamma/(mBeta*GetDeltaTime());
        rModelPart.GetProcessInfo()[DT_PRESSURE_COEFFICIENT] = 1.0/(mTheta*GetDeltaTime());

        KRATOS_CATCH("")
    }

    virtual inline void UpdateVariablesDerivatives(ModelPart& rModelPart)
    {
        KRATOS_TRY

        //Update Acceleration, Velocity and DtPressure
        block_for_each(rModelPart.Nodes(), [&](Node& rNode) {
            noalias(rNode.FastGetSolutionStepValue(ACCELERATION)) =  ((  rNode.FastGetSolutionStepValue(DISPLACEMENT)
                                                                      - rNode.FastGetSolutionStepValue(DISPLACEMENT, 1))
                                                                   - GetDeltaTime() * rNode.FastGetSolutionStepValue(VELOCITY, 1)
                                                                   - (0.5-mBeta) * GetDeltaTime() * GetDeltaTime()
                                                                   * rNode.FastGetSolutionStepValue(ACCELERATION, 1) )
                                                                   / (mBeta*GetDeltaTime()*GetDeltaTime());

            noalias(rNode.FastGetSolutionStepValue(VELOCITY)) =  rNode.FastGetSolutionStepValue(VELOCITY, 1)
                                                                + (1.0-mGamma)*GetDeltaTime()
                                                                * rNode.FastGetSolutionStepValue(ACCELERATION, 1)
                                                                + mGamma * GetDeltaTime()
                                                                * rNode.FastGetSolutionStepValue(ACCELERATION);

            const double DeltaPressure =  rNode.FastGetSolutionStepValue(WATER_PRESSURE)
                                        - rNode.FastGetSolutionStepValue(WATER_PRESSURE, 1);

            rNode.FastGetSolutionStepValue(DT_WATER_PRESSURE) = ( DeltaPressure
                                                                 - (1.0-mTheta) * GetDeltaTime()
                                                                  * rNode.FastGetSolutionStepValue(DT_WATER_PRESSURE, 1))
                                                                / (mTheta*GetDeltaTime());
        });

        KRATOS_CATCH( "" )
    }
}; // Class NewmarkQuasistaticUPwScheme

}