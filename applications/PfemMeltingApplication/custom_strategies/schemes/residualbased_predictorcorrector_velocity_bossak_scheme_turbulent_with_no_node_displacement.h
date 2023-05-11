//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//


#if !defined(KRATOS_RESIDUALBASED_PREDICTOR_CORRECTOR_VELOCITY_BOSSAK_TURBULENT_SCHEME_WITH_NO_NODE_DISPLACEMENT )
#define  KRATOS_RESIDUALBASED_PREDICTOR_CORRECTOR_VELOCITY_BOSSAK_TURBULENT_SCHEME_WITH_NO_NODE_DISPLACEMENT


/* System includes */


/* External includes */
#include "boost/smart_ptr.hpp"

/* Project includes */
#include "includes/define.h"

#include "../../../FluidDynamicsApplication/custom_strategies/schemes/residualbased_predictorcorrector_velocity_bossak_scheme_turbulent.h"

namespace Kratos {

    template<class TSparseSpace, class TDenseSpace>
    class ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulentWithNoNodeDisplacement : public ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent<TSparseSpace, TDenseSpace> {
    public:

        KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulentWithNoNodeDisplacement);

        typedef Scheme<TSparseSpace, TDenseSpace> BaseType;
        typedef typename BaseType::TDataType TDataType;
        typedef typename BaseType::DofsArrayType DofsArrayType;
        typedef typename Element::DofsVectorType DofsVectorType;
        typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
        typedef typename BaseType::TSystemVectorType TSystemVectorType;
        typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
        typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;
        typedef Element::GeometryType  GeometryType;

        ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulentWithNoNodeDisplacement(
            double NewAlphaBossak,
            double MoveMeshStrategy,
            unsigned int DomainSize)
        : ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent<TSparseSpace, TDenseSpace>(NewAlphaBossak,
                                                                                                MoveMeshStrategy,
                                                                                                DomainSize) {}
        ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulentWithNoNodeDisplacement(
            double NewAlphaBossak,
            unsigned int DomainSize,
            const Variable<int>& rPeriodicIdVar)
        : ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent<TSparseSpace, TDenseSpace>(NewAlphaBossak,
                                                                                                DomainSize,
                                                                                                rPeriodicIdVar) {}
        ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulentWithNoNodeDisplacement(
            double NewAlphaBossak,
            double MoveMeshStrategy,
            unsigned int DomainSize,
            Kratos::Flags& rSlipFlag)
        : ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent<TSparseSpace, TDenseSpace>(NewAlphaBossak,
                                                                                                MoveMeshStrategy,
                                                                                                DomainSize,
                                                                                                rSlipFlag) {}
        ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulentWithNoNodeDisplacement(
            double NewAlphaBossak,
            double MoveMeshStrategy,
            unsigned int DomainSize,
            Process::Pointer pTurbulenceModel)
        : ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent<TSparseSpace, TDenseSpace>(NewAlphaBossak,
                                                                                                MoveMeshStrategy,
                                                                                                DomainSize,
                                                                                                pTurbulenceModel) {}
        ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulentWithNoNodeDisplacement(
            double NewAlphaBossak,
            double MoveMeshStrategy,
            unsigned int DomainSize,
            const double RelaxationFactor,
            Process::Pointer pTurbulenceModel)
        :  ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent<TSparseSpace, TDenseSpace>(NewAlphaBossak,
                                                                                                MoveMeshStrategy,
                                                                                                DomainSize,
                                                                                                RelaxationFactor,
                                                                                                pTurbulenceModel) {}
        //*********************************************************************************
        //Updating first time Derivative
        //*********************************************************************************
        void UpdateDisplacement(array_1d<double, 3 > & CurrentDisplacement,
                                const array_1d<double, 3 > & OldDisplacement,
                                const array_1d<double, 3 > & OldVelocity,
                                const array_1d<double, 3 > & OldAcceleration,
                                const array_1d<double, 3 > & CurrentAcceleration) override
        {
            return;
        }
    };

} /* namespace Kratos.*/

#endif /* KRATOS_RESIDUALBASED_PREDICTOR_CORRECTOR_VELOCITY_BOSSAK_TURBULENT_SCHEME_WITH_NO_NODE_DISPLACEMENT  defined */
