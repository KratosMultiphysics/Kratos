//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//                   Miguel Maso
//


#if !defined(KRATOS_RESIDUALBASED_PREDICTOR_CORRECTOR_VELOCITY_BOSSAK_ALE_SCHEME )
#define  KRATOS_RESIDUALBASED_PREDICTOR_CORRECTOR_VELOCITY_BOSSAK_ALE_SCHEME


/* System includes */


/* External includes */
#include "boost/smart_ptr.hpp"

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/deprecated_variables.h"
#include "solving_strategies/schemes/scheme.h"
#include "includes/variables.h"
#include "includes/cfd_variables.h"
#include "containers/array_1d.h"
#include "utilities/openmp_utils.h"
#include "utilities/dof_updater.h"
#include "utilities/coordinate_transformation_utilities.h"
#include "processes/process.h"

#include "../../applications/FluidDynamicsApplication/custom_strategies/schemes/residualbased_predictorcorrector_velocity_bossak_scheme_turbulent.h"

namespace Kratos {

    /**@name Kratos Globals */
    /*@{ */


    /*@} */
    /**@name Type Definitions */
    /*@{ */

    /*@} */


    /**@name  Enum's */
    /*@{ */


    /*@} */
    /**@name  Functions */
    /*@{ */



    /*@} */
    /**@name Kratos Classes */
    /*@{ */

    /// Bossak time scheme for the incompressible flow problem.

    // template<template <class TSparseSpace, class TDenseSpace> class TSchemeType >
    // class ResidualBasedPredictorCorrectorVelocityBossakAleScheme : public TSchemeType<TSparseSpace, TDenseSpace> {
    template<class TSparseSpace, class TDenseSpace >
    class ResidualBasedPredictorCorrectorVelocityBossakAleScheme : public ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent<TSparseSpace, TDenseSpace> {

    public:
        /**@name Type Definitions */
        /*@{ */

        KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedPredictorCorrectorVelocityBossakAleScheme);

        typedef Scheme<TSparseSpace, TDenseSpace> BaseType;

        typedef ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent<TSparseSpace, TDenseSpace> SchemeType;

        typedef typename BaseType::TDataType TDataType;

        typedef typename BaseType::DofsArrayType DofsArrayType;

        typedef typename Element::DofsVectorType DofsVectorType;

        typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

        typedef typename BaseType::TSystemVectorType TSystemVectorType;

        typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

        typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

        typedef Element::GeometryType  GeometryType;

        // using SchemeType::mRotationTool;


        /*@} */
        /**@name Life Cycle
         */
        /*@{ */

        /** Constructor without a turbulence model
         */
        ResidualBasedPredictorCorrectorVelocityBossakAleScheme(
            double NewAlphaBossak,
            double MoveMeshStrategy,
            unsigned int DomainSize)
        : ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent<TSparseSpace, TDenseSpace>(NewAlphaBossak,MoveMeshStrategy,DomainSize)
        {
        }


        /** Constructor without a turbulence model with periodic conditions
         */
        ResidualBasedPredictorCorrectorVelocityBossakAleScheme(
            double NewAlphaBossak,
            unsigned int DomainSize,
            const Variable<int>& rPeriodicIdVar)
        : ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent<TSparseSpace, TDenseSpace>(NewAlphaBossak,DomainSize,rPeriodicIdVar)
        {
        }


        /** Constructor without a turbulence model
         */
        ResidualBasedPredictorCorrectorVelocityBossakAleScheme(
            double NewAlphaBossak,
            double MoveMeshStrategy,
            unsigned int DomainSize,
            Variable<double>& rSlipVar)
        : ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent<TSparseSpace, TDenseSpace>(NewAlphaBossak,MoveMeshStrategy,rSlipVar)
        {
        }

        /** Constructor with a turbulence model
         */
        ResidualBasedPredictorCorrectorVelocityBossakAleScheme(
            double NewAlphaBossak,
            double MoveMeshStrategy,
            unsigned int DomainSize,
            Process::Pointer pTurbulenceModel)
        : ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent<TSparseSpace, TDenseSpace>(NewAlphaBossak,MoveMeshStrategy,DomainSize,pTurbulenceModel)
        {
        }

        /** Destructor.
         */
        ~ResidualBasedPredictorCorrectorVelocityBossakAleScheme() override {
        }


        /*@} */
        /**@name Operators
         */
        /*@{ */

        /**
                Performing the update of the solution.
         */
        //***************************************************************************

        void Update(ModelPart& r_model_part,
                            DofsArrayType& rDofSet,
                            TSystemMatrixType& A,
                            TSystemVectorType& Dv,
                            TSystemVectorType& b) override
        {
            KRATOS_TRY;

            SchemeType::Update(r_model_part,rDofSet,A,Dv,b);

            this->Pfem2AdditionalUpdateOperations(r_model_part, rDofSet, A, Dv, b);

            KRATOS_CATCH("")
        }

        //***************************************************************************

        void Pfem2AdditionalUpdateOperations(ModelPart& rModelPart,
                                        DofsArrayType& rDofSet,
                                        TSystemMatrixType& A,
                                        TSystemVectorType& Dv,
                                        TSystemVectorType& b)
        {
            KRATOS_TRY

            int NumThreads = OpenMPUtils::GetNumThreads();
            OpenMPUtils::PartitionVector NodePartition;
            OpenMPUtils::DivideInPartitions(rModelPart.Nodes().size(), NumThreads, NodePartition);

            //updating time derivatives (nodally for efficiency)
            #pragma omp parallel
            {
                int k = OpenMPUtils::ThisThread();

                ModelPart::NodeIterator NodesBegin = rModelPart.NodesBegin() + NodePartition[k];
                ModelPart::NodeIterator NodesEnd = rModelPart.NodesBegin() + NodePartition[k + 1];

                for (ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; itNode++) {

                    // Pfem2 ALE update to eliminate convective term
                    noalias(itNode->FastGetSolutionStepValue(MESH_VELOCITY)) = itNode->FastGetSolutionStepValue(VELOCITY);
                }
            }

            KRATOS_CATCH("")

        }


        //***************************************************************************
        //predicts the solution at the current step as
        // v = vold
        void Predict(ModelPart& rModelPart,
                             DofsArrayType& rDofSet,
                             TSystemMatrixType& A,
                             TSystemVectorType& Dv,
                             TSystemVectorType& b) override
        {
            // if (rModelPart.GetCommunicator().MyPID() == 0)
            //     std::cout << "prediction" << std::endl;

            int NumThreads = OpenMPUtils::GetNumThreads();
            OpenMPUtils::PartitionVector NodePartition;
            OpenMPUtils::DivideInPartitions(rModelPart.Nodes().size(), NumThreads, NodePartition);

            #pragma omp parallel
            {
                //array_1d<double, 3 > DeltaDisp;

                int k = OpenMPUtils::ThisThread();

                ModelPart::NodeIterator NodesBegin = rModelPart.NodesBegin() + NodePartition[k];
                ModelPart::NodeIterator NodesEnd = rModelPart.NodesBegin() + NodePartition[k + 1];

                for (ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; itNode++) {
                    array_1d<double, 3 > & OldVelocity = (itNode)->FastGetSolutionStepValue(VELOCITY, 1);
                    double& OldPressure = (itNode)->FastGetSolutionStepValue(PRESSURE, 1);

                    //predicting velocity
                    //ATTENTION::: the prediction is performed only on free nodes
                    array_1d<double, 3 > & CurrentVelocity = (itNode)->FastGetSolutionStepValue(VELOCITY);
                    double& CurrentPressure = (itNode)->FastGetSolutionStepValue(PRESSURE);

                    if ((itNode->pGetDof(VELOCITY_X))->IsFree())
                        (CurrentVelocity[0]) = OldVelocity[0];
                    if (itNode->pGetDof(VELOCITY_Y)->IsFree())
                        (CurrentVelocity[1]) = OldVelocity[1];
                    if (itNode->HasDofFor(VELOCITY_Z))
                        if (itNode->pGetDof(VELOCITY_Z)->IsFree())
                            (CurrentVelocity[2]) = OldVelocity[2];

                    if (itNode->pGetDof(PRESSURE)->IsFree())
                        CurrentPressure = OldPressure;

                    // updating time derivatives ::: please note that displacements and
                    // their time derivatives can not be consistently fixed separately
                    array_1d<double, 3 > DeltaVel;
                    noalias(DeltaVel) = CurrentVelocity - OldVelocity;
                    array_1d<double, 3 > & OldAcceleration = (itNode)->FastGetSolutionStepValue(ACCELERATION, 1);
                    array_1d<double, 3 > & CurrentAcceleration = (itNode)->FastGetSolutionStepValue(ACCELERATION);

                    SchemeType::UpdateAcceleration(CurrentAcceleration, DeltaVel, OldAcceleration);

                    // Pfem2 ALE update to eliminate convective term
                    noalias(itNode->FastGetSolutionStepValue(MESH_VELOCITY)) = itNode->FastGetSolutionStepValue(VELOCITY);
                }
            }
        }

        /*@} */
        /**@name Operations */
        /*@{ */


        /*@} */
        /**@name Access */
        /*@{ */


        /*@} */
        /**@name Inquiry */
        /*@{ */


        /*@} */
        /**@name Friends */
        /*@{ */


        /*@} */

    protected:
        /**@name Protected static Member Variables */
        /*@{ */


        /*@} */
        /**@name Protected member Variables */
        /*@{ */


        /*@} */
        /**@name Protected Operators*/
        /*@{ */


        /*@} */
        /**@name Protected Operations*/
        /*@{ */


        /*@} */
        /**@name Protected  Access */
        /*@{ */


        /*@} */
        /**@name Protected Inquiry */
        /*@{ */


        /*@} */
        /**@name Protected LifeCycle */
        /*@{ */



        /*@} */

    private:
        /**@name Static Member Variables */
        /*@{ */


        /*@} */
        /**@name Member Variables */
        /*@{ */

        /*@} */
        /**@name Private Operators*/
        /*@{ */

        /*@} */
        /**@name Private Operations*/
        /*@{ */


        /*@} */
        /**@name Private  Access */
        /*@{ */


        /*@} */
        /**@name Private Inquiry */
        /*@{ */


        /*@} */
        /**@name Un accessible methods */
        /*@{ */


        /*@} */

    }; /* Class ResidualBasedPredictorCorrectorVelocityBossakAleScheme */

    /*@} */

    /**@name Type Definitions */
    /*@{ */


    /*@} */

} /* namespace Kratos.*/

#endif /* KRATOS_RESIDUALBASED_PREDICTOR_CORRECTOR_VELOCITY_BOSSAK_ALE_SCHEME  defined */
