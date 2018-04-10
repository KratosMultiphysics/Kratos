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
//                   Riccardo Rossi
//

#if !defined(KRATOS_TRILINOS_PREDICTOR_CORRECTOR_VELOCITY_BOSSAK_SCHEME_TURBULENT )
#define  KRATOS_TRILINOS_PREDICTOR_CORRECTOR_VELOCITY_BOSSAK_SCHEME_TURBULENT


/* System includes */


/* External includes */
#include "boost/smart_ptr.hpp"


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "../../../FluidDynamicsApplication/custom_strategies/strategies/residualbased_predictorcorrector_velocity_bossak_scheme_turbulent.h"
#include "includes/variables.h"
#include "containers/array_1d.h"

namespace Kratos
{




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

/// Trilinos version of ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent.
template<class TSparseSpace,
         class TDenseSpace //= DenseSpace<double>
         >
class TrilinosPredictorCorrectorVelocityBossakSchemeTurbulent
    : public ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent<TSparseSpace, TDenseSpace>
{
public:
    /**@name Type Definitions */
    /*@{ */

    KRATOS_CLASS_POINTER_DEFINITION(TrilinosPredictorCorrectorVelocityBossakSchemeTurbulent );

    typedef ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent<TSparseSpace, TDenseSpace> BaseType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename Element::DofsVectorType DofsVectorType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;


    /*@} */
    /**@name Life Cycle
     */
    /*@{ */

    /** Constructor.
     */

    TrilinosPredictorCorrectorVelocityBossakSchemeTurbulent(double NewAlphaBossak, double MoveMeshStrategy, unsigned int DomainSize):
        ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent<TSparseSpace,TDenseSpace>(NewAlphaBossak, MoveMeshStrategy, DomainSize)
    {
        std::cout << "using the TRILINOS velocity Bossak Time Integration Scheme (with turbulence model)" << std::endl;
    }

    TrilinosPredictorCorrectorVelocityBossakSchemeTurbulent(double NewAlphaBossak, double MoveMeshStrategy, unsigned int DomainSize, Process::Pointer pTurbulenceModel):
        ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent<TSparseSpace,TDenseSpace>(NewAlphaBossak, MoveMeshStrategy, DomainSize, pTurbulenceModel)
    {

        std::cout << "using the TRILINOS velocity Bossak Time Integration Scheme (with turbulence model)" << std::endl;
    }

    TrilinosPredictorCorrectorVelocityBossakSchemeTurbulent(double NewAlphaBossak, double MoveMeshStrategy, unsigned int DomainSize, const Variable<int>& rPeriodicIdVar):
        ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent<TSparseSpace,TDenseSpace>(NewAlphaBossak, MoveMeshStrategy, DomainSize, rPeriodicIdVar)
    {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
        if (rank == 0)
            std::cout << "using the TRILINOS velocity Bossak Time Integration Scheme (with periodic conditions)" << std::endl;
    }
    /** Destructor.
     */
    virtual ~TrilinosPredictorCorrectorVelocityBossakSchemeTurbulent()
    {
    }


    /*@} */
    /**@name Operators
     */
    /*@{ */

    virtual int Check(ModelPart& rModelPart)
    {
        KRATOS_TRY

        int ErrorCode = BaseType::Check(rModelPart);
        if (ErrorCode != 0) return ErrorCode;

        // Check buffer size
        if (rModelPart.GetBufferSize() < 2)
            KRATOS_THROW_ERROR(std::logic_error, "GearScheme error: Insufficient buffer size for Bossak scheme, should be at least 2, got ",rModelPart.GetBufferSize());

        // Check that all required variables were registered
        if(DELTA_TIME.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"DELTA_TIME Key is 0. Check if all applications were correctly registered.","");
        if(OSS_SWITCH.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"OSS_SWITCH Key is 0. Check if all applications were correctly registered.","");

        if(DISPLACEMENT.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"DISPLACEMENT Key is 0. Check if all applications were correctly registered.","");
        if(VELOCITY.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"VELOCITY Key is 0. Check if all applications were correctly registered.","");
        if(MESH_VELOCITY.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"MESH_VELOCITY Key is 0. Check if all applications were correctly registered.","");
        if(ACCELERATION.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"ACCELERATION Key is 0. Check if all applications were correctly registered.","");

        // Checks on process info
//            const ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();

//            if(rCurrentProcessInfo.Has(DELTA_TIME) == 0)
//                KRATOS_THROW_ERROR(std::invalid_argument,"ProcessInfo does not contain a value for DELTA_TIME","");

        return 0;
        KRATOS_CATCH("");
    }


    /// Calculate OSS projections, properly taking into account periodic boundaries.
    virtual void FinalizeNonLinIteration(ModelPart &rModelPart, TSystemMatrixType &A, TSystemVectorType &Dx, TSystemVectorType &b)
    {
        ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

        //if orthogonal subscales are computed
        if (CurrentProcessInfo[OSS_SWITCH] == 1.0) {
            if (rModelPart.GetCommunicator().MyPID() == 0)
                std::cout << "Computing OSS projections" << std::endl;
            for (typename ModelPart::NodesContainerType::iterator ind = rModelPart.NodesBegin(); ind != rModelPart.NodesEnd(); ind++) {

                noalias(ind->FastGetSolutionStepValue(ADVPROJ)) = ZeroVector(3);

                ind->FastGetSolutionStepValue(DIVPROJ) = 0.0;

                ind->FastGetSolutionStepValue(NODAL_AREA) = 0.0;


            }//end of loop over nodes

            //loop on nodes to compute ADVPROJ   CONVPROJ NODALAREA
            array_1d<double, 3 > output;


            for (typename ModelPart::ElementsContainerType::iterator elem = rModelPart.ElementsBegin(); elem != rModelPart.ElementsEnd(); elem++)
            {
                elem->Calculate(ADVPROJ, output, CurrentProcessInfo);
            }

            rModelPart.GetCommunicator().AssembleCurrentData(NODAL_AREA);
            rModelPart.GetCommunicator().AssembleCurrentData(DIVPROJ);
            rModelPart.GetCommunicator().AssembleCurrentData(ADVPROJ);


            // Correction for periodic conditions
            this->PeriodicConditionProjectionCorrection(rModelPart);


            for (typename ModelPart::NodesContainerType::iterator ind = rModelPart.NodesBegin(); ind != rModelPart.NodesEnd(); ind++)
            {
                if (ind->FastGetSolutionStepValue(NODAL_AREA) == 0.0)
                {
                    ind->FastGetSolutionStepValue(NODAL_AREA) = 1.0;
                    //KRATOS_WATCH("*********ATTENTION: NODAL AREA IS ZERRROOOO************");
                }
                const double Area = ind->FastGetSolutionStepValue(NODAL_AREA);
                ind->FastGetSolutionStepValue(ADVPROJ) /= Area;
                ind->FastGetSolutionStepValue(DIVPROJ) /= Area;
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

}; /* Class Scheme */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

} /* namespace Kratos.*/

#endif /* KRATOS_TRILINOS_PREDICTOR_CORRECTOR_VELOCITY_BOSSAK_SCHEME_TURBULENT  defined */



