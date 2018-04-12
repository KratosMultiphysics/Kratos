//  KRATOS  _____     _ _ _
//         |_   _| __(_) (_)_ __   ___  ___
//           | || '__| | | | '_ \ / _ \/ __|
//           | || |  | | | | | | | (_) \__
//           |_||_|  |_|_|_|_| |_|\___/|___/ APPLICATION
//
//  License:             BSD License
//                                       Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//


#if !defined(KRATOS_TRILINOS_GEAR_SCHEME )
#define  KRATOS_TRILINOS_GEAR_SCHEME


/* System includes */


/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "../../../FluidDynamicsApplication/custom_strategies/strategies/gear_scheme.h"
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

/// Trilinos version of GearScheme.
template<class TSparseSpace,
         class TDenseSpace //= DenseSpace<double>
         >
class TrilinosGearScheme
    : public GearScheme<TSparseSpace, TDenseSpace>
{
public:
    /**@name Type Definitions */
    /*@{ */

    KRATOS_CLASS_POINTER_DEFINITION(TrilinosGearScheme );

    typedef GearScheme<TSparseSpace, TDenseSpace> BaseType;

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

    TrilinosGearScheme():
        GearScheme<TSparseSpace,TDenseSpace>(),
        mrPeriodicIdVar(Kratos::Variable<int>::StaticObject())
    {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
        if (rank == 0)
            std::cout << "Using the Trilinos BDF2 time scheme" << std::endl;
    }

    TrilinosGearScheme(Process::Pointer pTurbulenceModel):
        GearScheme<TSparseSpace,TDenseSpace>(pTurbulenceModel),
        mrPeriodicIdVar(Kratos::Variable<int>::StaticObject())
    {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
        if (rank == 0)
            std::cout << "Using the Trilinos BDF2 time scheme (with turbulence model)" << std::endl;
    }

    TrilinosGearScheme(const Variable<int>& rPeriodicIdVar):
        GearScheme<TSparseSpace,TDenseSpace>(),
        mrPeriodicIdVar(rPeriodicIdVar)
    {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
        if (rank == 0)
            std::cout << "Using the Trilinos BDF2 time scheme (with periodic conditions)" << std::endl;
    }
    /** Destructor.
     */
    ~TrilinosGearScheme() override
    {
    }


    /*@} */
    /**@name Operators
     */
    /*@{ */

    int Check(ModelPart& rModelPart) override
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
    void FinalizeNonLinIteration(ModelPart &rModelPart, TSystemMatrixType &A, TSystemVectorType &Dx, TSystemVectorType &b) override
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

    /** On periodic boundaries, the nodal area and the values to project need to take into account contributions from elements on
     * both sides of the boundary. This is done using the conditions and the non-historical nodal data containers as follows:\n
     * 1- The partition that owns the PeriodicCondition adds the values on both nodes to their non-historical containers.\n
     * 2- The non-historical containers are added across processes, transmiting the right value from the condition owner to all partitions.\n
     * 3- The value on all periodic nodes is replaced by the one received in step 2.
     */
    void PeriodicConditionProjectionCorrection(ModelPart& rModelPart)
    {
        if (mrPeriodicIdVar.Key() != 0)
        {
            for (typename ModelPart::ConditionIterator itCond = rModelPart.ConditionsBegin(); itCond != rModelPart.ConditionsEnd(); itCond++ )
            {
                ModelPart::ConditionType::GeometryType& rGeom = itCond->GetGeometry();
                if (rGeom.PointsNumber() == 2)
                {
                    Node<3>& rNode0 = rGeom[0];
                    int Node0Pair = rNode0.FastGetSolutionStepValue(mrPeriodicIdVar);

                    Node<3>& rNode1 = rGeom[1];
                    int Node1Pair = rNode1.FastGetSolutionStepValue(mrPeriodicIdVar);

                    // If the nodes are marked as a periodic pair (this is to avoid acting on two-noded conditions that are not PeriodicCondition)
                    if ( ( static_cast<int>(rNode0.Id()) == Node1Pair ) && (static_cast<int>(rNode1.Id()) == Node0Pair ) )
                    {
                        double NodalArea = rNode0.FastGetSolutionStepValue(NODAL_AREA) + rNode1.FastGetSolutionStepValue(NODAL_AREA);
                        array_1d<double,3> AdvProj = rNode0.FastGetSolutionStepValue(ADVPROJ) + rNode1.FastGetSolutionStepValue(ADVPROJ);
                        double DivProj = rNode0.FastGetSolutionStepValue(DIVPROJ) + rNode1.FastGetSolutionStepValue(DIVPROJ);

                        rNode0.GetValue(NODAL_AREA) = NodalArea;
                        rNode0.GetValue(ADVPROJ) = AdvProj;
                        rNode0.GetValue(DIVPROJ) = DivProj;

                        rNode1.GetValue(NODAL_AREA) = NodalArea;
                        rNode1.GetValue(ADVPROJ) = AdvProj;
                        rNode1.GetValue(DIVPROJ) = DivProj;
                    }
                }
            }

            rModelPart.GetCommunicator().AssembleNonHistoricalData(NODAL_AREA);
            rModelPart.GetCommunicator().AssembleNonHistoricalData(ADVPROJ);
            rModelPart.GetCommunicator().AssembleNonHistoricalData(DIVPROJ);

            for (typename ModelPart::NodeIterator itNode = rModelPart.NodesBegin(); itNode != rModelPart.NodesEnd(); itNode++)
            {
                if (itNode->GetValue(NODAL_AREA) != 0.0)
                {
                    itNode->FastGetSolutionStepValue(NODAL_AREA) = itNode->GetValue(NODAL_AREA);
                    itNode->FastGetSolutionStepValue(ADVPROJ) = itNode->GetValue(ADVPROJ);
                    itNode->FastGetSolutionStepValue(DIVPROJ) = itNode->GetValue(DIVPROJ);

                    // reset for next iteration
                    itNode->GetValue(NODAL_AREA) = 0.0;
                    itNode->GetValue(ADVPROJ) = array_1d<double,3>(3,0.0);
                    itNode->GetValue(DIVPROJ) = 0.0;
                }
            }
        }
    }

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

    const Kratos::Variable<int>& mrPeriodicIdVar;

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

#endif /* KRATOS_TRILINOS_GEAR_SCHEME  defined */


