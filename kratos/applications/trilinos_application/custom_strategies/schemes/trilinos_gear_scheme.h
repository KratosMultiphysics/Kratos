/*
==============================================================================
KratosTrilinosApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
 */

#if !defined(KRATOS_TRILINOS_GEAR_SCHEME )
#define  KRATOS_TRILINOS_GEAR_SCHEME


/* System includes */


/* External includes */
#include "boost/smart_ptr.hpp"


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
        mImporterIsInitialized(false),
        mrPeriodicIdVar(Kratos::Variable<int>::StaticObject())
    {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
        if (rank == 0)
            std::cout << "Using the Trilinos BDF2 time scheme" << std::endl;
    }

    TrilinosGearScheme(Process::Pointer pTurbulenceModel):
        GearScheme<TSparseSpace,TDenseSpace>(pTurbulenceModel),
        mImporterIsInitialized(false),
        mrPeriodicIdVar(Kratos::Variable<int>::StaticObject())
    {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
        if (rank == 0)
            std::cout << "Using the Trilinos BDF2 time scheme (with turbulence model)" << std::endl;
    }

    TrilinosGearScheme(const Variable<int>& rPeriodicIdVar):
        GearScheme<TSparseSpace,TDenseSpace>(),
        mImporterIsInitialized(false),
        mrPeriodicIdVar(rPeriodicIdVar)
    {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
        if (rank == 0)
            std::cout << "Using the Trilinos BDF2 time scheme (with periodic conditions)" << std::endl;
    }
    /** Destructor.
     */
    virtual ~TrilinosGearScheme()
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


    /// Update solution step data with the result of the last linear system solution
    /**
     * @param r_model_part Problem ModelPart.
     * @param rDofSet Array of degreees of freedom of the system.
     * @param A System matrix (unused).
     * @param Dx Vector of nodal unknowns (increments to be added to each unknown value).
     * @param b System right hand side vector (unused).
     */
//    virtual void BasicUpdateOperations(ModelPart& r_model_part,
//                                       DofsArrayType& rDofSet,
//                                       TSystemMatrixType& A,
//                                       TSystemVectorType& Dx,
//                                       TSystemVectorType& b)
    virtual void UpdateDofs(DofsArrayType& rDofSet,
                            TSystemVectorType& Dx)
    {
        KRATOS_TRY;

        if (!DofImporterIsInitialized())
            this->InitializeDofImporter(rDofSet,Dx);

        int system_size = TSparseSpace::Size(Dx);

        //defining a temporary vector to gather all of the values needed
        Epetra_Vector temp( mpDofImporter->TargetMap() );

        //importing in the new temp vector the values
        int ierr = temp.Import(Dx,*mpDofImporter,Insert) ;
        if(ierr != 0) KRATOS_THROW_ERROR(std::logic_error,"Epetra failure found","");

        double* temp_values;
        temp.ExtractView( &temp_values );

        Dx.Comm().Barrier();

        //performing the update
        for (typename DofsArrayType::iterator itDof = rDofSet.begin(); itDof != rDofSet.end(); itDof++)
        {
            int global_id = itDof->EquationId();
            if(global_id < system_size)
            {
                double aaa = temp[mpDofImporter->TargetMap().LID(global_id)];
                itDof->GetSolutionStepValue() += aaa;
            }
        }

        KRATOS_CATCH("");
    }

    bool DofImporterIsInitialized()
    {
        return mImporterIsInitialized;
    }

    virtual void Clear()
    {
        mpDofImporter.reset();
        mImporterIsInitialized = false;
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

    virtual void InitializeDofImporter(DofsArrayType& rDofSet,
                                       TSystemVectorType& Dx)
    {
        int system_size = TSparseSpace::Size(Dx);
        int number_of_dofs = rDofSet.size();
        std::vector< int > index_array(number_of_dofs);

        //filling the array with the global ids
        int counter = 0;
        for(typename DofsArrayType::iterator i_dof = rDofSet.begin() ; i_dof != rDofSet.end() ; ++i_dof)
        {
            int id = i_dof->EquationId();
            if( id < system_size )
            {
                index_array[counter] = id;
                counter += 1;
            }
        }

        std::sort(index_array.begin(),index_array.end());
        std::vector<int>::iterator NewEnd = std::unique(index_array.begin(),index_array.end());
        index_array.resize(NewEnd-index_array.begin());

        int check_size = -1;
        int tot_update_dofs = index_array.size();
        Dx.Comm().SumAll(&tot_update_dofs,&check_size,1);
        if ( (check_size < system_size) &&  (Dx.Comm().MyPID() == 0) )
        {
            std::stringstream Msg;
            Msg << "Dof count is not correct. There are less dofs then expected." << std::endl;
            Msg << "Expected number of active dofs = " << system_size << " dofs found = " << check_size << std::endl;
            KRATOS_THROW_ERROR(std::runtime_error,Msg.str(),"")
        }

        //defining a map as needed
        Epetra_Map dof_update_map(-1,index_array.size(), &(*(index_array.begin())),0,Dx.Comm() );

        //defining the importer class
        boost::shared_ptr<Epetra_Import> pDofImporter( new Epetra_Import(dof_update_map,Dx.Map()) );
        mpDofImporter.swap(pDofImporter);

        mImporterIsInitialized = true;
    }

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

    bool mImporterIsInitialized;

    boost::shared_ptr<Epetra_Import> mpDofImporter;

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


