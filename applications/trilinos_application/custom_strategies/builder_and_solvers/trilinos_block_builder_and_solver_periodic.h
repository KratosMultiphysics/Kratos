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

#ifndef KRATOS_TRILINOS_BLOCK_BUILDER_AND_SOLVER_PERIODIC_H
#define KRATOS_TRILINOS_BLOCK_BUILDER_AND_SOLVER_PERIODIC_H

/* System includes */
#include <set>
/* #include <omp.h> */

/* External includes */
#include "boost/timer.hpp"


/* Project includes */
#include "includes/define.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "custom_strategies/builder_and_solvers/trilinos_block_builder_and_solver.h"

//trilinos includes
#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_FECrsGraph.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"

//aztec solver includes
#include "AztecOO.h"

#include "Amesos.h"
// #include "AmesosClassType.h"
#include "Epetra_LinearProblem.h"
#include "ml_MultiLevelPreconditioner.h"

namespace Kratos
{

///@addtogroup TrilinosApplication
///@{

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

template<class TSparseSpace,
         class TDenseSpace,
         class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         >
class TrilinosBlockBuilderAndSolverPeriodic
        : public TrilinosBlockBuilderAndSolver< TSparseSpace,TDenseSpace, TLinearSolver >
{
public:
    /**@name Type Definitions */
    /*@{ */
    KRATOS_CLASS_POINTER_DEFINITION( TrilinosBlockBuilderAndSolverPeriodic );


    typedef TrilinosBlockBuilderAndSolver<TSparseSpace,TDenseSpace, TLinearSolver > BaseType;

    typedef TSparseSpace SparseSpaceType;

    typedef typename BaseType::TSchemeType TSchemeType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;
    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;
    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;

    typedef typename BaseType::NodesArrayType NodesArrayType;
    typedef typename BaseType::ElementsArrayType ElementsArrayType;
    typedef typename BaseType::ConditionsArrayType ConditionsArrayType;

    typedef typename BaseType::ElementsContainerType ElementsContainerType;

    /*@} */
    /**@name Life Cycle
    */
    /*@{ */

    /** Constructor.
    */
    TrilinosBlockBuilderAndSolverPeriodic(Epetra_MpiComm& Comm,
                                          int guess_row_size,
                                          typename TLinearSolver::Pointer pNewLinearSystemSolver,
                                          const Kratos::Variable<int>& PeriodicIdVar):
        TrilinosBlockBuilderAndSolver< TSparseSpace,TDenseSpace,TLinearSolver >(Comm,guess_row_size,pNewLinearSystemSolver),
        mPeriodicIdVar(PeriodicIdVar)
    {}



    /** Destructor.
    */
    virtual ~TrilinosBlockBuilderAndSolverPeriodic() {}

    /*@} */
    /**@name Operations */
    /*@{ */


    /// Assign an Equation Id to all degrees of freedom in the system.
    /**
     * To properly set up the periodic conditions, it is assumed that all nodes on
     * a periodic boundary have the Id of their periodic "image" stored as Node->GetSolutionStepValue(mPeriodicIdVar,0).
     * @note It is assumed that the problem's partition takes the presence of periodic conditions into account, that is,
     * all processes that own periodic conditions know the nodes on both ends (either as local or ghost nodes).
     * @param rModelPart The problem's ModelPart
     */
    virtual void SetUpSystem(ModelPart &rModelPart) override
    {
        KRATOS_TRY;

        // This sort helps us preserve consistency of dof ordering across processors (as searching for a dof in
        // the nodal list of dofs will sort the list first).
        // If/when the DofList becomes a static array, this step can be skipped.
        for (ModelPart::NodeIterator itNode = rModelPart.NodesBegin(); itNode != rModelPart.NodesEnd(); ++itNode)
            itNode->GetDofs().Sort();

        unsigned int Rank = this->mrComm.MyPID();

        // Count the Dofs on this partition (on periodic node pairs, only the dofs on the node with higher Id are counted)
        int DofCount = 0;

        for (typename DofsArrayType::iterator itDof = this->mDofSet.begin(); itDof != this->mDofSet.end(); ++itDof)
            if ( (itDof->GetSolutionStepValue(PARTITION_INDEX) == Rank) && ( itDof->GetSolutionStepValue(mPeriodicIdVar) < static_cast<int>(itDof->Id())) )
                DofCount++;

        // Periodic conditions on corners are counted separately
        unsigned int ExtraDofs = SetUpEdgeDofs(rModelPart);
        DofCount += ExtraDofs;

        // Comunicate the Dof counts to set a global numbering
        int DofOffset;

        int ierr = this->mrComm.ScanSum(&DofCount,&DofOffset,1);
        if (ierr != 0)
            KRATOS_THROW_ERROR(std::runtime_error,"In TrilinosBlockBuilderAndSolverPeriodic::SetUpSystem: Found Epetra_MpiComm::ScanSum failure with error code ",ierr);

        DofOffset -= DofCount;
        this->mFirstMyId = DofOffset;

        // Assing Id to all local nodes except for the second nodes on a periodic pair
        for (typename DofsArrayType::iterator itDof = this->mDofSet.begin(); itDof != this->mDofSet.end(); ++itDof)
            if ( (itDof->GetSolutionStepValue(PARTITION_INDEX) == Rank) && ( itDof->GetSolutionStepValue(mPeriodicIdVar) < static_cast<int>(itDof->Id()) ) )
                itDof->SetEquationId(DofOffset++);
            else
                // If this rank is not responsible for assigning an EquationId to this Dof, reset the current value
                // This is necessary when dealing with changing meshes, as we risk inadvertently reusing the previous EqIdValue
                itDof->SetEquationId(0);

        // Synchronize Dof Ids across processes
        rModelPart.GetCommunicator().SynchronizeDofs();

        // Once the non-repeated EquationId are assigned and synchronized, copy the EquationId to the second node on each periodic pair
        for (ModelPart::ConditionIterator itCond = rModelPart.ConditionsBegin(); itCond != rModelPart.ConditionsEnd(); ++itCond)
        {
            ModelPart::ConditionType::GeometryType& rGeom = itCond->GetGeometry();
            if (rGeom.PointsNumber() == 2) // Regular 2-noded conditions
            {
                int Node0 = rGeom[0].Id();
                int Node0Pair = rGeom[0].FastGetSolutionStepValue(mPeriodicIdVar);

                int Node1 = rGeom[1].Id();
                int Node1Pair = rGeom[1].FastGetSolutionStepValue(mPeriodicIdVar);
                if ( Node0Pair == 0 || Node1Pair == 0 )
                    KRATOS_THROW_ERROR(std::runtime_error,"ERROR: a periodic condition has no periodic pair ids assigned. Condition Id is: ",itCond->Id());

                // If the nodes are marked as a periodic pair (this is to avoid acting on two-noded conditions that are not PeriodicCondition)
                if ( ( Node0 == Node1Pair ) && ( Node1 == Node0Pair ) )
                {
                    if ( Node0 < Node0Pair ) // If Node0 is the one with lower Id (the one that does not have an EquationId yet)
                        CopyEquationId(rGeom[1],rGeom[0]);
                    else
                        CopyEquationId(rGeom[0],rGeom[1]);
                }
            }
            else if (rGeom.PointsNumber() == 4) // Special treatment for edge nodes
            {
                for (unsigned int i = 1; i < 4; i++)
                    CopyEquationId(rGeom[0],rGeom[i]);
            }
            else if (rGeom.PointsNumber() == 8) // Special treatment for edge nodes
            {
                for (unsigned int i = 1; i < 8; i++)
                    CopyEquationId(rGeom[0],rGeom[i]);
            }
        }

        // Synchronize Dof Ids across processes a second time to transfer the Ids of periodic nodes.
        // Note that a custom communication routine is used, as the one provided by the communicator synchronizes all ghost values to the
        // value on the owner process, while here the EquationId value is given by the owner of the PeriodicCondition,
        // which is not necessarily the same as the owner of the node.
        CommunicateEquationId(rModelPart.GetCommunicator());

        // Clean temporary values to ensure that future calls to SetUpDofSet work as intended
        this->CleanEdgeDofData(rModelPart);

        // Store local and gloabal system size
        int TotalDofNum;
        ierr = this->mrComm.SumAll(&DofCount,&TotalDofNum,1);
        if (ierr != 0)
            KRATOS_THROW_ERROR(std::runtime_error,"In TrilinosBlockBuilderAndSolverPeriodic::SetUpSystem: Found Epetra_MpiComm::SumAll failure with error code ",ierr);

        this->mLocalSystemSize = DofCount;
        this->mEquationSystemSize = TotalDofNum;
        this->mLastMyId = DofOffset;

        /* some prints useful for debug
        //std::cout << Rank << ": mLocalSystemSize " << this->mLocalSystemSize << " mEquationSystemSize " << this->mEquationSystemSize << " mLastMyId " << this->mLastMyId << " mFirstMyId " << this->mFirstMyId << std::endl;
        for (typename DofsArrayType::iterator itDof = this->mDofSet.begin(); itDof != this->mDofSet.end(); ++itDof)
            //if ( (itDof->GetSolutionStepValue(PARTITION_INDEX) == Rank) && ( itDof->GetSolutionStepValue(mPeriodicIdVar) < static_cast<int>(itDof->Id())) )
            if (itDof->EquationId() == 0)
                std::cout << "Rank: " << Rank << " node id " << itDof->Id() << " eq Id " << itDof->EquationId() << " var " << itDof->GetVariable().Name() << " partition index " << itDof->GetSolutionStepValue(PARTITION_INDEX) << " flag variable " << itDof->GetSolutionStepValue(FLAG_VARIABLE) << " periodic pair index " <<  itDof->GetSolutionStepValue(mPeriodicIdVar) << std::endl;
        //std::cout << Rank << ": mLocalSystemSize " << this->mLocalSystemSize << " mEquationSystemSize " << this->mEquationSystemSize << " mLastMyId " << this->mLastMyId << " mFirstMyId " << this->mFirstMyId << std::endl;
        */

        KRATOS_CATCH("");
    }

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


    /*@} */
    /**@name Operators
    */
    /*@{ */

private:
    /**@name Static Member Variables */
    /*@{ */


    /*@} */
    /**@name Member Variables */
    /*@{ */

    const Kratos::Variable<int>& mPeriodicIdVar;

    /*@} */
    /**@name Private Operators*/
    /*@{ */


    /*@} */
    /**@name Private Operations*/
    /*@{ */

    /// Duplicate EquationIds to the second node on each periodic pair
    void CopyEquationId(ModelPart::NodeType& rOrigin,
                        ModelPart::NodeType& rDest)
    {
        for ( Node<3>::DofsContainerType::iterator itDof = rOrigin.GetDofs().begin(); itDof != rOrigin.GetDofs().end(); itDof++)
            rDest.pGetDof( itDof->GetVariable() )->SetEquationId( itDof->EquationId() );
    }

    /// Send the Equation Id of periodic nodes to the owner of the node, so it can be sync'ed across processes.
    void CommunicateEquationId(Communicator& rComm)
    {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        int destination = 0;

        Communicator::NeighbourIndicesContainerType& neighbours_indices = rComm.NeighbourIndices();

        for (unsigned int i_color = 0; i_color < neighbours_indices.size(); i_color++)
            if ((destination = neighbours_indices[i_color]) >= 0)
            {
                Communicator::NodesContainerType& r_origin_nodes = rComm.GhostMesh(i_color).Nodes();
                Communicator::NodesContainerType& r_dest_nodes = rComm.LocalMesh(i_color).Nodes();

                // Calculating send and received buffer size
                unsigned int send_buffer_size = 0;
                unsigned int receive_buffer_size = 0;

                for (Communicator::NodesContainerType::iterator i_node = r_origin_nodes.begin(); i_node != r_origin_nodes.end(); ++i_node)
                    send_buffer_size += i_node->GetDofs().size();

                for (Communicator::NodesContainerType::iterator i_node = r_dest_nodes.begin(); i_node != r_dest_nodes.end(); ++i_node)
                    receive_buffer_size += i_node->GetDofs().size();

                unsigned int position = 0;
                int* send_buffer = new int[send_buffer_size];
                int* receive_buffer = new int[receive_buffer_size];


                // Filling the buffer
                for (ModelPart::NodeIterator i_node = r_origin_nodes.begin(); i_node != r_origin_nodes.end(); ++i_node)
                    for (ModelPart::NodeType::DofsContainerType::iterator i_dof = i_node->GetDofs().begin(); i_dof != i_node->GetDofs().end(); i_dof++)
                        send_buffer[position++] = i_dof->EquationId();


                MPI_Status status;

                if (position > send_buffer_size)
                    std::cout << rank << " Error in estimating send buffer size...." << std::endl;


                int send_tag = i_color;
                int receive_tag = i_color;

                MPI_Sendrecv(send_buffer, send_buffer_size, MPI_INT, destination, send_tag, receive_buffer, receive_buffer_size, MPI_INT, destination, receive_tag,
                             MPI_COMM_WORLD, &status);

                // Updating nodes
                position = 0;
                for (ModelPart::NodeIterator i_node = r_dest_nodes.begin();
                     i_node != r_dest_nodes.end(); i_node++)
                    for (ModelPart::NodeType::DofsContainerType::iterator i_dof = i_node->GetDofs().begin(); i_dof != i_node->GetDofs().end(); i_dof++)
                    {
                        unsigned int NewId = static_cast<unsigned int>(receive_buffer[position++]);
                        if (NewId > i_dof->EquationId()) // Note: in a general case, only one rank will have assinged an EquationId, the others will send 0s
                            i_dof->SetEquationId(NewId);
                    }

                if (position > receive_buffer_size)
                    std::cout << rank << " Error in estimating receive buffer size...." << std::endl;

                delete [] send_buffer;
                delete [] receive_buffer;
            }

        // Once we are sure that the owner process of all nodes knows its equation id, replicate its value to all ghost images
        rComm.SynchronizeDofs();
    }

    /// Additional operations required to assign a correct EquationId to periodic nodes corresponding to an edge (four images of the Dof exist in the domain).
    /** Note that the EquationId is not assigned by this function (it will be assigned in SetUpSystem)
      * @return Number of additional degrees of freedom in this partition, due to the presence of edge nodes
      */
    unsigned int SetUpEdgeDofs(ModelPart& rModelPart)
    {
        int LocalNodesNum = rModelPart.GetCommunicator().LocalMesh().Nodes().size();
        int GlobalNodesNum = 0;
        this->mrComm.SumAll(&LocalNodesNum,&GlobalNodesNum,1);

        int NumProcs = this->mrComm.NumProc();
        int* ExtraDofs = new int[NumProcs];
        for (int i = 0; i < NumProcs; i++) ExtraDofs[i] = 0;

        Condition::DofsVectorType DofList;
        ProcessInfo& rProcessInfo = rModelPart.GetProcessInfo();
        for (ModelPart::ConditionIterator itCond = rModelPart.ConditionsBegin(); itCond != rModelPart.ConditionsEnd(); ++itCond)
        {
            Condition::GeometryType& rGeom = itCond->GetGeometry();
            const unsigned int NumNodes = rGeom.PointsNumber();
            // mPeriodicIdVar > GlobalNodesNum is an imposible value for regular periodic nodes, which we use to signal edges
            if ( ( NumNodes == 4 || NumNodes == 8 ) && rGeom[0].FastGetSolutionStepValue(mPeriodicIdVar) > GlobalNodesNum )
            {
                unsigned int FirstNode = rGeom[0].Id();
                // KLUDGE: the following call should go through the scheme!!
                itCond->GetDofList(DofList,rProcessInfo);
                for(typename Condition::DofsVectorType::iterator iDof = DofList.begin() ; iDof != DofList.end() ; ++iDof)
                    if ( (*iDof)->Id() == FirstNode)
                        ExtraDofs[ (unsigned int)( (*iDof)->GetSolutionStepValue(PARTITION_INDEX) ) ]++;

                rGeom[0].GetValue(mPeriodicIdVar) = rGeom[0].FastGetSolutionStepValue(mPeriodicIdVar);
            }
        }

        int* TotalExtraDofs = new int[NumProcs];
        for (int i = 0; i < NumProcs; i++) TotalExtraDofs[i] = 0;
        
        this->mrComm.SumAll(ExtraDofs,TotalExtraDofs,NumProcs);
        
        // free memory
        int LocalExtraDofs = TotalExtraDofs[this->mrComm.MyPID()];
        delete [] ExtraDofs;
        delete [] TotalExtraDofs;

        // Prepare edge dofs so only one of the duplicates gets an equation Id
        rModelPart.GetCommunicator().AssembleNonHistoricalData(mPeriodicIdVar);
        for ( ModelPart::NodeIterator iNode = rModelPart.NodesBegin(); iNode != rModelPart.NodesEnd(); iNode++)
            if ( iNode->GetValue(mPeriodicIdVar) != 0)
                iNode->FastGetSolutionStepValue(mPeriodicIdVar) = 0;

        return LocalExtraDofs;
    }

    /// Clean temporary values written in previous calls to SetUpEdgeDofs
    /** This is required for all calls after the first if the DofSet and the system need rebuilding
      */
    void CleanEdgeDofData(ModelPart& rModelPart)
    {
        for ( ModelPart::NodeIterator iNode = rModelPart.NodesBegin(); iNode != rModelPart.NodesEnd(); iNode++)
            if ( iNode->GetValue(mPeriodicIdVar) != 0)
            {
                iNode->FastGetSolutionStepValue(mPeriodicIdVar) = iNode->GetValue(mPeriodicIdVar);
                iNode->GetValue(mPeriodicIdVar) = 0;
            }
    }
    

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

}; /* Class TrilinosBlockBuilderAndSolverPeriodic */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

///@} //addtogroup

}  /* namespace Kratos.*/


#endif // KRATOS_TRILINOS_BLOCK_BUILDER_AND_SOLVER_PERIODIC_H
