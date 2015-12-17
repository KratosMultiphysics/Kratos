/*
Kratos Multi-Physics

Copyright (c) 2015, Pooyan Dadvand, Riccardo Rossi, CIMNE (International Center for Numerical Methods in Engineering)
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

    -	Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
    -	Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer
        in the documentation and/or other materials provided with the distribution.
    -	All advertising materials mentioning features or use of this software must display the following acknowledgement:
            This product includes Kratos Multi-Physics technology.
    -	Neither the name of the CIMNE nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED ANDON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: jcotela $
//   Date:                $Date: 2015-12-14
//
//


#if !defined(KRATOS_TRILINOS_VARIATIONAL_DISTANCE_CALCULATION_PROCESS_INCLUDED )
#define  KRATOS_TRILINOS_VARIATIONAL_DISTANCE_CALCULATION_PROCESS_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes
#include "mpi.h"
#include "Epetra_MpiComm.h"

// Project includes
#include "includes/communicator.h"
#include "processes/variational_distance_calculation_process.h"
#include "custom_strategies/builder_and_solvers/trilinos_block_builder_and_solver.h"
#include "custom_strategies/schemes/trilinos_residualbased_incrementalupdate_static_scheme.h"


namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{


///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{




/// Short class definition.
/**takes a model part full of SIMPLICIAL ELEMENTS (triangles and tetras) and recomputes a signed distance function
mantaining as much as possible the position of the zero of the function prior to the call.

This is achieved by minimizing the function  ( 1 - norm( gradient( distance ) )**2
with the restriction that "distance" is a finite elment function
*/

template< unsigned int TDim,
          class TSparseSpace,
          class TDenseSpace,
          class TLinearSolver >
class TrilinosVariationalDistanceCalculationProcess
    : public VariationalDistanceCalculationProcess<TDim,TSparseSpace,TDenseSpace,TLinearSolver>
{
public:

    ///@name Type Definitions
    ///@{

    typedef VariationalDistanceCalculationProcess<TDim,TSparseSpace,TDenseSpace,TLinearSolver> BaseType;


    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of VariationalDistanceCalculationProcess
    KRATOS_CLASS_POINTER_DEFINITION(TrilinosVariationalDistanceCalculationProcess);




    ///@}
    ///@name Life Cycle
    ///@{


    TrilinosVariationalDistanceCalculationProcess(Epetra_MpiComm& rComm,
                                                  ModelPart& base_model_part,
                                                  typename TLinearSolver::Pointer plinear_solver,
                                                  unsigned int max_iterations = 10):
        VariationalDistanceCalculationProcess<TDim,TSparseSpace,TDenseSpace,TLinearSolver>(base_model_part,max_iterations),
        mrComm(rComm)
    {
        KRATOS_TRY

        //check that there is at least one element and node in the model
        int NNode = base_model_part.Nodes().size();
        int NElem = base_model_part.Elements().size();

        if(NNode > 0 && base_model_part.NodesBegin()->SolutionStepsDataHas(DISTANCE) == false )
            KRATOS_THROW_ERROR(std::invalid_argument,"missing DISTANCE variable on solution step data","");

        if(NElem > 0)
        {
            if(TDim == 2)
            {
                if(base_model_part.ElementsBegin()->GetGeometry().GetGeometryFamily() != GeometryData::Kratos_Triangle)
                    KRATOS_THROW_ERROR(std::logic_error, "In 2D the element type is expected to be a triangle","");
            }
            else if(TDim == 3)
            {
                if(base_model_part.ElementsBegin()->GetGeometry().GetGeometryFamily() != GeometryData::Kratos_Tetrahedra)
                    KRATOS_THROW_ERROR(std::logic_error, "In 3D the element type is expected to be a tetrahedra","");
            }
        }


        base_model_part.GetCommunicator().SumAll(NNode);
        base_model_part.GetCommunicator().SumAll(NElem);

        if( base_model_part.GetCommunicator().MyPID() == 0 )
        {
            if( NNode == 0 ) KRATOS_THROW_ERROR(std::logic_error, "the model has no Nodes","");
            if( NElem == 0 ) KRATOS_THROW_ERROR(std::logic_error, "the model has no Elements","");
        }

        //generate an auxilary model part and populate it by elements of type DistanceCalculationElementSimplex
        this->ReGenerateDistanceModelPart(base_model_part);

        //generate a linear strategy

        // Scheme
        typename BaseType::SchemeType::Pointer pscheme = typename BaseType::SchemeType::Pointer( new TrilinosResidualBasedIncrementalUpdateStaticScheme< TSparseSpace,TDenseSpace >() );

        // Builder and Solver
        int RowSizeGuess = (TDim == 2 ? 15 : 40);
        typedef typename BuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>::Pointer BuilderSolverTypePointer;
        BuilderSolverTypePointer pBuilderSolver = BuilderSolverTypePointer(new TrilinosBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver >(mrComm,RowSizeGuess,plinear_solver) );

        // Solution strategy
        bool CalculateReactions = false;
        bool ReformDofAtEachIteration = false;
        bool CalculateNormDxFlag = false;

        this->mp_solving_strategy = typename BaseType::SolvingStrategyType::Pointer( new ResidualBasedLinearStrategy<TSparseSpace,TDenseSpace,TLinearSolver>(*(this->mp_distance_model_part),pscheme,plinear_solver,pBuilderSolver,CalculateReactions,ReformDofAtEachIteration,CalculateNormDxFlag) );

        //TODO: check flag DO_EXPENSIVE_CHECKS
        this->mp_solving_strategy->Check();

        KRATOS_CATCH("")
    }

    /// Destructor.
    virtual ~TrilinosVariationalDistanceCalculationProcess()
    {
    }


    ///@}
    ///@name Operators
    ///@{



    ///@}
    ///@name Operations
    ///@{


    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    virtual void ReGenerateDistanceModelPart(ModelPart& base_model_part)
    {
        KRATOS_TRY

        //generate
        ModelPart::Pointer pAuxModelPart = ModelPart::Pointer( new ModelPart("DistancePart",1) );

        ModelPart::Pointer& p_distance_model_part = this->mp_distance_model_part;
        p_distance_model_part.swap(pAuxModelPart);

        p_distance_model_part->Nodes().clear();
        p_distance_model_part->Conditions().clear();
        p_distance_model_part->Elements().clear();

        p_distance_model_part->SetProcessInfo(  base_model_part.pGetProcessInfo() );
        p_distance_model_part->SetBufferSize(base_model_part.GetBufferSize());
        p_distance_model_part->SetProperties(base_model_part.pProperties());
        p_distance_model_part->Tables() = base_model_part.Tables();

        //assigning the nodes to the new model part
        p_distance_model_part->Nodes() = base_model_part.Nodes();

        //ensure that the nodes have distance as a DOF
        for (ModelPart::NodesContainerType::iterator iii = base_model_part.NodesBegin(); iii != base_model_part.NodesEnd(); iii++)
        {
            iii->AddDof(DISTANCE);
        }

        // For MPI: copy communication data
        Communicator& rRefComm = base_model_part.GetCommunicator();
        Communicator::Pointer pNewComm = rRefComm.Create();

        pNewComm->SetNumberOfColors( rRefComm.GetNumberOfColors() ) ;
        pNewComm->NeighbourIndices() = rRefComm.NeighbourIndices();
        pNewComm->LocalMesh().SetNodes( rRefComm.LocalMesh().pNodes() );
        pNewComm->InterfaceMesh().SetNodes( rRefComm.InterfaceMesh().pNodes() );
        pNewComm->GhostMesh().SetNodes( rRefComm.GhostMesh().pNodes() );
        for (unsigned int i = 0; i < rRefComm.GetNumberOfColors(); i++)
        {
            pNewComm->pInterfaceMesh(i)->SetNodes( rRefComm.pInterfaceMesh(i)->pNodes() );
            pNewComm->pLocalMesh(i)->SetNodes( rRefComm.pLocalMesh(i)->pNodes() );
            pNewComm->pGhostMesh(i)->SetNodes( rRefComm.pGhostMesh(i)->pNodes() );
        }

        p_distance_model_part->SetCommunicator(pNewComm);

        //generating the elements
        p_distance_model_part->Elements().reserve(base_model_part.Elements().size());
        for (ModelPart::ElementsContainerType::iterator iii = base_model_part.ElementsBegin(); iii != base_model_part.ElementsEnd(); iii++)
        {
            Properties::Pointer properties = iii->pGetProperties();
            Element::Pointer p_element = Element::Pointer(new DistanceCalculationElementSimplex<TDim>(
                                             iii->Id(),
                                             iii->pGetGeometry(),
                                             iii->pGetProperties() ) );

            //assign EXACTLY THE SAME GEOMETRY, so that memory is saved!!
            p_element->pGetGeometry() = iii->pGetGeometry();

            p_distance_model_part->Elements().push_back(p_element);
            pNewComm->LocalMesh().Elements().push_back(p_element);
        }


        //using the conditions to mark the boundary with the flag boundary
        //note that we DO NOT add the conditions to the model part
        for (ModelPart::NodesContainerType::iterator iii = p_distance_model_part->NodesBegin(); iii != p_distance_model_part->NodesEnd(); iii++)
        {
            iii->Set(BOUNDARY,false);
        }
        for (ModelPart::ConditionsContainerType::iterator iii = base_model_part.ConditionsBegin(); iii != base_model_part.ConditionsEnd(); iii++)
        {
            Geometry< Node<3> >& geom = iii->GetGeometry();
            for(unsigned int i=0; i<geom.size(); i++) geom[i].Set(BOUNDARY,true);
        }

        // Communicate BOUNDARY status to all partitions
        this->CommunicateBoundaryFlagToOwner(*p_distance_model_part);
        this->CommunicateBoundaryFlagFromOwner(*p_distance_model_part);

        this->mdistance_part_is_initialized = true;

        //KRATOS_WATCH(this->mp_distance_model_part)

        KRATOS_CATCH("")
    }


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    Epetra_MpiComm& mrComm;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void CommunicateBoundaryFlagToOwner(ModelPart& rModelPart)
    {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        int destination = 0;

        Communicator::NeighbourIndicesContainerType& neighbours_indices = rModelPart.GetCommunicator().NeighbourIndices();

        std::vector<int*> recv_buffers(neighbours_indices.size());
        std::vector<int>  recv_sizes(neighbours_indices.size());

        for (unsigned int i_color = 0; i_color < neighbours_indices.size(); i_color++)
            if ((destination = neighbours_indices[i_color]) >= 0)
            {
                Communicator::NodesContainerType& r_local_nodes = rModelPart.GetCommunicator().InterfaceMesh(i_color).Nodes();
                Communicator::NodesContainerType& r_ghost_nodes = rModelPart.GetCommunicator().InterfaceMesh(i_color).Nodes();

                unsigned int send_size = r_ghost_nodes.size();
                unsigned int recv_size = r_local_nodes.size();

                if ( (send_size == 0) && (recv_size == 0) )
                    continue; // Nothing to transfer

                int* send_buffer = new int[send_size];
                recv_buffers[i_color] = new int[recv_size];
                recv_sizes[i_color] = recv_size;

                // Fill the send buffer
                unsigned int i = 0;
                for (ModelPart::NodeIterator i_node = r_ghost_nodes.begin(); i_node != r_ghost_nodes.end(); ++i_node)
                {
                    send_buffer[i++] = i_node->Is(BOUNDARY); // bool to int! (should be safe)
                }

                if (i > send_size)
                    std::cout << rank << " Error in estimating send buffer size...." << std::endl;

                MPI_Status status;

                int send_tag = i_color;
                int receive_tag = i_color;

                MPI_Sendrecv(send_buffer, send_size, MPI_INT, destination, send_tag,
                             recv_buffers[i_color], recv_size, MPI_INT, destination, receive_tag,
                             MPI_COMM_WORLD, &status);

                delete [] send_buffer;
            }

        // Write in nodes
        for (unsigned int i_color = 0; i_color < neighbours_indices.size(); i_color++)
            if ((destination = neighbours_indices[i_color]) >= 0)
            {
                Communicator::NodesContainerType& r_local_nodes = rModelPart.GetCommunicator().InterfaceMesh(i_color).Nodes();

                int* recv_buffer = recv_buffers[i_color];

                unsigned int i = 0;
                for (ModelPart::NodeIterator i_node = r_local_nodes.begin(); i_node != r_local_nodes.end(); ++i_node)
                {
                    i_node->Set(BOUNDARY, bool(recv_buffer[i++]) || i_node->Is(BOUNDARY) ); // OR with received value
                }

                if (i > recv_sizes[i_color])
                    std::cout << rank << " Error in estimating receive buffer size...." << std::endl;

                delete[] recv_buffer;
            }
    }

    void CommunicateBoundaryFlagFromOwner(ModelPart& rModelPart)
    {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        int destination = 0;

        Communicator::NeighbourIndicesContainerType& neighbours_indices = rModelPart.GetCommunicator().NeighbourIndices();

        std::vector<int*> recv_buffers(neighbours_indices.size());
        std::vector<int>  recv_sizes(neighbours_indices.size());

        for (unsigned int i_color = 0; i_color < neighbours_indices.size(); i_color++)
            if ((destination = neighbours_indices[i_color]) >= 0)
            {
                Communicator::NodesContainerType& r_local_nodes = rModelPart.GetCommunicator().InterfaceMesh(i_color).Nodes();
                Communicator::NodesContainerType& r_ghost_nodes = rModelPart.GetCommunicator().InterfaceMesh(i_color).Nodes();

                unsigned int send_size = r_local_nodes.size();
                unsigned int recv_size = r_ghost_nodes.size();

                if ( (send_size == 0) && (recv_size == 0) )
                    continue; // Nothing to transfer

                int* send_buffer = new int[send_size];
                recv_buffers[i_color] = new int[recv_size];
                recv_sizes[i_color] = recv_size;

                // Fill the send buffer
                unsigned int i = 0;
                for (ModelPart::NodeIterator i_node = r_local_nodes.begin(); i_node != r_local_nodes.end(); ++i_node)
                {
                    send_buffer[i++] = i_node->Is(BOUNDARY); // bool to int! (should be safe)
                }

                if (i > send_size)
                    std::cout << rank << " Error in estimating send buffer size...." << std::endl;

                MPI_Status status;

                int send_tag = i_color;
                int receive_tag = i_color;

                MPI_Sendrecv(send_buffer, send_size, MPI_INT, destination, send_tag,
                             recv_buffers[i_color], recv_size, MPI_INT, destination, receive_tag,
                             MPI_COMM_WORLD, &status);

                delete [] send_buffer;
            }

        // Write in nodes
        for (unsigned int i_color = 0; i_color < neighbours_indices.size(); i_color++)
            if ((destination = neighbours_indices[i_color]) >= 0)
            {
                Communicator::NodesContainerType& r_ghost_nodes = rModelPart.GetCommunicator().InterfaceMesh(i_color).Nodes();

                int* recv_buffer = recv_buffers[i_color];

                unsigned int i = 0;
                for (ModelPart::NodeIterator i_node = r_ghost_nodes.begin(); i_node != r_ghost_nodes.end(); ++i_node)
                {
                    i_node->Set(BOUNDARY, bool(recv_buffer[i++]) );
                }

                if (i > recv_sizes[i_color])
                    std::cout << rank << " Error in estimating receive buffer size...." << std::endl;

                delete[] recv_buffer;
            }
    }

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    TrilinosVariationalDistanceCalculationProcess& operator=(TrilinosVariationalDistanceCalculationProcess const& rOther);

    /// Copy constructor.
    //VariationalDistanceCalculationProcess(VariationalDistanceCalculationProcess const& rOther);


    ///@}

}; // Class TrilinosVariationalDistanceCalculationProcess


///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // KRATOS_TRILINOS_VARIATIONAL_DISTANCE_CALCULATION_PROCESS_INCLUDED  defined


