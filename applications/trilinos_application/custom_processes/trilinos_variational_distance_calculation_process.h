//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotea
//                   Riccardo Rossi
//                   Ruben Zorrilla
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
#include "custom_strategies/builder_and_solvers/trilinos_block_builder_and_solver.h"
#include "processes/variational_distance_calculation_process.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"

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
template< unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver >
class TrilinosVariationalDistanceCalculationProcess
    : public VariationalDistanceCalculationProcess<TDim,TSparseSpace,TDenseSpace,TLinearSolver>
{
public:

    ///@name Type Definitions
    ///@{

    typedef VariationalDistanceCalculationProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
    typedef typename TLinearSolver::Pointer LinearSolverPointerType;
    typedef typename BaseType::SchemeType::Pointer SchemePointerType;
    typedef typename BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer BuilderSolverPointerType;

    ///@}
    ///@name Pointer Definitions
    ///@{

    /// Pointer definition of TrilinosVariationalDistanceCalculationProcess
    KRATOS_CLASS_POINTER_DEFINITION(TrilinosVariationalDistanceCalculationProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    TrilinosVariationalDistanceCalculationProcess(
        Epetra_MpiComm &rComm,
        ModelPart &base_model_part,
        LinearSolverPointerType plinear_solver,
        unsigned int max_iterations = 10) 
        : VariationalDistanceCalculationProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>(base_model_part, max_iterations),
        mrComm(rComm)
    {

        KRATOS_TRY

        // Check that there is at least one element and node in the model
        int NNode = base_model_part.Nodes().size();
        int NElem = base_model_part.Elements().size();

        if (NNode > 0){
            VariableUtils().CheckVariableExists<Variable<double > >(DISTANCE, base_model_part.Nodes());
            VariableUtils().CheckVariableExists<Variable<double > >(FLAG_VARIABLE, base_model_part.Nodes());
        }

        if(NElem > 0){
            if(TDim == 2){
                KRATOS_ERROR_IF(base_model_part.ElementsBegin()->GetGeometry().GetGeometryFamily() != GeometryData::Kratos_Triangle) << 
                    "In 2D the element type is expected to be a triangle" << std::endl;
            } else if(TDim == 3) {
                KRATOS_ERROR_IF(base_model_part.ElementsBegin()->GetGeometry().GetGeometryFamily() != GeometryData::Kratos_Tetrahedra) <<
                    "In 3D the element type is expected to be a tetrahedra" << std::endl;
            }
        }

        base_model_part.GetCommunicator().SumAll(NNode);
        base_model_part.GetCommunicator().SumAll(NElem);

        KRATOS_ERROR_IF(NNode == 0) << "The model has no nodes" << std::endl;
        KRATOS_ERROR_IF(NElem == 0) << "The model has no elements" << std::endl;

        // Generate an auxilary model part and populate it by elements of type DistanceCalculationElementSimplex
        this->ReGenerateDistanceModelPart(base_model_part);

        // Generate a linear strategy
        // Scheme
        SchemePointerType pscheme = Kratos::make_shared<ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace, TDenseSpace> >();

        // Builder and Solver
        const int RowSizeGuess = (TDim == 2 ? 15 : 40);
        BuilderSolverPointerType pBuilderSolver = Kratos::make_shared<TrilinosBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver > >(
            mrComm, 
            RowSizeGuess, 
            plinear_solver);

        // Solution strategy
        const bool CalculateReactions = false;
        const bool ReformDofAtEachIteration = false;
        const bool CalculateNormDxFlag = false;

        ModelPart& r_distance_model_part = base_model_part.GetOwnerModel().GetModelPart("DistancePart");
        (this->mp_solving_strategy) = Kratos::make_unique<ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver> >(
            r_distance_model_part, 
            pscheme, 
            plinear_solver, 
            pBuilderSolver, 
            CalculateReactions, 
            ReformDofAtEachIteration, 
            CalculateNormDxFlag);

        //TODO: check flag DO_EXPENSIVE_CHECKS
        (this->mp_solving_strategy)->Check();

        KRATOS_CATCH("")
    }

    /// Destructor.
    virtual ~TrilinosVariationalDistanceCalculationProcess() {};

    ///@}
    ///@name Operators
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

    void ReGenerateDistanceModelPart(ModelPart& base_model_part) override
    {
        KRATOS_TRY

        // Generate distance model part
        if(!base_model_part.GetOwnerModel().HasModelPart("DistancePart"))
            base_model_part.GetOwnerModel().CreateModelPart("DistancePart",2);

        ModelPart& r_distance_model_part = base_model_part.GetOwnerModel().GetModelPart("DistancePart");
        
        r_distance_model_part.Nodes().clear();
        r_distance_model_part.Conditions().clear();
        r_distance_model_part.Elements().clear();

        r_distance_model_part.SetProcessInfo(base_model_part.pGetProcessInfo());
        r_distance_model_part.SetBufferSize(base_model_part.GetBufferSize());
        r_distance_model_part.SetProperties(base_model_part.pProperties());
        r_distance_model_part.Tables() = base_model_part.Tables();

        // Assigning the nodes to the new model part
        r_distance_model_part.Nodes() = base_model_part.Nodes();

        // Ensure that the nodes have distance as a DOF
        VariableUtils().AddDof<Variable<double> >(DISTANCE, base_model_part);

        // For MPI: copy communication data
        Communicator& rRefComm = base_model_part.GetCommunicator();
        Communicator::Pointer pNewComm = rRefComm.Create();

        pNewComm->SetNumberOfColors(rRefComm.GetNumberOfColors());
        pNewComm->NeighbourIndices() = rRefComm.NeighbourIndices();
        pNewComm->LocalMesh().SetNodes(rRefComm.LocalMesh().pNodes());
        pNewComm->InterfaceMesh().SetNodes(rRefComm.InterfaceMesh().pNodes());
        pNewComm->GhostMesh().SetNodes(rRefComm.GhostMesh().pNodes());
        for (unsigned int i = 0; i < rRefComm.GetNumberOfColors(); ++i){
            pNewComm->pInterfaceMesh(i)->SetNodes(rRefComm.pInterfaceMesh(i)->pNodes());
            pNewComm->pLocalMesh(i)->SetNodes(rRefComm.pLocalMesh(i)->pNodes());
            pNewComm->pGhostMesh(i)->SetNodes(rRefComm.pGhostMesh(i)->pNodes());
        }

        r_distance_model_part.SetCommunicator(pNewComm);

        // Generating the elements
        r_distance_model_part.Elements().reserve(base_model_part.Elements().size());
        for (auto it_elem = base_model_part.ElementsBegin(); it_elem != base_model_part.ElementsEnd(); ++it_elem){
            Element::Pointer p_element = Kratos::make_shared<DistanceCalculationElementSimplex<TDim> >(
                it_elem->Id(),
                it_elem->pGetGeometry(),
                it_elem->pGetProperties());

            // Assign EXACTLY THE SAME GEOMETRY, so that memory is saved!!
            p_element->pGetGeometry() = it_elem->pGetGeometry();

            r_distance_model_part.Elements().push_back(p_element);
            pNewComm->LocalMesh().Elements().push_back(p_element);
        }

        // Using the conditions to mark the boundary with the flag boundary
        // Note that we DO NOT add the conditions to the model part
        VariableUtils().SetFlag<ModelPart::NodesContainerType>(BOUNDARY, false, r_distance_model_part.Nodes());
        // Note that above we have assigned the same geometry. Thus the flag is 
        // set in the distance model part despite we are iterating the base one
        for (auto it_cond = base_model_part.ConditionsBegin(); it_cond != base_model_part.ConditionsEnd(); ++it_cond){
            auto &r_geom = it_cond->GetGeometry();
            for (unsigned int i = 0; i < r_geom.size(); ++i){
                r_geom[i].Set(BOUNDARY, true);
            }
        }

        // Communicate BOUNDARY status to all partitions
        this->CommunicateBoundaryFlagToOwner(r_distance_model_part);
        this->CommunicateBoundaryFlagFromOwner(r_distance_model_part);

        this->mdistance_part_is_initialized = true;

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

                if ((int)i > recv_sizes[i_color])
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

                if ((int)i > recv_sizes[i_color])
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
    TrilinosVariationalDistanceCalculationProcess(TrilinosVariationalDistanceCalculationProcess const &rOther);

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


