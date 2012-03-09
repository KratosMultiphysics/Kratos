/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

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

//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2008-04-30 08:02:32 $
//   Revision:            $Revision: 1.12 $
//
//


#if !defined(KRATOS_MPI_COMMUNICATOR_H_INCLUDED )
#define  KRATOS_MPI_COMMUNICATOR_H_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <sstream>
#include <cstddef>


// External includes 



// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "mpi.h"

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

    /** Detail class definition.
     */
    class MPICommunicator : public Communicator
    {
    public:
        ///@name  Enum's
        ///@{


        ///@}
        ///@name Type Definitions
        ///@{

        /// Pointer definition of MPICommunicator
        KRATOS_CLASS_POINTER_DEFINITION(MPICommunicator);

        typedef Communicator BaseType;

        typedef BaseType::IndexType IndexType;

        typedef BaseType::SizeType SizeType;

        typedef BaseType::NodeType NodeType;

        typedef BaseType::PropertiesType PropertiesType;

        typedef BaseType::ElementType ElementType;

        typedef BaseType::ConditionType ConditionType;

        typedef BaseType::NeighbourIndicesContainerType NeighbourIndicesContainerType;

        typedef BaseType::MeshType MeshType;

        typedef BaseType::MeshesContainerType MeshesContainerType;

        /// Nodes container. Which is a vector set of nodes with their Id's as key.
        typedef MeshType::NodesContainerType NodesContainerType;

        /** Iterator over the nodes. This iterator is an indirect
            iterator over Node::Pointer which turn back a reference to
            node by * operator and not a pointer for more convenient
            usage. */
        typedef MeshType::NodeIterator NodeIterator;

        /** Const iterator over the nodes. This iterator is an indirect
            iterator over Node::Pointer which turn back a reference to
            node by * operator and not a pointer for more convenient
            usage. */
        typedef MeshType::NodeConstantIterator NodeConstantIterator;

        /** Iterator over the properties. This iterator is an indirect
            iterator over Properties::Pointer which turn back a reference to
            properties by * operator and not a pointer for more convenient
            usage. */

        /// Properties container. Which is a vector set of Properties with their Id's as key.
        typedef MeshType::PropertiesContainerType PropertiesContainerType;

        /** Iterator over the Properties. This iterator is an indirect
            iterator over Node::Pointer which turn back a reference to
            node by * operator and not a pointer for more convenient
            usage. */
        typedef MeshType::PropertiesIterator PropertiesIterator;

        /** Const iterator over the Properties. This iterator is an indirect
            iterator over Properties::Pointer which turn back a reference to
            Properties by * operator and not a pointer for more convenient
            usage. */
        typedef MeshType::PropertiesConstantIterator PropertiesConstantIterator;

        /** Iterator over the properties. This iterator is an indirect
            iterator over Properties::Pointer which turn back a reference to
            properties by * operator and not a pointer for more convenient
            usage. */

        /// Element container. A vector set of Elements with their Id's as key.
        typedef MeshType::ElementsContainerType ElementsContainerType;

        /** Iterator over the Elements. This iterator is an indirect
            iterator over Elements::Pointer which turn back a reference to
            Element by * operator and not a pointer for more convenient
            usage. */
        typedef MeshType::ElementIterator ElementIterator;

        /** Const iterator over the Elements. This iterator is an indirect
            iterator over Elements::Pointer which turn back a reference to
            Element by * operator and not a pointer for more convenient
            usage. */
        typedef MeshType::ElementConstantIterator ElementConstantIterator;

        /// Condintions container. A vector set of Conditions with their Id's as key.
        typedef MeshType::ConditionsContainerType ConditionsContainerType;

        /** Iterator over the Conditions. This iterator is an indirect
           iterator over Conditions::Pointer which turn back a reference to
           Condition by * operator and not a pointer for more convenient
           usage. */
        typedef MeshType::ConditionIterator ConditionIterator;

        /** Const iterator over the Conditions. This iterator is an indirect
            iterator over Conditions::Pointer which turn back a reference to
            Condition by * operator and not a pointer for more convenient
            usage. */
        typedef MeshType::ConditionConstantIterator ConditionConstantIterator;

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.

        MPICommunicator() : BaseType()
        {
        }

        /// Copy constructor.

        MPICommunicator(MPICommunicator const& rOther) : BaseType(rOther)
        {
        }


        /// Destructor.

        virtual ~MPICommunicator()
        {
        }

        virtual Communicator::Pointer Create()
        {
            KRATOS_TRY

            return Communicator::Pointer(new MPICommunicator);

            KRATOS_CATCH("");
        }


        ///@}
        ///@name Operators
        ///@{

        /// Assignment operator.

        MPICommunicator & operator=(MPICommunicator const& rOther)
        {
            BaseType::operator=(rOther);
            return *this;
        }

        int MyPID()
        {
            int rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            return rank;
        }

        int TotalProcesses()
        {
            int nproc;
            MPI_Comm_size(MPI_COMM_WORLD, &nproc);
            return nproc;

        }

        ///@}
        ///@name Access
        ///@{


        ///@}
        ///@name Operations
        ///@{

        void Barrier()
        {
            MPI_Barrier(MPI_COMM_WORLD);
        }

        virtual bool SumAll(int& rValue)
        {
            int local_value = rValue;
            MPI_Allreduce(&local_value, &rValue, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            return true;
        }

        virtual bool SumAll(double& rValue)
        {
            double local_value = rValue;
            MPI_Allreduce(&local_value, &rValue, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            return true;
        }

        virtual bool MinAll(int& rValue)
        {
            int local_value = rValue;
            MPI_Allreduce(&local_value, &rValue, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
            return true;
        }

        virtual bool MinAll(double& rValue)
        {
            double local_value = rValue;
            MPI_Allreduce(&local_value, &rValue, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
            return true;
        }

        virtual bool MaxAll(int& rValue)
        {
            int local_value = rValue;
            MPI_Allreduce(&local_value, &rValue, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
            return true;
        }

        virtual bool MaxAll(double& rValue)
        {
            double local_value = rValue;
            MPI_Allreduce(&local_value, &rValue, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            return true;
        }

        virtual bool SynchronizeNodalSolutionStepsData()
        {
            int rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);

            int destination = 0;

            NeighbourIndicesContainerType& neighbours_indices = NeighbourIndices();

            for (unsigned int i_color = 0; i_color < neighbours_indices.size(); i_color++)
                if ((destination = neighbours_indices[i_color]) >= 0)
                {
                    NodesContainerType& r_local_nodes = LocalMesh(i_color).Nodes();
                    NodesContainerType& r_ghost_nodes = GhostMesh(i_color).Nodes();

                    // Calculating send and received buffer size
                    // NOTE: This part works ONLY when all nodes have the same variables list size!
                    unsigned int nodal_data_size = 0;
                    unsigned int local_nodes_size = r_local_nodes.size();
                    unsigned int ghost_nodes_size = r_ghost_nodes.size();
                    unsigned int send_buffer_size = 0;
                    unsigned int receive_buffer_size = 0;

                    if (local_nodes_size == 0)
                    {
                        if (ghost_nodes_size == 0)
                            continue; // nothing to transfer!
                        else
                        {
                            nodal_data_size = r_ghost_nodes.begin()->SolutionStepData().TotalSize();
                            receive_buffer_size = ghost_nodes_size * nodal_data_size;
                        }
                    } else
                    {
                        nodal_data_size = r_local_nodes.begin()->SolutionStepData().TotalSize();
                        send_buffer_size = local_nodes_size * nodal_data_size;
                        if (ghost_nodes_size != 0)
                            receive_buffer_size = ghost_nodes_size * nodal_data_size;
                    }

                    unsigned int position = 0;
                    double* send_buffer = new double[send_buffer_size];
                    double* receive_buffer = new double[receive_buffer_size];


                    // Filling the buffer
                    for (ModelPart::NodeIterator i_node = r_local_nodes.begin(); i_node != r_local_nodes.end(); ++i_node)
                    {
                        std::memcpy(send_buffer + position, i_node->SolutionStepData().Data(), nodal_data_size * sizeof (double));
                        position += nodal_data_size;
                    }

                    MPI_Status status;

                    if (position > send_buffer_size)
                        std::cout << rank << " Error in estimating send buffer size...." << std::endl;


                    int send_tag = i_color;
                    int receive_tag = i_color;


                    MPI_Sendrecv(send_buffer, send_buffer_size, MPI_DOUBLE, destination, send_tag, receive_buffer, receive_buffer_size, MPI_DOUBLE, destination, receive_tag,
                            MPI_COMM_WORLD, &status);

                    // Updating nodes
                    position = 0;
                    for (ModelPart::NodeIterator i_node = GhostMesh(i_color).NodesBegin();
                            i_node != GhostMesh(i_color).NodesEnd(); i_node++)
                    {
                        std::memcpy(i_node->SolutionStepData().Data(), receive_buffer + position, nodal_data_size * sizeof (double));
                        position += nodal_data_size;
                    }

                    if (position > receive_buffer_size)
                        std::cout << rank << " Error in estimating receive buffer size...." << std::endl;

                    delete [] send_buffer;
                    delete [] receive_buffer;
                }

            return true;
        }

        virtual bool SynchronizeDofs()
        {
            int rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);

            int destination = 0;

            NeighbourIndicesContainerType& neighbours_indices = NeighbourIndices();

            for (unsigned int i_color = 0; i_color < neighbours_indices.size(); i_color++)
                if ((destination = neighbours_indices[i_color]) >= 0)
                {
                    NodesContainerType& r_local_nodes = LocalMesh(i_color).Nodes();
                    NodesContainerType& r_ghost_nodes = GhostMesh(i_color).Nodes();

                    // Calculating send and received buffer size
                    unsigned int send_buffer_size = 0;
                    unsigned int receive_buffer_size = 0;

                    for (NodesContainerType::iterator i_node = r_local_nodes.begin(); i_node != r_local_nodes.end(); ++i_node)
                        send_buffer_size += i_node->GetDofs().size();

                    for (NodesContainerType::iterator i_node = r_ghost_nodes.begin(); i_node != r_ghost_nodes.end(); ++i_node)
                        receive_buffer_size += i_node->GetDofs().size();

                    unsigned int position = 0;
                    int* send_buffer = new int[send_buffer_size];
                    int* receive_buffer = new int[receive_buffer_size];


                    // Filling the buffer
                    for (ModelPart::NodeIterator i_node = r_local_nodes.begin(); i_node != r_local_nodes.end(); ++i_node)
                        for (ModelPart::NodeType::DofsContainerType::iterator i_dof = i_node->GetDofs().begin(); i_dof != i_node->GetDofs().end(); i_dof++)
                        {
                            send_buffer[position++] = i_dof->EquationId();
                        }


                    MPI_Status status;

                    if (position > send_buffer_size)
                        std::cout << rank << " Error in estimating send buffer size...." << std::endl;


                    int send_tag = i_color;
                    int receive_tag = i_color;



                    MPI_Sendrecv(send_buffer, send_buffer_size, MPI_INT, destination, send_tag, receive_buffer, receive_buffer_size, MPI_INT, destination, receive_tag,
                            MPI_COMM_WORLD, &status);

                    // Updating nodes
                    position = 0;
                    for (ModelPart::NodeIterator i_node = GhostMesh(i_color).NodesBegin();
                            i_node != GhostMesh(i_color).NodesEnd(); i_node++)
                        for (ModelPart::NodeType::DofsContainerType::iterator i_dof = i_node->GetDofs().begin(); i_dof != i_node->GetDofs().end(); i_dof++)
                            i_dof->SetEquationId(receive_buffer[position++]);

                    if (position > receive_buffer_size)
                        std::cout << rank << " Error in estimating receive buffer size...." << std::endl;

                    delete [] send_buffer;
                    delete [] receive_buffer;
                }

            return true;
        }

        virtual bool AssembleCurrentData(Variable<double> const& ThisVariable)
        {
            AssembleThisVariable(ThisVariable);
            return true;
        }

        virtual bool AssembleCurrentData(Variable<array_1d<double, 3 > > const& ThisVariable)
        {
            AssembleThisVariable(ThisVariable);
            return true;
        }

        virtual bool AssembleCurrentData(Variable<Vector> const& ThisVariable)
        {
            AssembleThisVariable(ThisVariable);
            return true;
        }

        virtual bool AssembleCurrentData(Variable<Matrix> const& ThisVariable)
        {
            AssembleThisVariable(ThisVariable);
            return true;
        }

        ///@}
        ///@name Access
        ///@{


        ///@}
        ///@name Inquiry
        ///@{


        ///@}
        ///@name Input and output
        ///@{

        /// Turn back information as a string.

        virtual std::string Info() const
        {
            return "MPICommunicator";
        }

        /// Print information about this object.

        virtual void PrintInfo(std::ostream& rOStream) const
        {
            rOStream << Info();
        }

        /// Print object's data.

        virtual void PrintData(std::ostream& rOStream) const
        {
            for (IndexType i = 0; i < mLocalMeshes.size(); i++)
            {
                rOStream << "    Local Mesh " << i << " : " << std::endl;
                LocalMesh(i).PrintData(rOStream);
                rOStream << "    Ghost Mesh " << i << " : " << std::endl;
                GhostMesh(i).PrintData(rOStream);
                rOStream << "    Interface Mesh " << i << " : " << std::endl;
                InterfaceMesh(i).PrintData(rOStream);
            }
        }


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
        ///@name Protected Operators
        ///@{


        ///@}
        ///@name Protected Operations
        ///@{


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


        //      SizeType mNumberOfColors;

        //      NeighbourIndicesContainerType mNeighbourIndices;
        //
        //      // To store all local entities
        //      MeshType::Pointer mpLocalMesh;
        //
        //      // To store all ghost entities
        //      MeshType::Pointer mpGhostMesh;
        //
        //      // To store all interface entities
        //      MeshType::Pointer mpInterfaceMesh;
        //
        //      // To store interfaces local entities
        //      MeshesContainerType mLocalMeshes;

        //      // To store interfaces ghost entities
        //      MeshesContainerType mGhostMeshes;
        //
        //      // To store interfaces ghost+local entities
        //      MeshesContainerType mInterfaceMeshes;

        ///@}
        ///@name Private Operators
        ///@{


        ///@}
        ///@name Private Operations
        ///@{

        void PrintNodesId()
        {
            NeighbourIndicesContainerType& neighbours_indices = NeighbourIndices();

            int nproc;
            MPI_Comm_size(MPI_COMM_WORLD, &nproc);
            int rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);

            MPI_Barrier(MPI_COMM_WORLD);
            for (int proc_id = 0; proc_id < nproc; proc_id++)
            {
                if (proc_id == rank)
                {

                    for (int i_color = 0; i_color < static_cast<int>(neighbours_indices.size()); i_color++)
                    {
                        if ((neighbours_indices[i_color]) >= 0)
                        {
                            NodesContainerType& r_local_nodes = LocalMesh(i_color).Nodes();
                            NodesContainerType& r_ghost_nodes = GhostMesh(i_color).Nodes();
                            std::string tag = "Local nodes in rank ";
                            PrintNodesId(r_local_nodes, tag, i_color);
                            tag = "Ghost nodes in rank ";
                            PrintNodesId(r_ghost_nodes, tag, i_color);
                            tag = "Interface nodes in rank ";
                            PrintNodesId(InterfaceMesh(i_color).Nodes(), tag, i_color);

                        }
                    }
                }
                MPI_Barrier(MPI_COMM_WORLD);
            }
        }

        template<class TNodesArrayType>
        void PrintNodesId(TNodesArrayType& rNodes, std::string Tag, int color)
        {
            int rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            std::cout << Tag << rank << " with color " << color << ":";
            for (typename TNodesArrayType::iterator i_node = rNodes.begin(); i_node != rNodes.end(); i_node++)
                std::cout << i_node->Id() << ", ";

            std::cout << std::endl;
        }

        template<class TDataType>
        bool AssembleThisVariable(Variable<TDataType> const& ThisVariable)
        {
            // PrintNodesId();
            /*	KRATOS_WATCH("AssembleThisVariable")
                    KRATOS_WATCH(ThisVariable)*/
            int rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);

            int destination = 0;

            NeighbourIndicesContainerType& neighbours_indices = NeighbourIndices();
            std::vector<double*> receive_buffer(neighbours_indices.size());
            std::vector<int> receive_buffer_size(neighbours_indices.size());

            //first of all gather everything to the owner node
            for (unsigned int i_color = 0; i_color < neighbours_indices.size(); i_color++)
                if ((destination = neighbours_indices[i_color]) >= 0)
                {
                    NodesContainerType& r_local_nodes = InterfaceMesh(i_color).Nodes();
                    NodesContainerType& r_ghost_nodes = InterfaceMesh(i_color).Nodes();

                    // Calculating send and received buffer size
                    // NOTE: This part can be optimized getting the offset from variables list and using pointers.
                    unsigned int nodal_data_size = sizeof (TDataType) / sizeof (double);
                    unsigned int local_nodes_size = r_local_nodes.size();
                    unsigned int ghost_nodes_size = r_ghost_nodes.size();
                    unsigned int send_buffer_size = local_nodes_size * nodal_data_size;
                    receive_buffer_size[i_color] = ghost_nodes_size * nodal_data_size;

                    if ((local_nodes_size == 0) && (ghost_nodes_size == 0))
                        continue; // nothing to transfer!

                    unsigned int position = 0;
                    double* send_buffer = new double[send_buffer_size];
                    receive_buffer[i_color] = new double[receive_buffer_size[i_color]];

                    // Filling the buffer
                    for (ModelPart::NodeIterator i_node = r_local_nodes.begin(); i_node != r_local_nodes.end(); ++i_node)
                    {
                        *(TDataType*) (send_buffer + position) = i_node->FastGetSolutionStepValue(ThisVariable);
                        position += nodal_data_size;
                    }

                    MPI_Status status;

                    if (position > send_buffer_size)
                        std::cout << rank << " Error in estimating send buffer size...." << std::endl;

                    int send_tag = i_color;
                    int receive_tag = i_color;

                    MPI_Sendrecv(send_buffer, send_buffer_size, MPI_DOUBLE, destination, send_tag,
                            receive_buffer[i_color], receive_buffer_size[i_color], MPI_DOUBLE, destination, receive_tag,
                            MPI_COMM_WORLD, &status);

                    delete [] send_buffer;
                }

            for (unsigned int i_color = 0; i_color < neighbours_indices.size(); i_color++)
                if ((destination = neighbours_indices[i_color]) >= 0)
                {
                    // Updating nodes
                    int position = 0;
                    unsigned int nodal_data_size = sizeof (TDataType) / sizeof (double);
                    NodesContainerType& r_ghost_nodes = InterfaceMesh(i_color).Nodes();

                    for (ModelPart::NodeIterator i_node = r_ghost_nodes.begin(); i_node != r_ghost_nodes.end(); ++i_node)
                    {
                        i_node->FastGetSolutionStepValue(ThisVariable) += *reinterpret_cast<TDataType*> (receive_buffer[i_color] + position);
                        position += nodal_data_size;
                    }

                    if (position > receive_buffer_size[i_color])
                        std::cout << rank << " Error in estimating receive buffer size...." << std::endl;

                    delete [] receive_buffer[i_color];
                }

            //MPI_Barrier(MPI_COMM_WORLD);


            //SynchronizeNodalSolutionStepsData();
            SynchronizeVariable(ThisVariable);




            return true;
        }

        template<class TDataType>
        bool SynchronizeVariable(Variable<TDataType> const& ThisVariable)
        {
            int rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);

            int destination = 0;

            NeighbourIndicesContainerType& neighbours_indices = NeighbourIndices();

            for (unsigned int i_color = 0; i_color < neighbours_indices.size(); i_color++)
                if ((destination = neighbours_indices[i_color]) >= 0)
                {
                    NodesContainerType& r_local_nodes = LocalMesh(i_color).Nodes();
                    NodesContainerType& r_ghost_nodes = GhostMesh(i_color).Nodes();

                    unsigned int nodal_data_size = sizeof (TDataType) / sizeof (double);
                    unsigned int local_nodes_size = r_local_nodes.size();
                    unsigned int ghost_nodes_size = r_ghost_nodes.size();
                    unsigned int send_buffer_size = local_nodes_size * nodal_data_size;
                    unsigned int receive_buffer_size = ghost_nodes_size * nodal_data_size;

                    if ((local_nodes_size == 0) && (ghost_nodes_size == 0))
                        continue; // nothing to transfer!

                    unsigned int position = 0;
                    double* send_buffer = new double[send_buffer_size];
                    double* receive_buffer = new double[receive_buffer_size];

                    // Filling the send buffer
                    for (ModelPart::NodeIterator i_node = r_local_nodes.begin(); i_node != r_local_nodes.end(); ++i_node)
                    {
                        *(TDataType*) (send_buffer + position) = i_node->FastGetSolutionStepValue(ThisVariable);
                        position += nodal_data_size;
                    }

                    MPI_Status status;

                    int send_tag = i_color;
                    int receive_tag = i_color;

                    MPI_Sendrecv(send_buffer, send_buffer_size, MPI_DOUBLE, destination, send_tag, receive_buffer, receive_buffer_size, MPI_DOUBLE, destination, receive_tag,
                            MPI_COMM_WORLD, &status);

                    position = 0;
                    for (ModelPart::NodeIterator i_node = r_ghost_nodes.begin(); i_node != r_ghost_nodes.end(); ++i_node)
                    {
                        i_node->FastGetSolutionStepValue(ThisVariable) = *reinterpret_cast<TDataType*> (receive_buffer + position);
                        position += nodal_data_size;
                    }

                    if (position > receive_buffer_size)
                        std::cout << rank << " Error in estimating receive buffer size...." << std::endl;

                    delete [] send_buffer;
                    delete [] receive_buffer;
                }

            return true;
        }



        //       friend class boost::serialization::access;

        //       template<class TArchive>
        // 	  void serialize(TArchive & ThisArchive, const unsigned int ThisVersion)
        // 	  {
        // /* 	      ThisArchive & mName & mBufferSize & mCurrentIndex; */
        // 	  }

        //       void RemoveSolutionStepData(IndexType SolutionStepIndex, MeshType& ThisMesh)
        // 	{
        // 	  for(NodeIterator i_node = ThisMesh.NodesBegin() ; i_node != ThisMesh.NodesEnd() ; ++i_node)
        // 	    i_node->RemoveSolutionStepNodalData(SolutionStepIndex);
        // 	}

        ///@}
        ///@name Private  Access
        ///@{


        ///@}
        ///@name Private Inquiry
        ///@{


        ///@}
        ///@name Un accessible methods
        ///@{


        ///@}

    }; // Class MPICommunicator 

    ///@}

    ///@name Type Definitions
    ///@{


    ///@}
    ///@name Input and output
    ///@{


    /// input stream function
    inline std::istream & operator >>(std::istream& rIStream,
            MPICommunicator& rThis);

    /// output stream function

    inline std::ostream & operator <<(std::ostream& rOStream,
            const MPICommunicator& rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }
    ///@}


} // namespace Kratos.

#endif // KRATOS_MPI_COMMUNICATOR_H_INCLUDED  defined 


