//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Jordi Cotela
//                   Carlos Roig
//


#if !defined(KRATOS_METIS_CONTACT_PARTITIONING_PROCESS_INCLUDED )
#define  KRATOS_METIS_CONTACT_PARTITIONING_PROCESS_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <algorithm>
#include <fstream>

// External includes
#include <parmetis.h>


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "processes/graph_coloring_process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/mpi_communicator.h"

extern "C" {
    //extern void METIS_PartMeshDual(int*, int*, idxtype*, int*, int*, int*, int*, idxtype*, idxtype*);
    extern int METIS_PartMeshDual(int*, int*, int*, int*, int*, int*, int*, int*, int*);
}



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

class MetisContactPartitioningProcess
    : public MetisPartitioningProcess
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MetisContactPartitioningProcess
    KRATOS_CLASS_POINTER_DEFINITION(MetisContactPartitioningProcess);

    typedef std::size_t SizeType;
    typedef std::size_t IndexType;
    typedef matrix<int> GraphType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MetisContactPartitioningProcess(ModelPart& rModelPart, IO& rIO, SizeType NumberOfPartitions, std::vector<int> const& ContactNodes, int Dimension = 3)
        : MetisPartitioningProcess(rModelPart, rIO, NumberOfPartitions, Dimension) , mContactNodes(ContactNodes)
    {
        KRATOS_TRY
        int rank = GetRank();

        std::stringstream log_filename;
        log_filename << "kratos_metis_" << rank << ".log";
//		    mLogFile.open(log_filename.str().c_str());
        KRATOS_CATCH("")
    }

    /// Copy constructor.
    MetisContactPartitioningProcess(MetisContactPartitioningProcess const& rOther)
        : MetisPartitioningProcess(rOther) , mContactNodes(rOther.mContactNodes)
    {
        KRATOS_TRY
        int rank = GetRank();

        std::stringstream log_filename;
        log_filename << "kratos_metis_" << rank << ".log";
//		    mLogFile.open(log_filename.str().c_str());
        KRATOS_CATCH("")
    }

    /// Destructor.
    virtual ~MetisContactPartitioningProcess()
    {
    }


    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
        Execute();
    }


    ///@}
    ///@name Operations
    ///@{

    void Execute() override
    {
        KRATOS_TRY;

        int number_of_processes;
        MPI_Comm_size (MPI_COMM_WORLD,&number_of_processes);

        if(number_of_processes < 2) // There is no need to partition it and just reading the input
        {
            mrIO.ReadModelPart(mrModelPart);
            return;
        }

        // Set MPICommunicator as modelpart's communicator
        VariablesList * mVariables_List = &mrModelPart.GetNodalSolutionStepVariablesList();
        mrModelPart.SetCommunicator(Communicator::Pointer(new MPICommunicator(mVariables_List)));

        // if mNumberOfPartitions is not defined we set it to the number_of_processes
        if (mNumberOfPartitions == 0)
            mNumberOfPartitions = static_cast<SizeType>(number_of_processes);

        KRATOS_WATCH(mNumberOfPartitions);

        int rank = GetRank();


        // Reading connectivities
        IO::ConnectivitiesContainerType elements_connectivities;
        int number_of_elements =  mrIO.ReadElementsConnectivities(elements_connectivities);

        mLogFile << rank << " :Reading nodes" << std::endl;
        ModelPart::NodesContainerType temp_nodes;
        mrIO.ReadNodes(temp_nodes);

        int number_of_nodes = temp_nodes.size(); // considering sequencial numbering!!

        idxtype* epart = new idxtype[number_of_elements];
        idxtype* npart = new idxtype[number_of_nodes];

        GraphType domains_graph = zero_matrix<int>(mNumberOfPartitions, mNumberOfPartitions);
        GraphType domains_colored_graph;
        int* coloring_send_buffer = NULL;

        // Adding interface meshes
        mrModelPart.GetMeshes().push_back(Kratos::make_shared<ModelPart::MeshType>());

        int colors_number;


        if(rank == 0)
        {
            CallingMetis(number_of_nodes, number_of_elements, elements_connectivities, npart, epart);
            for(unsigned int i = 0 ; i < mContactNodes.size() ; i++)
            {
                npart[mContactNodes[i] - 1] = mNumberOfPartitions - 1;
            }
            CalculateDomainsGraph(domains_graph, number_of_elements, elements_connectivities, npart, epart);
            GraphColoringProcess(mNumberOfPartitions, domains_graph,domains_colored_graph, colors_number).Execute();
// 		      colors_number = GraphColoring(domains_graph, domains_colored_graph);
            KRATOS_WATCH(colors_number);
            KRATOS_WATCH(domains_graph);
            KRATOS_WATCH(domains_colored_graph);

            // Filling the sending buffer
            int buffer_index = 0;
            coloring_send_buffer = new int[mNumberOfPartitions*colors_number];
            for(unsigned int i = 0 ; i <  mNumberOfPartitions ; i++)
                for(int j = 0 ; j <  colors_number ; j++)
                    coloring_send_buffer[buffer_index++] = domains_colored_graph(i,j);

            mLogFile << rank << "  : colors_number = " << colors_number << std::endl;

            mLogFile << rank << " : coloring_send_buffer = [";
            for(SizeType j = 0 ; j < mNumberOfPartitions*colors_number ; j++)
                mLogFile << coloring_send_buffer[j] << " ,";
            mLogFile << "]" << std::endl;
        }

        // Broadcasting partioning information
        MPI_Bcast(&colors_number, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(npart, number_of_nodes, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(epart, number_of_elements, MPI_INT, 0, MPI_COMM_WORLD);


        /* 		  vector<int> neighbours_indices(colors_number); */





        Communicator::NeighbourIndicesContainerType& neighbours_indices = mrModelPart.GetCommunicator().NeighbourIndices();
        if(neighbours_indices.size() != static_cast<unsigned int>(colors_number))
            neighbours_indices.resize(colors_number,false);
        for(int i = 0 ; i  < colors_number ; i++)
            neighbours_indices[i] = 0;

        mLogFile << rank << " : neighbours_indices = [";
        for(unsigned int j = 0 ; j < neighbours_indices.size()  ; j++)
            mLogFile << neighbours_indices[j] << " ,";
        mLogFile << "]" << std::endl;

        mLogFile << rank << "  : colors_number = " << colors_number << std::endl;
        mrModelPart.GetCommunicator().SetNumberOfColors(colors_number);


        MPI_Scatter(coloring_send_buffer,colors_number,MPI_INT,&(neighbours_indices[0]),colors_number,MPI_INT,0,MPI_COMM_WORLD);



        mLogFile << rank << " : ";
        for(int j = 0 ; j <  colors_number ; j++)
            mLogFile << mrModelPart.GetCommunicator().NeighbourIndices()[j] << " ,";
        mLogFile << "]" << std::endl;

        // Adding local, ghost and interface meshes to ModelPart if is necessary
        int number_of_meshes =  ModelPart::Kratos_Ownership_Size + colors_number; // (all + local + ghost) + (colors_number for interfaces)
        if(mrModelPart.GetMeshes().size() < static_cast<unsigned int>(number_of_meshes))
            for(int i = mrModelPart.GetMeshes().size() ; i < number_of_meshes ; i++)
	      mrModelPart.GetMeshes().push_back(Kratos::make_shared<ModelPart::MeshType>());

        for(ModelPart::NodeIterator i_node = temp_nodes.begin() ;
                i_node != temp_nodes.end() ; i_node++)
            i_node->SetSolutionStepVariablesList(&(mrModelPart.GetNodalSolutionStepVariablesList()));

        // Adding nodes to modelpart
        AddingNodes(temp_nodes, number_of_elements, elements_connectivities, npart, epart);

        mLogFile << rank << " : Start reading Properties " << std::endl;

        // Adding properties to modelpart
        mrIO.ReadProperties(mrModelPart.rProperties());

        mLogFile << rank << " : End adding Properties " << std::endl;

        // Adding elements to each partition mesh
        AddingElements(temp_nodes, npart, epart);


        // Adding conditions to each partition mesh
        ModelPart::ConditionsContainerType temp_conditions;
        AddingConditions(temp_nodes, npart, epart, temp_conditions);



        mrIO.ReadInitialValues(temp_nodes, mrModelPart.Elements(), temp_conditions);

        mLogFile << rank << " : start cleaning memory " << std::endl;


        delete[] epart;
        delete[] npart;
        if(rank == 0)
        {
            mLogFile << rank << " : deleting coloring_send_buffer " << std::endl;
            delete[] coloring_send_buffer;
        }

        mLogFile << rank << " : cleaning memory Finished" << std::endl;

        KRATOS_CATCH("")
    }

    /*		void CalculateDomainsGraph(GraphType& rDomainsGraph, SizeType NumberOfElements, IO::ConnectivitiesContainerType& ElementsConnectivities, idxtype* NPart, idxtype* EPart )
    		  {
    		    for(SizeType i_element = 0 ; i_element < NumberOfElements ; i_element++)
    		      for(std::vector<std::size_t>::iterator i_node = ElementsConnectivities[i_element].begin() ;
    			  i_node != ElementsConnectivities[i_element].end() ; i_node++)
    			{
    			  SizeType node_rank = NPart[*i_node-1];
    			  SizeType element_rank = EPart[i_element];
    			  if(node_rank != element_rank)
    			  {
    			    rDomainsGraph(node_rank, element_rank) = 1;
    			    rDomainsGraph(element_rank, node_rank) = 1;
    			  }
    			}

    KRATOS_WATCH(rDomainsGraph)

    		  }
    */
// 		int GraphColoring(GraphType& rDomainsGraph, GraphType& rDomainsColoredGraph)
// 		  {
// 		    int max_color = 0;
// 		    // Initializing the coloered graph. -1 means no connection
// 		    rDomainsColoredGraph = ScalarMatrix(mNumberOfPartitions, mNumberOfPartitions, -1.00);

// 		    // Start coloring...
// 		    for(SizeType i = 0 ; i < rDomainsGraph.size1() ; i++) // for each domain
// 		      for(SizeType j = i + 1 ; j < rDomainsGraph.size2() ; j++) // finding neighbor domains
// 			if(rDomainsGraph(i,j) != 0.00) // domain i has interface with domain j
// 			  for(SizeType color = 0 ; color < rDomainsColoredGraph.size2() ; color++) // finding color
// 			    if((rDomainsColoredGraph(i,color) == -1.00) && (rDomainsColoredGraph(j,color) == -1.00)) // the first unused color
// 			    {
// 			      rDomainsColoredGraph(i,color) = j;
// 			      rDomainsColoredGraph(j,color) = i;
// 			      if(max_color < static_cast<int>(color + 1))
// 				max_color = color + 1;
// 			      break;
// 			    }
// 		    return max_color;
// 		  }

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
    std::string Info() const override
    {
        return "MetisContactPartitioningProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MetisContactPartitioningProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
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

    void CallingMetis(SizeType NumberOfNodes, SizeType NumberOfElements, IO::ConnectivitiesContainerType& ElementsConnectivities, idxtype* NPart, idxtype* EPart)
    {
        int rank = GetRank();

        // calculating total size of connectivity vector
        int connectivity_size = 0;
        for( IO::ConnectivitiesContainerType::iterator i_connectivities = ElementsConnectivities.begin(); i_connectivities != ElementsConnectivities.end(); i_connectivities++ )
            connectivity_size += i_connectivities->size();

        int number_of_element_nodes = ElementsConnectivities.begin()->size(); // here assuming that all elements are the same!!

        int ne = NumberOfElements;
        int nn = NumberOfNodes;

        int quadratic_type = 0;
        int etype;
        if(number_of_element_nodes == 3) // triangles
            etype = 1;
        else if(number_of_element_nodes == 4) // tetrahedra or quadilateral
        {
            if(mDimension == 2) // quadilateral
                etype = 4;
            else  // tetrahedra
                etype = 2;
        }
        else if(number_of_element_nodes == 10) //quadratic tetrahedra
        {
            etype = 2;
            quadratic_type = 4;
            connectivity_size = ne*quadratic_type;
        }
        else if(number_of_element_nodes == 8) // hexahedra
            etype = 3;
        else if(number_of_element_nodes == 20 || number_of_element_nodes == 27) //quadratic hexahedra
        {
            etype = 3;
            quadratic_type = 8;
            connectivity_size = ne*quadratic_type;
        }
        else
            KRATOS_THROW_ERROR(std::invalid_argument, "invalid element type with number of nodes : ", number_of_element_nodes);

        int numflag = 0;
        int number_of_partitions = static_cast<int>(mNumberOfPartitions) - 1; // we reserve one partition for contact
        int edgecut;

        idxtype* elmnts = new idxtype[connectivity_size];

        mLogFile << rank << " : Preparing Data for metis..." << std::endl;
        int i = 0;

        //handle quadratic elements
        if( quadratic_type != 0 )
        {
            idxtype index = 0;
            idxtype* new_node_index = new idxtype[NumberOfNodes];
            for( unsigned int k=0; k<NumberOfNodes; k++ )
                new_node_index[k] = NumberOfNodes+1;
            for(IO::ConnectivitiesContainerType::iterator i_connectivities = ElementsConnectivities.begin() ; i_connectivities != ElementsConnectivities.end() ; i_connectivities++)
            {
                for(int j = 0 ; j < quadratic_type ; j++)
                {
                    if( static_cast<SizeType>(new_node_index[(*i_connectivities)[j]-1]) == NumberOfNodes+1 )
                        new_node_index[(*i_connectivities)[j]-1] = index++;
                    elmnts[i++] = new_node_index[(*i_connectivities)[j]-1]; // transforming to zero base indexing
                }
            }
            nn = index;

            mLogFile << rank << " : Calling metis..." << std::endl;
            // Calling Metis to partition
            METIS_PartMeshDual(&ne, &nn, elmnts, &etype, &numflag, &number_of_partitions, &edgecut, EPart, NPart);
            mLogFile << rank << " : Metis Finished!!!" << std::endl;
            mLogFile << rank << " :     edgecut = " << edgecut << std::endl;

            //distribution of nodeal partition indices by elemental partition indices
            for(unsigned int i=0; i<NumberOfElements; i++)
            {
                for(unsigned int j = 0 ; j < ElementsConnectivities[i].size() ; j++)
                    NPart[(ElementsConnectivities[i][j]-1)] = EPart[i]; // transforming to zero base indexing
            }
        }//end quadratic branch
        else //linear branch
        {
            // Creating the elmnts array for Metis
            for(IO::ConnectivitiesContainerType::iterator i_connectivities = ElementsConnectivities.begin(); i_connectivities != ElementsConnectivities.end(); i_connectivities++)
                for(unsigned int j = 0 ; j < i_connectivities->size() ; j++)
                    elmnts[i++] = (*i_connectivities)[j] - 1; // transforming to zero base indexing

            mLogFile << rank << " : Calling metis..." << std::endl;
            // Calling Metis to partition
            METIS_PartMeshDual(&ne, &nn, elmnts, &etype, &numflag, &number_of_partitions, &edgecut, EPart, NPart);
            mLogFile << rank << " : Metis Finished!!!" << std::endl;
            mLogFile << rank << " :     edgecut = " << edgecut << std::endl;
        }
        //deallocate memory
        delete[] elmnts;
    }

    /*	        void AddingNodes(ModelPart::NodesContainerType& AllNodes, SizeType NumberOfElements, IO::ConnectivitiesContainerType& ElementsConnectivities, idxtype* NPart, idxtype* EPart)
    	        {
    		  int rank = GetRank();

    		  mLogFile << rank << " : Adding nodes to modelpart" << std::endl;
    		  // first adding the partition's nodes
    		  for(ModelPart::NodeIterator i_node = AllNodes.begin() ;
    		      i_node != AllNodes.end() ; i_node++)
    		    if(NPart[i_node->Id()-1] == rank)
    		      {
    			mrModelPart.AssignNode(*(i_node.base()));
     			mrModelPart.GetCommunicator().LocalMesh(). AddNode(*(i_node.base()));
    //  			mrModelPart.AssignNode(*(i_node.base()), ModelPart::Kratos_Local);
    			i_node->GetSolutionStepValue(PARTITION_INDEX) = rank;
    		      }

      		  std::vector<int> interface_indices(mNumberOfPartitions, -1);
    		  vector<int>& neighbours_indices = mrModelPart.GetCommunicator().NeighbourIndices();

     		  for(SizeType i = 0 ; i <  neighbours_indices.size() ; i++)
     		    if(SizeType(neighbours_indices[i]) < interface_indices.size())
     		     interface_indices[neighbours_indices[i]] = i;



    		  mLogFile << rank << " : Adding interface nodes to modelpart" << std::endl;
    		  // now adding interface nodes which belongs to other partitions
     		  for(SizeType i_element = 0 ; i_element < NumberOfElements ; i_element++)
     		    if(EPart[i_element] == rank)
    		    {
     		      for(std::vector<std::size_t>::iterator i_node = ElementsConnectivities[i_element].begin() ;
     			  i_node != ElementsConnectivities[i_element].end() ; i_node++)
    			{
    			  int node_partition = NPart[*i_node-1];
    			  if(node_partition != rank)
    			  {
    			    mrModelPart.AssignNode(AllNodes((*i_node)));
    			    mrModelPart.GetCommunicator().GhostMesh().AddNode(AllNodes((*i_node)));
    			    mrModelPart.GetCommunicator().GhostMesh(interface_indices[node_partition]).AddNode(AllNodes((*i_node)));
    			    mrModelPart.GetCommunicator().InterfaceMesh().AddNode(AllNodes((*i_node)));
    			    mrModelPart.GetCommunicator().InterfaceMesh(interface_indices[node_partition]).AddNode(AllNodes((*i_node)));
    //   			   mrModelPart.AssignNode(AllNodes((*i_node)),  ModelPart::Kratos_Ghost);
    			    AllNodes((*i_node))->GetSolutionStepValue(PARTITION_INDEX) = NPart[*i_node-1];
    			  }
    			}
    		    }
    		    else  // adding the owened interface nodes
    		    {
      		      for(std::vector<std::size_t>::iterator i_node = ElementsConnectivities[i_element].begin() ;
      			  i_node != ElementsConnectivities[i_element].end() ; i_node++)
      			if(NPart[*i_node-1] == rank)
     			  {
    			    SizeType mesh_index = interface_indices[EPart[i_element]];
    			    if(mesh_index > mNumberOfPartitions) // Means the neighbour domain is not registered!!
    			      KRATOS_THROW_ERROR(std::logic_error, "Cannot find the neighbour domain : ", EPart[i_element]);

    			    mrModelPart.GetCommunicator().LocalMesh().AddNode(AllNodes((*i_node)));
    			    mrModelPart.GetCommunicator().LocalMesh(mesh_index).AddNode(AllNodes((*i_node)));
    			    mrModelPart.GetCommunicator().InterfaceMesh(mesh_index).AddNode(AllNodes((*i_node)));
    			    mrModelPart.GetCommunicator().InterfaceMesh().AddNode(AllNodes((*i_node)));

    // 			    SizeType mesh_index = interface_indices[EPart[i_element]] +  ModelPart::Kratos_Ownership_Size;
    // 			    if(mesh_index > mNumberOfPartitions  +  ModelPart::Kratos_Ownership_Size) // Means the neighbour domain is not registered!!
    // 			      KRATOS_THROW_ERROR(std::logic_error, "Cannot find the neighbour domain : ", EPart[i_element]);
    //  			    mrModelPart.AssignNode(AllNodes((*i_node)), mesh_index);
    // // 		      mLogFile << rank << " : Adding interface node # " << *i_node << std::endl;
    // 			    mLogFile << rank << " : Adding interface node # " << *i_node << " to mesh: " << mesh_index  << std::endl;
    // 			AllNodes((*i_node))->GetSolutionStepValue(PARTITION_INDEX) = rank;
    // // 			mLogFile << rank << " : Adding interface node # " << AllNodes[(*i_node)] << "with partition index " << AllNodes[(*i_node)].GetSolutionStepValue(PARTITION_INDEX)  << " in " <<  AllNodes((*i_node)) << std::endl;
     			  }
    		    }
    		  mLogFile << rank << " : Nodes added to modelpart" << std::endl;
    		}


    	        void AddingElements(ModelPart::NodesContainerType& AllNodes, idxtype* NPart, idxtype* EPart)
    	        {
    		  int rank = GetRank();

    		  mLogFile << rank << " : Reading elements" << std::endl;
    		  ModelPart::ElementsContainerType temp_elements;
    		  mrIO.ReadElements(AllNodes, mrModelPart.rProperties(), temp_elements);

    		  mLogFile << rank << " : Adding elements to modelpart" << std::endl;
    		  idxtype* epart_position = EPart;
    		  for(Kratos::ModelPart::ElementIterator i_element = temp_elements.begin() ;
    		      i_element != temp_elements.end() ; i_element++)
    		    {
    		      if(*epart_position == rank)
    			{
    			  mrModelPart.AddElement(*(i_element.base()));
     			  mrModelPart.GetCommunicator().LocalMesh().AddElement(*(i_element.base()));
    //   			  mrModelPart.AddElement(*(i_element.base()),  ModelPart::Kratos_Local);
    			}
    		      epart_position++;
    		    }
    		  mLogFile << rank << " : Elements added" << std::endl;
    		}

    	  void AddingConditions(ModelPart::NodesContainerType& AllNodes, idxtype* NPart, idxtype* EPart, ModelPart::ConditionsContainerType& AllConditions)
    	        {
    		  int rank = GetRank();

    		  mLogFile << rank << " : Reading conditions" << std::endl;
    		  mrIO.ReadConditions(AllNodes, mrModelPart.rProperties(), AllConditions);

    		  mLogFile << rank << " : Adding conditions to modelpart" << std::endl;
    		  for(Kratos::ModelPart::ConditionIterator i_condition = AllConditions.begin() ;
    		      i_condition != AllConditions.end() ; i_condition++)
    		    {
    		      bool is_local = 1;
    		      // See if all of the condition nodes are in this partition as a local or even as a ghost
    		      // TODO: THIS IS DANGEROUSE AND MAY FAILE DUE TO THE MESH!!! MUST BE CHANGED!!
    		      for(ModelPart::ConditionType::GeometryType::iterator i_node = i_condition->GetGeometry().begin() ;
    			  i_node != i_condition->GetGeometry().end() ; i_node++)
    			if(mrModelPart.Nodes().find(i_node->Id()) == mrModelPart.Nodes().end())
    			  is_local = 0;
    		      if(is_local)
    			{
    			  mrModelPart.AddCondition(*(i_condition.base()));
     			  mrModelPart.GetCommunicator().LocalMesh().AddCondition(*(i_condition.base()));
    //   			  mrModelPart.AddCondition(*(i_condition.base()),  ModelPart::Kratos_Local);
    			}
    		    }
    		  mLogFile << rank << " : Conditions added" << std::endl;
    		}

    */
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


//		ModelPart& mrModelPart;

//	        IO& mrIO;

//		SizeType mNumberOfPartitions;

    std::vector<int> mContactNodes;

//		std::ofstream mLogFile;

//		SizeType mDimension;


    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{


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
    MetisContactPartitioningProcess& operator=(MetisContactPartitioningProcess const& rOther);

    /// Copy constructor.
    //MetisContactPartitioningProcess(MetisContactPartitioningProcess const& rOther);


    ///@}

}; // Class MetisContactPartitioningProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  MetisContactPartitioningProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const MetisContactPartitioningProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_METIS_CONTACTmContactNodes[i]_PARTITIONING_PROCESS_INCLUDED defined


