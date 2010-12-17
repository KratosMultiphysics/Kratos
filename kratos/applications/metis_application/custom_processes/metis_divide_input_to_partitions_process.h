/*
==============================================================================
KratosPFEMApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


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
 
//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2009-01-15 11:11:35 $
//   Revision:            $Revision: 1.6 $
//
//


#if !defined(KRATOS_METIS_DIVIDE_INPUT_TO_PARTITIONS_PROCESS_INCLUDED )
#define  KRATOS_METIS_DIVIDE_INPUT_TO_PARTITIONS_PROCESS_INCLUDED



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
#include "custom_processes/metis_graph_partitioning_process.h"



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

	class MetisDivideInputToPartitionsProcess 
		: public Process
	{
	public:
		///@name Type Definitions
		///@{

		/// Pointer definition of MetisDivideInputToPartitionsProcess
		KRATOS_CLASS_POINTER_DEFINITION(MetisDivideInputToPartitionsProcess);

		typedef std::size_t SizeType;
		typedef std::size_t IndexType;
	        typedef matrix<int> GraphType;

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
	MetisDivideInputToPartitionsProcess(ModelPart& rModelPart, IO& rIO, SizeType NumberOfPartitions, int Dimension = 3)
	  :  mrIO(rIO), mNumberOfPartitions(NumberOfPartitions), mDimension(Dimension)
		{
		}

		/// Copy constructor.
		MetisDivideInputToPartitionsProcess(MetisDivideInputToPartitionsProcess const& rOther)
		  : mrIO(rOther.mrIO), mNumberOfPartitions(rOther.mNumberOfPartitions), mDimension(rOther.mDimension)
		{
		}

		/// Destructor.
		virtual ~MetisDivideInputToPartitionsProcess()
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

		virtual void Execute()
		{
			KRATOS_TRY;

			if(mNumberOfPartitions < 2) // There is no need to partition it and just reading the input
			{
			  return;
			}
			
			// Reading connectivities
			IO::ConnectivitiesContainerType elements_connectivities;
			IO::ConnectivitiesContainerType conditions_connectivities;
			int number_of_elements =  mrIO.ReadElementsConnectivities(elements_connectivities);
			mrIO.ReadConditionsConnectivities(conditions_connectivities);
			
			MetisGraphPartitioningProcess::PartitionIndicesType nodes_partitions;
			MetisGraphPartitioningProcess::PartitionIndicesType elements_partitions;
			MetisGraphPartitioningProcess::PartitionIndicesType conditions_partitions;
			MetisGraphPartitioningProcess metis_graph_partitioning_process(elements_connectivities, nodes_partitions, elements_partitions, mNumberOfPartitions, mDimension);
			metis_graph_partitioning_process.Execute();
			
			int number_of_nodes = nodes_partitions.size(); 

			GraphType domains_graph = zero_matrix<int>(mNumberOfPartitions, mNumberOfPartitions);
			GraphType domains_colored_graph;

			int colors_number;
			
			CalculateDomainsGraph(domains_graph, number_of_elements, elements_connectivities, nodes_partitions, elements_partitions);
			GraphColoringProcess(mNumberOfPartitions, domains_graph,domains_colored_graph, colors_number).Execute();
			// 		      colors_number = GraphColoring(domains_graph, domains_colored_graph);
			KRATOS_WATCH(colors_number);
			KRATOS_WATCH(domains_colored_graph);

// 			std::vector<DomainEntitiesIdContainer> domains_nodes;
			IO::PartitionIndicesContainerType nodes_all_partitions;
			IO::PartitionIndicesContainerType elements_all_partitions;
			IO::PartitionIndicesContainerType conditions_all_partitions;
			
			ConditionsPartitioning(conditions_connectivities, nodes_partitions, conditions_partitions);

			// Dividing nodes 
 			DividingNodes(nodes_all_partitions, elements_connectivities, conditions_connectivities, nodes_partitions, elements_partitions, conditions_partitions);

			// Dividing elements
 			DividingElements(elements_all_partitions, elements_partitions);

			// Dividing conditions
 			DividingConditions(conditions_all_partitions, conditions_partitions);

			
/*			for(SizeType i = 0 ; i < number_of_nodes ; i++)
			{
			  std::cout << "Node #" << i << "->";
			  for(std::vector<std::size_t>::iterator j = nodes_all_partitions[i].begin() ; j != nodes_all_partitions[i].end() ; j++)
			    std::cout << *j << ",";
			  std::cout << std::endl;
			}*/

			IO::PartitionIndicesType io_nodes_partitions(nodes_partitions.begin(), nodes_partitions.end());
			IO::PartitionIndicesType io_elements_partitions(elements_partitions.begin(), elements_partitions.end());
			IO::PartitionIndicesType io_conditions_partitions(conditions_partitions.begin(), conditions_partitions.end());
			
			// Now dividing the input file
			mrIO.DivideInputToPartitions(mNumberOfPartitions, domains_colored_graph,
						     io_nodes_partitions, io_elements_partitions, io_conditions_partitions,
						     nodes_all_partitions, elements_all_partitions, conditions_all_partitions);

			return;
			
			
			
			
			
			
			
			
			// Adding properties to modelpart
// 			mrIO.ReadProperties(mrModelPart.rProperties());

			// Adding elements to each partition mesh
// 			AddingElements(temp_nodes, npart, epart);


			// Adding conditions to each partition mesh
// 			ModelPart::ConditionsContainerType temp_conditions;
// 			AddingConditions(temp_nodes, npart, epart, temp_conditions);



// 			mrIO.ReadInitialValues(temp_nodes, mrModelPart.Elements(), temp_conditions);

			KRATOS_CATCH("")
		}

		void CalculateDomainsGraph(GraphType& rDomainsGraph, SizeType NumberOfElements, IO::ConnectivitiesContainerType& ElementsConnectivities, MetisGraphPartitioningProcess::PartitionIndicesType const& NPart, MetisGraphPartitioningProcess::PartitionIndicesType const&  EPart )
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
			return "MetisDivideInputToPartitionsProcess";
		}

		/// Print information about this object.
		virtual void PrintInfo(std::ostream& rOStream) const
		{
			rOStream << "MetisDivideInputToPartitionsProcess";
		}

		/// Print object's data.
		virtual void PrintData(std::ostream& rOStream) const
		{
		}


		///@}      
		///@name Friends
		///@{


		///@}

	protected:

		class DomainEntitiesIdContainer
		{
		public:
			DomainEntitiesIdContainer(std::size_t NumberOfNeighbours)
			{
				mLocalIds.resize(NumberOfNeighbours);
				mGhostsIds.resize(NumberOfNeighbours);
				mInterfacesIds.resize(NumberOfNeighbours);
			}

			std::vector<std::size_t>& AllIds()
			{
				return mAllIds;
			}

			std::vector<std::vector<std::size_t> >& LocalIds()
			{
				return mLocalIds;
			}

			std::vector<std::vector<std::size_t> >& GhostsIds()
			{
				return mGhostsIds;
			}

			std::vector<std::vector<std::size_t> >& InterfacesIds()
			{
				return mInterfacesIds;
			}
		private:
			std::vector<std::size_t> mAllIds;
			std::vector<std::vector<std::size_t> > mLocalIds;
			std::vector<std::vector<std::size_t> > mGhostsIds;
			std::vector<std::vector<std::size_t> > mInterfacesIds;

		};

		///@name Protected static Member Variables 
		///@{ 


		///@} 
		///@name Protected member Variables 
		///@{ 

// 		ModelPart& mrModelPart;

	        IO& mrIO;
		
		SizeType mNumberOfPartitions;

		SizeType mDimension;


		///@} 
		///@name Protected Operators
		///@{ 


		///@} 
		///@name Protected Operations
		///@{ 

/*
		void CallingMetis(SizeType NumberOfNodes, SizeType NumberOfElements, IO::ConnectivitiesContainerType& ElementsConnectivities, idxtype* NPart, idxtype* EPart)
		{
			// calculating total size of connectivity vector 
			int connectivity_size = 0;
			for(IO::ConnectivitiesContainerType::iterator i_connectivities = ElementsConnectivities.begin() ;
				i_connectivities != ElementsConnectivities.end() ; i_connectivities++)
				connectivity_size += i_connectivities->size();

			int number_of_element_nodes = ElementsConnectivities.begin()->size(); // here assuming that all elements are the same!!

			int ne = NumberOfElements;
			int nn = NumberOfNodes;



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
			else if(number_of_element_nodes == 8) // hexahedra
				etype = 3;
			else
				KRATOS_ERROR(std::invalid_argument, "invalid element type with number of nodes : ", number_of_element_nodes);


			int numflag = 0;
			int number_of_partitions = static_cast<int>(mNumberOfPartitions);
			int edgecut;

			idxtype* elmnts = new idxtype[connectivity_size];

			int i = 0;
			// Creating the elmnts array for Metis
			for(IO::ConnectivitiesContainerType::iterator i_connectivities = ElementsConnectivities.begin() ; 
				i_connectivities != ElementsConnectivities.end() ; i_connectivities++)
				for(unsigned int j = 0 ; j < i_connectivities->size() ; j++)
					elmnts[i++] = (*i_connectivities)[j] - 1; // transforming to zero base indexing

			// Calling Metis to partition
			METIS_PartMeshDual(&ne, &nn, elmnts, &etype, &numflag, &number_of_partitions, &edgecut, EPart, NPart);

			delete[] elmnts;

		}*/
			
		void ConditionsPartitioning(IO::ConnectivitiesContainerType& ConditionsConnectivities, 
					    MetisGraphPartitioningProcess::PartitionIndicesType const& NodesPartitions, 
					    MetisGraphPartitioningProcess::PartitionIndicesType& ConditionsPartitions)
		{
		  SizeType number_of_conditions = ConditionsConnectivities.size();
		  
		  ConditionsPartitions.resize(number_of_conditions);
		  
		  // getting the average of the partion indices of the condition nodes and take the nearest partition indices of the nodes to the average.
		  for(SizeType i_condition = 0 ; i_condition < number_of_conditions ; i_condition++)
		  {
		    if(ConditionsConnectivities[i_condition].size() > 0)
		    {
		      double average_index = 0.00;
		      for(IO::ConnectivitiesContainerType::value_type::iterator i_node = ConditionsConnectivities[i_condition].begin() ;  
		      i_node != ConditionsConnectivities[i_condition].end() ; i_node++)
		      {
			//get global id. We assume that node ids are began with one
			const int my_gid = *i_node-1;
			
			average_index += NodesPartitions[my_gid];
			
		      }
		      average_index /= ConditionsConnectivities[i_condition].size();
		      
		      double difference = mNumberOfPartitions + 10; // Just to be sure! ;-)
		      
		      for(IO::ConnectivitiesContainerType::value_type::iterator i_node = ConditionsConnectivities[i_condition].begin() ;  
		      i_node != ConditionsConnectivities[i_condition].end() ; i_node++)
		      {
			//get global id. We assume that node ids are began with one
			const int my_gid = *i_node-1;
			
			//get the partition index for the node i am interested in
			const int node_partition = NodesPartitions[my_gid];
			
			if(difference > fabs(average_index - node_partition))
			{
			  difference = fabs(average_index - node_partition);
			  ConditionsPartitions[i_condition] = node_partition;	
			}
		      }
		      
		    }
		  }
		}
		
 		void DividingNodes(IO::PartitionIndicesContainerType& rNodesAllPartitions, 
				   IO::ConnectivitiesContainerType& ElementsConnectivities, 
				   IO::ConnectivitiesContainerType& ConditionsConnectivities, 
				   MetisGraphPartitioningProcess::PartitionIndicesType const& NodesPartitions, 
				   MetisGraphPartitioningProcess::PartitionIndicesType const& ElementsPartitions, 
				   MetisGraphPartitioningProcess::PartitionIndicesType const& ConditionsPartitions)
 		{
		  SizeType number_of_nodes = NodesPartitions.size();
		  SizeType number_of_elements = ElementsPartitions.size();
		  SizeType number_of_conditions = ConditionsPartitions.size();
 		  
		  rNodesAllPartitions.resize(number_of_nodes);
		  
		  for(SizeType i_element = 0 ; i_element < number_of_elements ; i_element++) 
		  {
		    const int element_partition = ElementsPartitions[i_element];
		    
		    //for each element in the model loop over its connectivities
		    for(IO::ConnectivitiesContainerType::value_type::iterator i_node = ElementsConnectivities[i_element].begin() ;  
		    i_node != ElementsConnectivities[i_element].end() ; i_node++)
		    {
		      //get global id. We assume that node ids are began with one
		      const int my_gid = *i_node-1;
		      
		      //get the partition index for the node i am interested in
		      const int node_partition = NodesPartitions[my_gid];

		      // adding the partition of the element to its nodes
		      if(element_partition != node_partition) // we will add the node_partition once afterward
			 rNodesAllPartitions[my_gid].push_back(element_partition);
		    }
		  }
		  
		  for(SizeType i_condition = 0 ; i_condition < number_of_conditions ; i_condition++) 
		  {
		    const int condition_partition = ConditionsPartitions[i_condition];
		    
		    //for each element in the model loop over its connectivities
		    for(IO::ConnectivitiesContainerType::value_type::iterator i_node = ConditionsConnectivities[i_condition].begin() ;  
		    i_node != ConditionsConnectivities[i_condition].end() ; i_node++)
		    {
		      //get global id. We assume that node ids are began with one
		      const int my_gid = *i_node-1;
		      
		      //get the partition index for the node i am interested in
		      const int node_partition = NodesPartitions[my_gid];

		      // adding the partition of the element to its nodes
		      if(condition_partition != node_partition) // we will add the node_partition once afterward
			 rNodesAllPartitions[my_gid].push_back(condition_partition);
		    }
		  }
		  
		  // adding the nodes partition to their array of partitions and clear the repeated ones
		  for(SizeType i_node = 0 ; i_node < number_of_nodes ; i_node++)
		  {
		    IO::PartitionIndicesContainerType::value_type& node_partitions = rNodesAllPartitions[i_node];
		    node_partitions.push_back(NodesPartitions[i_node]);
	    
		    std::sort(node_partitions.begin(), node_partitions.end());
		    IO::PartitionIndicesContainerType::value_type::iterator new_end=std::unique(node_partitions.begin(), node_partitions.end());
		    node_partitions.resize(new_end - node_partitions.begin());
		  }
		}
		
 		void DividingElements(IO::PartitionIndicesContainerType& rElementsAllPartitions, MetisGraphPartitioningProcess::PartitionIndicesType const& ElementsPartitions)
 		{
		  SizeType number_of_elements = ElementsPartitions.size();
 		  
		  rElementsAllPartitions.resize(number_of_elements);
		  
		  // adding the elements partition to their array of partitions 
		  for(SizeType i_element = 0 ; i_element < number_of_elements ; i_element++)
		  {
		    rElementsAllPartitions[i_element].push_back(ElementsPartitions[i_element]);
		  }
		}
		
		void DividingConditions(IO::PartitionIndicesContainerType& rConditionsAllPartitions, MetisGraphPartitioningProcess::PartitionIndicesType const& ConditionsPartitions)
		{
		  SizeType number_of_conditions = ConditionsPartitions.size();
		  
		  rConditionsAllPartitions.resize(number_of_conditions);
		  
		  // adding the condition partition to their array of partitions 
		  for(SizeType i_condition = 0 ; i_condition < number_of_conditions ; i_condition++)
		  {
		    rConditionsAllPartitions[i_condition].push_back(ConditionsPartitions[i_condition]);
		  }
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
		MetisDivideInputToPartitionsProcess& operator=(MetisDivideInputToPartitionsProcess const& rOther);

		/// Copy constructor.
		//MetisDivideInputToPartitionsProcess(MetisDivideInputToPartitionsProcess const& rOther);


		///@}    

	}; // Class MetisDivideInputToPartitionsProcess 

	///@} 

	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 


	/// input stream function
	inline std::istream& operator >> (std::istream& rIStream, 
		MetisDivideInputToPartitionsProcess& rThis);

	/// output stream function
	inline std::ostream& operator << (std::ostream& rOStream, 
		const MetisDivideInputToPartitionsProcess& rThis)
	{
		rThis.PrintInfo(rOStream);
		rOStream << std::endl;
		rThis.PrintData(rOStream);

		return rOStream;
	}
	///@} 


}  // namespace Kratos.

#endif // KRATOS_METIS_DIVIDE_INPUT_TO_PARTITIONS_PROCESS_INCLUDED defined 


