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


#if !defined(KRATOS_METIS_GRAPH_PARTITIONING_PROCESS_INCLUDED )
#define  KRATOS_METIS_GRAPH_PARTITIONING_PROCESS_INCLUDED



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

 extern "C" { 
   //extern void METIS_PartMeshDual(int*, int*, idxtype*, int*, int*, int*, int*, idxtype*, idxtype*); 
 extern int METIS_PartMeshDual(int*, int*, int*, int*, int*, int*, int*, int*, int*); 
 }; 



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

	class MetisGraphPartitioningProcess 
		: public Process
	{
	public:
		///@name Type Definitions
		///@{

		/// Pointer definition of MetisGraphPartitioningProcess
		KRATOS_CLASS_POINTER_DEFINITION(MetisGraphPartitioningProcess);

		typedef std::size_t SizeType;
		typedef std::size_t IndexType;
		typedef std::vector<idxtype> PartitionIndicesType;

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
		MetisGraphPartitioningProcess(IO::ConnectivitiesContainerType& rElementsConnectivities, 
									  PartitionIndicesType& rNodesPartitions, 
									  PartitionIndicesType& rElementsPartitions, 
									  SizeType NumberOfPartitions, int Dimension = 3)
			: mrElementsConnectivities(rElementsConnectivities),
			  mrNodesPartitions(rNodesPartitions), 
			  mrElementsPartitions(rElementsPartitions), 
			  mNumberOfPartitions(NumberOfPartitions), 
			  mDimension(Dimension)
		{
		}

		/// Copy constructor.
		MetisGraphPartitioningProcess(MetisGraphPartitioningProcess const& rOther)
			: mrElementsConnectivities(rOther.mrElementsConnectivities),
			  mrNodesPartitions(rOther.mrNodesPartitions), 
			  mrElementsPartitions(rOther.mrElementsPartitions), 
			  mNumberOfPartitions(rOther.mNumberOfPartitions), 
			  mDimension(rOther.mDimension)
		{
		}

		/// Destructor.
		virtual ~MetisGraphPartitioningProcess()
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


			int number_of_elements = mrElementsConnectivities.size();

			int number_of_nodes = 0; 
			// calculating number of nodes considering sequencial numbering!! We can get it from input for no sequencial one. Pooyan.
			for(IO::ConnectivitiesContainerType::iterator i_element = mrElementsConnectivities.begin() ; i_element != mrElementsConnectivities.end() ; i_element++)
				for(IO::ConnectivitiesContainerType::value_type::iterator i_node_id = i_element->begin() ; i_node_id != i_element->end() ; i_node_id++)
					if(static_cast<int>(*i_node_id) > number_of_nodes)
						number_of_nodes = *i_node_id;

			mrElementsPartitions.resize(number_of_elements);
			mrNodesPartitions.resize(number_of_nodes);


			idxtype* epart = &(*(mrElementsPartitions.begin()));
			idxtype* npart = &(*(mrNodesPartitions.begin()));

			CallingMetis(number_of_nodes, number_of_elements, mrElementsConnectivities, npart, epart);

			KRATOS_CATCH("")
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
			return "MetisGraphPartitioningProcess";
		}

		/// Print information about this object.
		virtual void PrintInfo(std::ostream& rOStream) const
		{
			rOStream << "MetisGraphPartitioningProcess";
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
		///@name Protected static Member Variables 
		///@{ 


		///@} 
		///@name Protected member Variables 
		///@{ 

		IO::ConnectivitiesContainerType& mrElementsConnectivities;
		PartitionIndicesType& mrNodesPartitions; 
		PartitionIndicesType& mrElementsPartitions; 		
		SizeType mNumberOfPartitions;

		SizeType mDimension;


		///@} 
		///@name Protected Operators
		///@{ 


		///@} 
		///@name Protected Operations
		///@{ 


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
		MetisGraphPartitioningProcess& operator=(MetisGraphPartitioningProcess const& rOther);

		/// Copy constructor.
		//MetisGraphPartitioningProcess(MetisGraphPartitioningProcess const& rOther);


		///@}    

	}; // Class MetisGraphPartitioningProcess 

	///@} 

	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 


	/// input stream function
	inline std::istream& operator >> (std::istream& rIStream, 
		MetisGraphPartitioningProcess& rThis);

	/// output stream function
	inline std::ostream& operator << (std::ostream& rOStream, 
		const MetisGraphPartitioningProcess& rThis)
	{
		rThis.PrintInfo(rOStream);
		rOStream << std::endl;
		rThis.PrintData(rOStream);

		return rOStream;
	}
	///@} 


}  // namespace Kratos.

#endif // KRATOS_METIS_GRAPH_PARTITIONING_PROCESS_INCLUDED defined 


