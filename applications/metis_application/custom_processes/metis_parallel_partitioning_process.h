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
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2008-04-21 13:33:42 $
//   Revision:            $Revision: 1.1 $
//
//


#if !defined(KRATOS_METIS_PARALLEL_PARTITIONING_PROCESS_INCLUDED )
#define  KRATOS_METIS_PARALLEL_PARTITIONING_PROCESS_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <algorithm>

// External includes 
#include <parmetislib.h>


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"





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

	class MetisParallelPartitioningProcess 
		: public Process
	{
	public:
		///@name Type Definitions
		///@{

		/// Pointer definition of MetisParallelPartitioningProcess
		KRATOS_CLASS_POINTER_DEFINITION(MetisParallelPartitioningProcess);

		typedef std::size_t SizeType;
		typedef std::size_t IndexType;

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
		MetisParallelPartitioningProcess(ModelPart& model_part, SizeType NumberOfPartitions = 0)
		    : mrModelPart(model_part), mNumberOfPartitions(NumberOfPartitions)
		{
			KRATOS_TRY
			KRATOS_CATCH("")
		}

		/// Destructor.
		virtual ~MetisParallelPartitioningProcess()
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
			
			if (mNumberOfPartitions == 0)
			  MPI_Comm_size (MPI_COMM_WORLD,&mNumberOfPartitions);

			const IndexType number_of_nodes_per_element = mrModelPart.ElementsBegin()->GetGeometry().size();
			int i = 0;
			
			KRATOS_WATCH(number_of_nodes_per_element);

			int ne = mrModelPart.NumberOfElements();
			int nn = mrModelPart.NumberOfNodes();
			idxtype* elmnts = new idxtype[mrModelPart.NumberOfElements() * number_of_nodes_per_element];
			int etype = 1;
			int numflag = 0;
			int number_of_partitions = static_cast<int>(mNumberOfPartitions);
			int edgecut;
			idxtype* epart = new idxtype[mrModelPart.NumberOfElements()];
			idxtype* npart = new idxtype[mrModelPart.NumberOfNodes()];
			
			std::cout << "Preparing Data for metis..." << std::endl;
			// Creating the elmnts array for Metis
			for(Kratos::ModelPart::ElementIterator i_element = mrModelPart.ElementsBegin() ;
			    i_element != mrModelPart.ElementsEnd() ; i_element++)
			    for(int j = 0 ; j < number_of_nodes_per_element ; j++)
				elmnts[i++] = i_element->GetGeometry()[j].Id()-1;
			
			std::cout << "Calling metis..." << std::endl;
			// Calling Metis to partition
			ParMETIS_V3_PartMeshKway(&ne, &nn, elmnts, &etype, &numflag, &number_of_partitions, &edgecut, epart, npart,MPI_COMM_WORLD);
			std::cout << "Metis Finished!!!" << std::endl;
			
			
			// Updating the calculated partition index for nodes
/* 			for(IndexType h = 0 ; h < nn ; h++) */
/* 			    model_part.Nodes()[h+1].GetValue(PARTITION_INDEX) = npart[h]; */

			
/* 			for(Kratos::ModelPart::NodeIterator i_node = mrModelPart.NodesBegin() ; */
/* 			    i_node != mrModelPart.NodesEnd() ; i_node++) */
/* 			    model_part.AddNode(npart */

			std::cout << "Resizing modelpart meshes" << std::endl;
			for(int i = mrModelPart.GetMeshes().size() ; i < number_of_partitions + 2  ; i++)
			    mrModelPart.GetMeshes().push_back(ModelPart::MeshType());


			std::cout << "Adding nodes to each partition mesh" << std::endl;
			// Adding node to each partition mesh
 			for(IndexType h = 0 ; h < nn ; h++) 
			    mrModelPart.AddNode(mrModelPart.Nodes()(h+1),npart[h]+1);

			std::cout << "Adding elements to each partition mesh" << std::endl;
			// Adding elements to each partition mesh
			idxtype* elmnt_position = elmnts;
			idxtype* epart_position = epart;
			for(Kratos::ModelPart::ElementIterator i_element = mrModelPart.ElementsBegin() ;
			    i_element != mrModelPart.ElementsEnd() ; i_element++)
			{
				mrModelPart.AddElement(*(i_element.base()), *epart_position + 1);
			    epart_position++;
			}

				


			delete[] elmnts;
			delete[] epart;
			delete[] npart;

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
			return "MetisParallelPartitioningProcess";
		}

		/// Print information about this object.
		virtual void PrintInfo(std::ostream& rOStream) const
		{
			rOStream << "MetisParallelPartitioningProcess";
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
		ModelPart& mrModelPart;
		
		SizeType mNumberOfPartitions;


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
		MetisParallelPartitioningProcess& operator=(MetisParallelPartitioningProcess const& rOther);

		/// Copy constructor.
		//MetisParallelPartitioningProcess(MetisParallelPartitioningProcess const& rOther);


		///@}    

	}; // Class MetisParallelPartitioningProcess 

	///@} 

	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 


	/// input stream function
	inline std::istream& operator >> (std::istream& rIStream, 
		MetisParallelPartitioningProcess& rThis);

	/// output stream function
	inline std::ostream& operator << (std::ostream& rOStream, 
		const MetisParallelPartitioningProcess& rThis)
	{
		rThis.PrintInfo(rOStream);
		rOStream << std::endl;
		rThis.PrintData(rOStream);

		return rOStream;
	}
	///@} 


}  // namespace Kratos.

#endif // KRATOS_METIS_PARALLEL_PARTITIONING_PROCESS_INCLUDED defined 


