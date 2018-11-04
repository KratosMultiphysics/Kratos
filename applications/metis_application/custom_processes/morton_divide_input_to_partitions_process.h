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

#if !defined(KRATOS_MORTON_DIVIDE_INPUT_TO_PARTITIONS_PROCESS_INCLUDED )
#define  KRATOS_MORTON_DIVIDE_INPUT_TO_PARTITIONS_PROCESS_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <algorithm>
#include <fstream>

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "processes/graph_coloring_process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "custom_processes/morton_partitioning_process.h"

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

class MortonDivideInputToPartitionsProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    #ifdef KRATOS_USE_METIS_5
      typedef idx_t idxtype;
    #endif

    /// Pointer definition of MortonDivideInputToPartitionsProcess
    KRATOS_CLASS_POINTER_DEFINITION(MortonDivideInputToPartitionsProcess);

    typedef std::size_t SizeType;
    typedef std::size_t IndexType;
    typedef GraphColoringProcess::GraphType GraphType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MortonDivideInputToPartitionsProcess(IO& rIO, SizeType NumberOfPartitions, int Dimension = 3)
        :  mrIO(rIO), mNumberOfPartitions(NumberOfPartitions), mDimension(Dimension)
    {
    }

    /// Copy constructor.
    MortonDivideInputToPartitionsProcess(MortonDivideInputToPartitionsProcess const& rOther)
        : mrIO(rOther.mrIO), mNumberOfPartitions(rOther.mNumberOfPartitions), mDimension(rOther.mDimension)
    {
    }

    /// Destructor.
    virtual ~MortonDivideInputToPartitionsProcess()
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

        if(mNumberOfPartitions < 2) // There is no need to partition it and just reading the input
        {
            return;
        }

        // Reading connectivities
        IO::ConnectivitiesContainerType elements_connectivities;
        IO::ConnectivitiesContainerType conditions_connectivities;
        IO::NodesContainerType nodes_container;

        int number_of_elements = mrIO.ReadElementsConnectivities(elements_connectivities);

        mrIO.ReadConditionsConnectivities(conditions_connectivities);
        mrIO.ReadNodes(nodes_container);

        MortonPartitioningProcess::PartitionIndicesType nodes_partitions;
        MortonPartitioningProcess::PartitionIndicesType elements_partitions;
        MortonPartitioningProcess::PartitionIndicesType conditions_partitions;

        MortonPartitioningProcess morton_partitioning_process(elements_connectivities, nodes_container, nodes_partitions, elements_partitions, mNumberOfPartitions, mDimension);
        morton_partitioning_process.Execute();

//         GraphType domains_graph = zero_matrix<int>(mNumberOfPartitions, mNumberOfPartitions);
        GraphType domains_graph(mNumberOfPartitions,mNumberOfPartitions);
        domains_graph = ScalarMatrix(mNumberOfPartitions,mNumberOfPartitions,1);

        GraphType domains_colored_graph;

        int colors_number;

        CalculateDomainsGraph(domains_graph, number_of_elements, elements_connectivities, nodes_partitions, elements_partitions);

        GraphColoringProcess(mNumberOfPartitions, domains_graph, domains_colored_graph, colors_number).Execute();

        //KRATOS_WATCH(domains_graph);

// 			std::vector<DomainEntitiesIdContainer> domains_nodes;
        IO::PartitionIndicesContainerType nodes_all_partitions;
        IO::PartitionIndicesContainerType elements_all_partitions;
        IO::PartitionIndicesContainerType conditions_all_partitions;

        ConditionsPartitioning(conditions_connectivities, nodes_partitions, conditions_partitions);
        KRATOS_WATCH("ConditionsPartitioning finished")
        // Dividing nodes
        DividingNodes(nodes_all_partitions, elements_connectivities, conditions_connectivities, nodes_partitions, elements_partitions, conditions_partitions);
        KRATOS_WATCH("DividingNodes finished")
        // Dividing elements
        DividingElements(elements_all_partitions, elements_partitions);
        KRATOS_WATCH("DividingElements finished")
        // Dividing conditions
        DividingConditions(conditions_all_partitions, conditions_partitions);
        KRATOS_WATCH("DividingConditions finished")

        IO::PartitionIndicesType io_nodes_partitions(nodes_partitions.begin(), nodes_partitions.end());
        IO::PartitionIndicesType io_elements_partitions(elements_partitions.begin(), elements_partitions.end());
        IO::PartitionIndicesType io_conditions_partitions(conditions_partitions.begin(), conditions_partitions.end());

        // Now dividing the input file
        mrIO.DivideInputToPartitions(mNumberOfPartitions, domains_colored_graph,
                                     io_nodes_partitions, io_elements_partitions, io_conditions_partitions,
                                     nodes_all_partitions, elements_all_partitions, conditions_all_partitions);
        KRATOS_WATCH("DivideInputToPartitions finished")
        return;

        KRATOS_CATCH("")
    }

    void CalculateDomainsGraph(GraphType& rDomainsGraph, SizeType NumberOfElements, IO::ConnectivitiesContainerType& ElementsConnectivities, MortonPartitioningProcess::PartitionIndicesType const& NPart, MortonPartitioningProcess::PartitionIndicesType const&  EPart )
    {
        for(SizeType i_element = 0 ; i_element < NumberOfElements ; i_element++)
        {
            for(std::vector<std::size_t>::iterator i_node = ElementsConnectivities[i_element].begin() ;
                    i_node != ElementsConnectivities[i_element].end() ; i_node++)
            {
                SizeType node_rank = NPart[*i_node-1];
                SizeType element_rank = EPart[i_element];
                if(node_rank != element_rank)
                {
                    //matrix!
                    rDomainsGraph(node_rank, element_rank) = 1;
                    rDomainsGraph(element_rank, node_rank) = 1;
                }
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
    std::string Info() const override
    {
        return "MortonDivideInputToPartitionsProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MortonDivideInputToPartitionsProcess";
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


    void ConditionsPartitioning(IO::ConnectivitiesContainerType& ConditionsConnectivities,
                                MortonPartitioningProcess::PartitionIndicesType const& NodesPartitions,
                                MortonPartitioningProcess::PartitionIndicesType& ConditionsPartitions)
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
                       MortonPartitioningProcess::PartitionIndicesType const& NodesPartitions,
                       MortonPartitioningProcess::PartitionIndicesType const& ElementsPartitions,
                       MortonPartitioningProcess::PartitionIndicesType const& ConditionsPartitions)
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
                {
                    std::cout << "Partiton element differnet form partition node" << std::endl;
                    rNodesAllPartitions[my_gid].push_back(element_partition);
                }
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

    void DividingElements(IO::PartitionIndicesContainerType& rElementsAllPartitions, MortonPartitioningProcess::PartitionIndicesType const& ElementsPartitions)
    {
        SizeType number_of_elements = ElementsPartitions.size();

        rElementsAllPartitions.resize(number_of_elements);

        // adding the elements partition to their array of partitions
        for(SizeType i_element = 0 ; i_element < number_of_elements ; i_element++)
        {
            rElementsAllPartitions[i_element].push_back(ElementsPartitions[i_element]);
        }
    }

    void DividingConditions(IO::PartitionIndicesContainerType& rConditionsAllPartitions, MortonPartitioningProcess::PartitionIndicesType const& ConditionsPartitions)
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
    MortonDivideInputToPartitionsProcess& operator=(MortonDivideInputToPartitionsProcess const& rOther);

    /// Copy constructor.
    //MortonDivideInputToPartitionsProcess(MortonDivideInputToPartitionsProcess const& rOther);


    ///@}

}; // Class MortonDivideInputToPartitionsProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  MortonDivideInputToPartitionsProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const MortonDivideInputToPartitionsProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_MORTON_DIVIDE_INPUT_TO_PARTITIONS_PROCESS_INCLUDED defined
