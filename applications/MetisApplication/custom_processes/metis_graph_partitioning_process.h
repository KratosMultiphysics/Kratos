//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//
//

#pragma once

// System includes

// External includes
#include "metis.h"

// Project includes
#include "includes/io.h"
#include "processes/process.h"

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

    using idxtype = idx_t;

    /// Pointer definition of MetisGraphPartitioningProcess
    KRATOS_CLASS_POINTER_DEFINITION(MetisGraphPartitioningProcess);

    typedef std::size_t SizeType;
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

    void Execute() override
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
// KRATOS_WATCH(number_of_nodes)

        //verify that all of the nodes exist
        std::vector< bool > aux(number_of_nodes,false);
        // calculating number of nodes considering sequencial numbering!! We can get it from input for no sequencial one. Pooyan.
        for(IO::ConnectivitiesContainerType::iterator i_element = mrElementsConnectivities.begin() ; i_element != mrElementsConnectivities.end() ; i_element++)
            for(IO::ConnectivitiesContainerType::value_type::iterator i_node_id = i_element->begin() ; i_node_id != i_element->end() ; i_node_id++)
            {
                aux[static_cast<int>(*i_node_id)-1] = true;
            }

        for(unsigned int i=0; i<aux.size(); i++)
            if(aux[i] != true)
            {
                KRATOS_THROW_ERROR(std::logic_error,"Isolated node found! The problematic node has Id  ",i+1);
            }
// KRATOS_WATCH("sequential numbering of nodes verified")


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
    std::string Info() const override
    {
        return "MetisGraphPartitioningProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MetisGraphPartitioningProcess";
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

        if(number_of_element_nodes != 3 && number_of_element_nodes != 4  && number_of_element_nodes != 8) {
            KRATOS_THROW_ERROR(std::invalid_argument, "invalid element type with number of nodes : ", number_of_element_nodes);
        }

        //int numflag = 0;
        //int number_of_partitions = static_cast<int>(mNumberOfPartitions);
        //int edgecut;


        idxtype* elmnts = new idxtype[connectivity_size];

        int i = 0;
        // Creating the elmnts array for Metis
        for(IO::ConnectivitiesContainerType::iterator i_connectivities = ElementsConnectivities.begin() ;
                i_connectivities != ElementsConnectivities.end() ; i_connectivities++)
            for(unsigned int j = 0 ; j < i_connectivities->size() ; j++)
                elmnts[i++] = (*i_connectivities)[j] - 1; // transforming to zero base indexing

        // Calling Metis to partition
        //METIS_PartMeshDual(&ne, &nn, elmnts, &etype, &numflag, &number_of_partitions, &edgecut, EPart, NPart);
        KRATOS_WATCH("not implemented!!!")

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
