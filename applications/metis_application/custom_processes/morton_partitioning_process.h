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


#if !defined(KRATOS_MORTON_PARTITIONING_PROCESS_INCLUDED )
#define  KRATOS_MORTON_PARTITIONING_PROCESS_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <list>

// Project includes
#include "includes/define.h"
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

class MortonPartitioningProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    #ifdef KRATOS_USE_METIS_5
      typedef idx_t idxtype;
    #endif

    /// Pointer definition of MortonPartitioningProcess
    KRATOS_CLASS_POINTER_DEFINITION(MortonPartitioningProcess);

    typedef std::size_t SizeType;
    typedef std::size_t IndexType;
    typedef std::vector<idxtype> PartitionIndicesType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MortonPartitioningProcess(IO::ConnectivitiesContainerType& rElementsConnectivities,
                              IO::NodesContainerType& rNodesContainer,
                                  PartitionIndicesType& rNodesPartitions,
                                  PartitionIndicesType& rElementsPartitions,
                                  SizeType NumberOfPartitions, int Dimension = 3)
        : mrElementsConnectivities(rElementsConnectivities),
          mrNodesContainer(rNodesContainer),
          mrNodesPartitions(rNodesPartitions),
          mrElementsPartitions(rElementsPartitions),
          mNumberOfPartitions(NumberOfPartitions),
          mDimension(Dimension)
    {
    }

    /// Copy constructor.
    MortonPartitioningProcess(MortonPartitioningProcess const& rOther)
        : mrElementsConnectivities(rOther.mrElementsConnectivities),
          mrNodesContainer(rOther.mrNodesContainer),
          mrNodesPartitions(rOther.mrNodesPartitions),
          mrElementsPartitions(rOther.mrElementsPartitions),
          mNumberOfPartitions(rOther.mNumberOfPartitions),
          mDimension(rOther.mDimension)
    {
    }

    /// Destructor.
    virtual ~MortonPartitioningProcess()
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
        int real_r_of_nodes = 0;
        // calculating number of nodes considering sequencial numbering!! We can get it from input for no sequencial one. Pooyan.
        for(IO::ConnectivitiesContainerType::iterator i_element = mrElementsConnectivities.begin() ; i_element != mrElementsConnectivities.end() ; i_element++)
            for(IO::ConnectivitiesContainerType::value_type::iterator i_node_id = i_element->begin() ; i_node_id != i_element->end() ; i_node_id++)
            {
                if(static_cast<int>(*i_node_id) > number_of_nodes)
                    number_of_nodes = *i_node_id;

                real_r_of_nodes++;
            }

        //verify that all of the nodes exist
        std::vector< bool > aux(number_of_nodes,false);
        // calculating number of nodes considering sequencial numbering!! We can get it from input for no sequencial one. Pooyan.
        for(IO::ConnectivitiesContainerType::iterator i_element = mrElementsConnectivities.begin() ; i_element != mrElementsConnectivities.end() ; i_element++)
            for(IO::ConnectivitiesContainerType::value_type::iterator i_node_id = i_element->begin() ; i_node_id != i_element->end() ; i_node_id++)
            {
                aux[static_cast<int>(*i_node_id)-1] = true;
            }

        mrElementsPartitions.resize(number_of_elements);
        mrNodesPartitions.resize(number_of_nodes);

        idxtype* epart = &(*(mrElementsPartitions.begin()));
        idxtype* npart = &(*(mrNodesPartitions.begin()));

        AssignPartition(number_of_nodes, real_r_of_nodes, npart, epart);

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
        return "MortonPartitioningProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MortonPartitioningProcess";
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
    IO::NodesContainerType& mrNodesContainer;
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

    void AssignPartition(SizeType NumberOfNodes, SizeType NumberOfElements, idxtype* NPart, idxtype* EPart)
    {
        for(SizeType i_element = 0 ; i_element < NumberOfElements ; i_element++)
        {
            for(IO::ConnectivitiesContainerType::value_type::iterator i_node = mrElementsConnectivities[i_element].begin() ;
                    i_node != mrElementsConnectivities[i_element].end() ; i_node++)
            {
                //IO::NodesContainerType::iterator i_nod = mrNodesContainer.begin() + (*i_node-1);

                NPart[(*i_node-1)] = i_element%mNumberOfPartitions;
                EPart[i_element] = i_element%mNumberOfPartitions;
            }
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
    MortonPartitioningProcess& operator=(MortonPartitioningProcess const& rOther);

    /// Copy constructor.
    //MortonGraphPartitioningProcess(MortonPartitioningProcess const& rOther);


    ///@}

}; // Class MortonGraphPartitioningProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  MortonPartitioningProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const MortonPartitioningProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_MORTON_PARTITIONING_PROCESS_INCLUDED defined


