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

        //verify that all of the nodes exist
        std::vector< bool > aux(number_of_nodes,false);
        // calculating number of nodes considering sequencial numbering!! We can get it from input for no sequencial one. Pooyan.
        for(IO::ConnectivitiesContainerType::iterator i_element = mrElementsConnectivities.begin() ; i_element != mrElementsConnectivities.end() ; i_element++)
            for(IO::ConnectivitiesContainerType::value_type::iterator i_node_id = i_element->begin() ; i_node_id != i_element->end() ; i_node_id++)
            {
                aux[static_cast<int>(*i_node_id)-1] = true;
            }

        // ? correspondecias 1 a 1 ara aixi que no pot pasar.
//         for(unsigned int i=0; i<aux.size(); i++)
//             if(aux[i] != true)
//             {
//                 KRATOS_ERROR(std::logic_error,"Isolated node found! The problematic node has Id  ",i+1);
//             }

        mrElementsPartitions.resize(number_of_elements);
        mrNodesPartitions.resize(number_of_nodes);

        idxtype* epart = &(*(mrElementsPartitions.begin()));
        idxtype* npart = &(*(mrNodesPartitions.begin()));

        AssignPartition(number_of_nodes, number_of_elements, npart, epart);

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
        return "MortonPartitioningProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "MortonPartitioningProcess";
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
      
    int max (int a, int b)
    {
        return a > b ? a : b;
    }
    
    //////////////////////////////////////////
    // Interleave 2x16 table
      
    // 0101 0101 0101 0101 0101 0101 0101 0101 << 1
    // 0011 0011 0011 0011 0011 0011 0011 0011 << 2
    // 0000 1111 0000 1111 0000 1111 0000 1111 << 4
    // 0000 0000 1111 1111 0000 0000 1111 1111 << 8
    
    //////////////////////////////////////////
    // Interleave 3x9 Table & BitHack
    
    // 0000 0001 0010 0100 1001 0010 0100 1001 << 2
    // 0000 0000 0100 1100 0011 0000 1100 0011 << 4
    // 0000 0000 0000 0001 1111 0000 0000 1111 << 8
   
    // X = (X | (X << 8)) & B[2];
    // X = (X | (X << 4)) & B[1];
    // X = (X | (X << 2)) & B[0];
    // 
    // Y = (Y | (Y << 8)) & B[2];
    // Y = (Y | (Y << 4)) & B[1];
    // Y = (Y | (Y << 2)) & B[0];
    // 
    // Z = (Z | (Z << 8)) & B[2];
    // Z = (Z | (Z << 4)) & B[1];
    // Z = (Z | (Z << 2)) & B[0];
    //
    // R = X | (Y << 1) | (Z << 2)
    
    //2 integer interleave using bithacks
    void morton_xy(int &x, int &y, int &r)
    {
        static const unsigned int B[] = {0x55555555, 0x33333333, 0x0F0F0F0F, 0x00FF00FF};
        static const unsigned int S[] = {1, 2, 4, 8};
        
        x &= 0x000FFFF;
        y &= 0x000FFFF;
        
        x = (x | (x << S[3])) & B[3];
        x = (x | (x << S[2])) & B[2];
        x = (x | (x << S[1])) & B[1];
        x = (x | (x << S[0])) & B[0];

        y = (y | (y << S[3])) & B[3];
        y = (y | (y << S[2])) & B[2];
        y = (y | (y << S[1])) & B[1];
        y = (y | (y << S[0])) & B[0];

        r = x | (y << 1);
    }
    
    //3 integer interleave using bithacks
    void morton_xyz(int &x, int &y, int &z, int &r)
    {        
        r = (((x & 0x1 ) | ((y & 0x1 ) << 1 ) | ((z & 0x1 ) << 2 ))            ) | 
            (((x & 0x2 ) | ((y & 0x2 ) << 1 ) | ((z & 0x2 ) << 2 )) << (3  - 1)) |
            (((x & 0x4 ) | ((y & 0x4 ) << 1 ) | ((z & 0x4 ) << 2 )) << (6  - 2)) |
            (((x & 0x8 ) | ((y & 0x8 ) << 1 ) | ((z & 0x8 ) << 2 )) << (9  - 3)) |
            (((x & 0x10) | ((y & 0x10) << 1 ) | ((z & 0x10) << 2 )) << (12 - 4));
    }

    void AssignPartition(SizeType NumberOfNodes, SizeType NumberOfElements, idxtype* NPart, idxtype* EPart)
    {  
        int result;
        int index;
        
        int x;
        int y;
        int z;
        
        double xD;
        double yD;
        double zD;
        
        int maxX = std::numeric_limits<int>::min();
        int maxY = std::numeric_limits<int>::min();
        int maxZ = std::numeric_limits<int>::min();
        
        double maxXD = std::numeric_limits<double>::min();
        double maxYD = std::numeric_limits<double>::min();
        double maxZD = std::numeric_limits<double>::min();
        double minXD = std::numeric_limits<double>::max();
        double minYD = std::numeric_limits<double>::max();
        double minZD = std::numeric_limits<double>::max();
        
        double CellSizeX = 1.08169f;
        double CellSizeY = 1.1914f;
        double CellSizeZ = 1.12675f;
        
        for(SizeType i_element = 0 ; i_element < NumberOfElements ; i_element++)
        {
            for(IO::ConnectivitiesContainerType::value_type::iterator i_node = mrElementsConnectivities[i_element].begin() ;
                    i_node != mrElementsConnectivities[i_element].end() ; i_node++)
            {
                IO::NodesContainerType::iterator i_nod = mrNodesContainer.begin() + (*i_node-1);
              
                xD = i_nod->X();
                yD = i_nod->Y();
                zD = i_nod->Y();
                
                maxXD = xD > maxXD ? xD : maxXD;
                maxYD = yD > maxYD ? yD : maxYD;
                maxZD = zD > maxZD ? yD : maxZD;
                
                minXD = xD < minXD ? xD : minXD;
                minYD = yD < minYD ? yD : minYD;
                minZD = zD < minZD ? zD : minZD;
            }
        }
   
        maxX = maxXD/CellSizeX;
        maxY = maxYD/CellSizeY;
        maxZ = maxZD/CellSizeZ;

//         morton_xy(maxX,maxY,result);
        morton_xyz(maxX,maxY,maxZ,result);
        
        std::cout << "Result: " <<  result;
        
        result/=(mNumberOfPartitions-1);
        
        std::cout << " " << result << std::endl;
        
        std::list<int> mid;

        for(SizeType i_element = 0 ; i_element < NumberOfElements ; i_element++)
        {
            for(IO::ConnectivitiesContainerType::value_type::iterator i_node = mrElementsConnectivities[i_element].begin() ;
                    i_node != mrElementsConnectivities[i_element].end() ; i_node++)
            { 
                IO::NodesContainerType::iterator i_nod = mrNodesContainer.begin() + (*i_node-1);

                x = i_nod->X()/CellSizeX;
                y = i_nod->Y()/CellSizeY;
                z = i_nod->Y()/CellSizeZ;
//                 morton_xy(x,y,index);
                morton_xyz(x,y,z,index);
                
                mid.push_back(index);
                
                NPart[(*i_node-1)] = index;
                EPart[i_element] = index;
            }
        }
        
        mid.sort();
        mid.unique();
        
        std::cout << mid.size() << std::endl;
        
        std::vector<int> mid_v(mid.begin(), mid.end());
        
        int mid_element[mNumberOfPartitions];
        
        for(unsigned int i = 0; i < mNumberOfPartitions; i++)
            mid_element[i] = mid_v[i*(mid.size()/mNumberOfPartitions-1)];

        for(SizeType i_element = 0 ; i_element < NumberOfElements ; i_element++)
        {
            for(IO::ConnectivitiesContainerType::value_type::iterator i_node = mrElementsConnectivities[i_element].begin() ;
                    i_node != mrElementsConnectivities[i_element].end() ; i_node++)
            { 
                IO::NodesContainerType::iterator i_nod = mrNodesContainer.begin() + (*i_node-1);
                
                int scaledIndexNode = NPart[(*i_node-1)];
                int scaledIndexElement = EPart[i_element];
                
                for(unsigned int i = 0; i < mNumberOfPartitions; i++)
                {
                    if(scaledIndexNode >= mid_element[i])
                        NPart[(*i_node-1)] = i;
                    if(scaledIndexElement >= mid_element[i])
                        EPart[i_element] = i;
                }
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


