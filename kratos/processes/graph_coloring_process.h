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

#if !defined(KRATOS_GRAPH_COLORING_PROCESS_H_INCLUDED )
#define  KRATOS_GRAPH_COLORING_PROCESS_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/ublas_interface.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/// Getting a graph of connectivities (of domains) and performs a coloring procedure.
/** This class takes a graph of connectivites in form of a dense matrix and
    calculates a coloring schem where no neighbours with the same color are
  allowed.
  This process uses a very simple algorithm that for each domain make a loop over its
  neighbours and finds the first unused color for this pair of domain. The cost of this
  procedure is quadratic and make it in efficient for large number of domains.

  @see Execute
*/
class GraphColoringProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of GraphColoringProcess
    KRATOS_CLASS_POINTER_DEFINITION(GraphColoringProcess);
    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

    /** Defining a dense matrix of integer as graph type
    */
    typedef DenseMatrix<int> GraphType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    /** This constructor takes all necessary data to perform the coloring.

        @param NumberOfPartitions Number of partitions (or domains)

      @param rDomainsGraph The graph of domains with their neighbours presented as a matrix where any connection between
             two domains is presented by putting 1 in corresponding row and columns

      @param rDomainsColoredGraph The result of the process. Each column represents a color and the neighbour domain id for each color
             is stored in this matrix. The rest of the values are set to -1 which shows no neighbour
    		 for the color.

      @param rMaxColor After executing the process, the maximum number of colors will be stored in this variable

        @see GraphType
    */
    GraphColoringProcess(int NumberOfPartitions, const GraphType& rDomainsGraph, GraphType& rDomainsColoredGraph, int& rMaxColor):
        mNumberOfPartitions(NumberOfPartitions), mrMaxColor(rMaxColor), mrDomainsGraph(rDomainsGraph), mrDomainsColoredGraph(rDomainsColoredGraph)
    {}

    /// Destructor.
    virtual ~GraphColoringProcess() {}


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
		//// Applying the Misra and Gries edge coloring algorithm.
  //      // Start coloring...
  //      for(SizeType i = 0 ; i < mrDomainsGraph.size1() ; i++) // for each domain
  //          for(SizeType j = i + 1 ; j < mrDomainsGraph.size2() ; j++) // finding neighbor domains
  //              if(mrDomainsGraph(i,j) != 0.00) // domain i has interface with domain j
  // I should continue the implementation. Pooyan.

        mrMaxColor = 0;
        // Initializing the coloered graph. -1 means no connection
        // TODO: I have to change this part and create this matrix using max color number
        mrDomainsColoredGraph.resize(mNumberOfPartitions, 2*mNumberOfPartitions,false);
        mrDomainsColoredGraph = ScalarMatrix(mNumberOfPartitions, 2*mNumberOfPartitions, -1.00);

        // Start coloring...
        for(SizeType i = 0 ; i < mrDomainsGraph.size1() ; i++) // for each domain
            for(SizeType j = i + 1 ; j < mrDomainsGraph.size2() ; j++) // finding neighbor domains
                if(mrDomainsGraph(i,j) != 0.00) // domain i has interface with domain j
                    for(SizeType color = 0 ; color < mrDomainsColoredGraph.size2() ; color++) // finding color
                        if((mrDomainsColoredGraph(i,color) == -1.00) && (mrDomainsColoredGraph(j,color) == -1.00)) // the first unused color
                        {
                            mrDomainsColoredGraph(i,color) = j;
                            mrDomainsColoredGraph(j,color) = i;
                            if(mrMaxColor < static_cast<int>(color + 1))
                                mrMaxColor = color + 1;
                            break;
                        }
    }



    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "GraphColoringProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "GraphColoringProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}


private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    int mNumberOfPartitions;
    int& mrMaxColor;
    const GraphType& mrDomainsGraph;
    GraphType& mrDomainsColoredGraph;


    /// Assignment operator.
    GraphColoringProcess& operator=(GraphColoringProcess const& rOther);

    /// Copy constructor.
    //GraphColoringProcess(GraphColoringProcess const& rOther);


    ///@}

}; // Class GraphColoringProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  GraphColoringProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const GraphColoringProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_GRAPH_COLORING_PROCESS_H_INCLUDED  defined
