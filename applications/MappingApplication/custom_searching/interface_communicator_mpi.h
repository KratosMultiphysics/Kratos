//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

#if !defined(KRATOS_INTERFACE_COMMUNICATOR_MPI_H_INCLUDED )
#define  KRATOS_INTERFACE_COMMUNICATOR_MPI_H_INCLUDED

// System includes

// External includes

// Project includes
#include "interface_communicator.h"


namespace Kratos
{
///@addtogroup MappingApplication
///@{

///@name Kratos Classes
///@{

/// Object for exchanging data on the Interface in MPI
/** It extends it's baseclass by remote-searching capabilities. I.e. before the local search,
 * data is exchanged among the partitions to be used by the local-search
*/
class InterfaceCommunicatorMPI : public InterfaceCommunicator
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of InterfaceCommunicatorMPI
    KRATOS_CLASS_POINTER_DEFINITION(InterfaceCommunicatorMPI);

    typedef std::vector<std::vector<double>> BufferTypeDouble;
    typedef std::vector<std::vector<char>> BufferTypeChar;

    ///@}
    ///@name Life Cycle
    ///@{

    InterfaceCommunicatorMPI(ModelPart& rModelPartOrigin,
                                MapperLocalSystemPointerVectorPointer pMapperLocalSystems,
                                Parameters SearchSettings);

    /// Destructor.
    virtual ~InterfaceCommunicatorMPI()
    {
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "InterfaceCommunicatorMPI" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "InterfaceCommunicatorMPI";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {}

    ///@}

protected:
    ///@name Protected Operations
    ///@{

    void InitializeSearch(const Kratos::Flags& rOptions,
                       const MapperInterfaceInfoUniquePointerType& rpRefInterfaceInfo,
                       InterfaceObject::ConstructionType InterfaceObjectTypeOrigin) override;

    // This function constructs the InterfaceObjects on the Destination
    // In serial it only does it once, whereas in MPI this involves Data-Exchange!
    // Imagine a sliding interface, there the partitions might change!
    void InitializeSearchIteration(const Kratos::Flags& rOptions,
                          const MapperInterfaceInfoUniquePointerType& rpRefInterfaceInfo) override;

    void FinalizeSearchIteration(const MapperInterfaceInfoUniquePointerType& rpRefInterfaceInfo) override;

    ///@}

private:
    ///@name Member Variables
    ///@{

    std::vector<double> mGlobalBoundingBoxes;
    // xmax, xmin,  ymax, ymin,  zmax, zmin

    int mCommRank;
    int mCommSize;

    std::vector<int> mSendSizes;
    std::vector<int> mRecvSizes;

    BufferTypeDouble mSendBufferDouble;
    BufferTypeDouble mRecvBufferDouble;

    BufferTypeChar mSendBufferChar;
    BufferTypeChar mRecvBufferChar;

    ///@}
    ///@name Private Operations
    ///@{

    std::size_t GetBufferSizeEstimate() const
    {
        return mpMapperLocalSystems->size() / mCommSize;
    }

    void ComputeGlobalBoundingBoxes();

    template< typename TDataType >
    int ExchangeDataAsync(
        const std::vector<std::vector<TDataType>>& rSendBuffer,
        std::vector<std::vector<TDataType>>& rRecvBuffer);

    ///@}

}; // Class InterfaceCommunicatorMPI

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_INTERFACE_COMMUNICATOR_MPI_H_INCLUDED  defined
