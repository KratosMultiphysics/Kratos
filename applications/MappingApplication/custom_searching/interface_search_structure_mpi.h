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

#if !defined(KRATOS_INTERFACE_SEARCH_STRUCTURE_MPI_H_INCLUDED )
#define  KRATOS_INTERFACE_SEARCH_STRUCTURE_MPI_H_INCLUDED

// System includes

// External includes

// Project includes
#include "interface_search_structure_base.h"


namespace Kratos
{
///@addtogroup MappingApplication
///@{

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

/// MPI-parallel searching
/** This class provides features for remote searching. It first computes candidate partitions
* (partitions in the vicinity of a point), where possible neighbors can be.
* It then send the objects for which neighbors are to be found to the corresponding candidate partitions
* In the candidate partitions (can be either local or remote), a local search is conducted (see BaseClass)
* The results are sent back to the partition where the object is local, and the best result is then chosen.
* Look into the class description of the MapperCommunicator to see how this Object is used in the application
*/
class InterfaceSearchStructureMPI : public InterfaceSearchStructureBase
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of InterfaceSearchStructureMPI
    KRATOS_CLASS_POINTER_DEFINITION(InterfaceSearchStructureMPI);

    using BaseType = InterfaceSearchStructureBase;

    using MapperLocalSystemPointer = typename BaseType::MapperLocalSystemPointer;
    using MapperLocalSystemPointerVector = typename BaseType::MapperLocalSystemPointerVector;
    using MapperLocalSystemPointerVectorPointer = typename BaseType::MapperLocalSystemPointerVectorPointer;

    using BufferTypeDouble = std::vector<std::vector<double>>;
    using BufferTypeChar   = std::vector<std::vector<char>>;

    ///@}
    ///@name Life Cycle
    ///@{

    InterfaceSearchStructureMPI(ModelPart& rModelPartOrigin,
                                MapperLocalSystemPointerVectorPointer pMapperLocalSystems,
                                Parameters SearchSettings);

    /// Destructor.
    virtual ~InterfaceSearchStructureMPI()
    {
    }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


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
        std::stringstream buffer;
        buffer << "InterfaceSearchStructureMPI" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "InterfaceSearchStructureMPI";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {}


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

    void PrepareSearch(const Kratos::Flags& rOptions,
                       const MapperInterfaceInfoUniquePointerType& rpRefInterfaceInfo,
                       InterfaceObject::ConstructionType InterfaceObjectTypeOrigin) override;

    // This function constructs the InterfaceObjects on the Destination
    // In serial it only does it once, whereas in MPI this involves Data-Exchange!
    // Imagine a sliding interface, there the partitions might change!
    void InitializeSearchIteration(const Kratos::Flags& rOptions,
                          const MapperInterfaceInfoUniquePointerType& rpRefInterfaceInfo) override;

    void FinalizeSearchIteration(const MapperInterfaceInfoUniquePointerType& rpRefInterfaceInfo) override;

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
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    std::size_t GetBufferSizeEstimate() const
    {
        return mpMapperLocalSystems->size() / mCommSize;
    }

    void ComputeGlobalBoundingBoxes();

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
    // InterfaceSearchStructureMPI& operator=(InterfaceSearchStructureMPI const& rOther) {}

    //   /// Copy constructor.
    //   InterfaceSearchStructureMPI(InterfaceSearchStructureMPI const& rOther){}


    ///@}

}; // Class InterfaceSearchStructureMPI

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_INTERFACE_SEARCH_STRUCTURE_MPI_H_INCLUDED  defined
