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

#if !defined(KRATOS_MAPPER_LOCAL_SYSTEM_H_INCLUDED )
#define  KRATOS_MAPPER_LOCAL_SYSTEM_H_INCLUDED


// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "mapper_interface_info.h"


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

/// This is the "Condition" of the Mappers.
/** Detail class definition.
In case the mapping matrix is build we need at some point the EquationIDs of the MappingMatrix
=> To get these we have to compute also the MappingWeights which might be expensive
(e.g. triangulation, integration, projection)
=> Therefore we save this to later only querry it
*/
class MapperLocalSystem
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MapperLocalSystem
    KRATOS_CLASS_POINTER_DEFINITION(MapperLocalSystem);

    typedef Kratos::shared_ptr<MapperInterfaceInfo> MapperInterfaceInfoPointerType;
    typedef Kratos::unique_ptr<MapperLocalSystem> MapperLocalSystemUniquePointer;

    typedef typename MapperInterfaceInfo::CoordinatesArrayType CoordinatesArrayType;

    typedef Matrix MatrixType;
    typedef std::vector<int> EquationIdVectorType; // int bcs of mpi

    typedef InterfaceObject::NodePointerType NodePointerType;
    typedef InterfaceObject::GeometryPointerType GeometryPointerType;

    ///@}
    ///@name  Enum's
    ///@{

    enum PairingStatus
    {
        NoInterfaceInfo,
        Approximation,
        InterfaceInfoFound
    };


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MapperLocalSystem() {}

    /// Destructor.
    virtual ~MapperLocalSystem() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    // Only one of the Create functions has to be implemented, thats why they cannot be pure virtual!
    virtual MapperLocalSystemUniquePointer Create(NodePointerType pNode) const
    {
        KRATOS_ERROR << "Base class function called!" << std::endl;
    }

    virtual MapperLocalSystemUniquePointer Create(GeometryPointerType pGeometry) const
    {
        KRATOS_ERROR << "Base class function called!" << std::endl;
    }

    void EquationIdVectors(EquationIdVectorType& rOriginIds,
                           EquationIdVectorType& rDestinationIds)
    {
        if (!mIsComputed) {
            CalculateAll(mLocalMappingMatrix, mOriginIds, mDestinationIds, mPairingStatus);
            mIsComputed = true;
        }

        rOriginIds      = mOriginIds;
        rDestinationIds = mDestinationIds;
    }

    void CalculateLocalSystem(MatrixType& rLocalMappingMatrix,
                              EquationIdVectorType& rOriginIds,
                              EquationIdVectorType& rDestinationIds) // TODO should be const if it werent for the PairingStatus
    {
        if (mIsComputed) {
            // This will be called if the EquationIdVectors have been querried before
            // i.e. matrix-based mapping
            rLocalMappingMatrix = mLocalMappingMatrix;
            rOriginIds      = mOriginIds;
            rDestinationIds = mDestinationIds;
        }
        else {
            // This will be called if the EquationIdVectors have NOT been querried before
            // i.e. matrix-free mapping
            CalculateAll(rLocalMappingMatrix, rOriginIds, rDestinationIds, mPairingStatus);
        }
    }

    /**
    * @brief Resizing the output if no InterfaceInfo is available
    * This function resizes the system vectors to zero and also sets that no valid
    * Information from the other side could be found to compute the local system
    * @param rLocalMappingMatrix The vector conatining the mapping weights
    * @param rOriginIds The vector containing the ids on the origin
    * @param rDestinationIds The vector containing the ids on the destination
    * @param rPairingStatus The pairingstatus of the MapperLocalSystem
    * @see CalculateAll
    * @author Philipp Bucher
    */
    void ResizeToZero(MatrixType& rLocalMappingMatrix,
                      EquationIdVectorType& rOriginIds,
                      EquationIdVectorType& rDestinationIds,
                      MapperLocalSystem::PairingStatus& rPairingStatus) const
    {
        rPairingStatus = MapperLocalSystem::PairingStatus::NoInterfaceInfo;

        rLocalMappingMatrix.resize(0, 0, false);
        rOriginIds.resize(0);
        rDestinationIds.resize(0);
    }


    // // This specifies if Nodes should be used for the construction
    // // => this is the case if the Geometry on the destination is not important
    virtual bool UseNodesAsBasis() const = 0;

    virtual CoordinatesArrayType& Coordinates() const = 0;


    void AddInterfaceInfo(MapperInterfaceInfoPointerType pInterfaceInfo) // TODO pass by const ref?
    {
        mInterfaceInfos.push_back(pInterfaceInfo);
    }

    bool HasInterfaceInfo() const
    {
        return mInterfaceInfos.size() > 0;
    }


    virtual void Clear()
    {
        mInterfaceInfos.clear();
        mLocalMappingMatrix.clear();
        mOriginIds.clear();
        mDestinationIds.clear();
        mIsComputed = false;
    }

    PairingStatus GetPairingStatus() const
    {
        return mPairingStatus;
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

    virtual std::string PairingInfo(const int EchoLevel, const int CommRank) const = 0;

    /// Turn back information as a string.
    virtual std::string Info() const {return "MapperLocalSystem";}

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}


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

    std::vector<MapperInterfaceInfoPointerType> mInterfaceInfos;

    bool mIsComputed = false;

    MatrixType mLocalMappingMatrix;
    EquationIdVectorType mOriginIds;
    EquationIdVectorType mDestinationIds;

    PairingStatus mPairingStatus = PairingStatus::NoInterfaceInfo;

    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    // This function calculates the components necessary for the mapping
    // Note that it is "const", therefore it can NOT modify its members
    // Whether members are to be saved is decided in other functions of this class
    virtual void CalculateAll(MatrixType& rLocalMappingMatrix,
                              EquationIdVectorType& rOriginIds,
                              EquationIdVectorType& rDestinationIds,
                              MapperLocalSystem::PairingStatus& rPairingStatus) const = 0;


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

    // /// Assignment operator.
    // MapperLocalSystem& operator=(MapperLocalSystem const& rOther);

    // /// Copy constructor.
    // MapperLocalSystem(MapperLocalSystem const& rOther);


    ///@}

}; // Class MapperLocalSystem

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


// /// input stream function
// inline std::istream& operator >> (std::istream& rIStream,
//                 MapperLocalSystem& rThis);

// /// output stream function
// inline std::ostream& operator << (std::ostream& rOStream,
//                 const MapperLocalSystem& rThis)
// {
//     rThis.PrintInfo(rOStream);
//     rOStream << std::endl;
//     rThis.PrintData(rOStream);

//     return rOStream;
// }
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MAPPER_LOCAL_SYSTEM_H_INCLUDED  defined
