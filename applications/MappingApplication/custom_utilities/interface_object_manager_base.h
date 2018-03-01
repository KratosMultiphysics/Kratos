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

#if !defined(KRATOS_INTERFACE_OBJECT_MANAGER_BASE_H_INCLUDED )
#define  KRATOS_INTERFACE_OBJECT_MANAGER_BASE_H_INCLUDED

// System includes
#include <unordered_map> // for CandidateManager

// External includes

// Project includes
#include "includes/define.h"

#include "interface_node.h"
#include "interface_geometry_object.h"
#include "custom_configures/interface_object_configure.h"


namespace Kratos
{
///@addtogroup ApplicationNameApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

typedef matrix<int> GraphType; // GraphColoringProcess

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

// This struct stores temporary variables, such that they don't have to be
// stored as members in the InterfaceObjectManager
struct CandidateManager
{
    std::unordered_map<int, std::vector<InterfaceObject::Pointer> > mCandidateSendObjects;
    std::unordered_map<int, std::vector<InterfaceObject::Pointer> > mCandidateReceiveObjects;
    std::unordered_map<int, std::vector<std::vector<double> > > mCandidateShapeFunctionValues;

    std::unordered_map<int, std::vector<int> > mMatchingInformation;
};

/// BaseClass for managing the InterfaceObjects
/** This class is the interface between the Searching objects and the InterfaceObjects. It is responsible
* for filling buffers, reconstructing things from buffers and the construction of the InterfaceObjects
* Look into the class description of the MapperCommunicator to see how this Object is used in the application
*/
class InterfaceObjectManagerBase
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of InterfaceObjectManagerBase
    KRATOS_CLASS_POINTER_DEFINITION(InterfaceObjectManagerBase);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Destructor.
    virtual ~InterfaceObjectManagerBase() { }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void Clear()
    {
        mSendObjects.clear();
        mReceiveObjects.clear();
        mShapeFunctionValues.clear();

        for (auto& interface_obj : mInterfaceObjects)
        {
            interface_obj->Reset(); // Set the object back to its initial state
            // i.e. reset its information whether it has been sent somewhere or
            // whether it already found a neighbor
        }
    }

    template <typename T>
    void MapInsertElement(std::unordered_map<int, T>& rMap, int Key, T& rValue)
    {
        KRATOS_DEBUG_ERROR_IF(rMap.count(Key) > 0) << "Key already present in Map!" << std::endl;

        rMap.emplace(Key, rValue);
    }

    // **********************************************************************
    // Side we want to find neighbors for aka destination *******************
    // **********************************************************************
    // ** InterfaceObjectManagerSerial and InterfaceObjectManagerParallel **

    bool AllNeighborsFound()
    {
        int all_neighbors_found = 1; // set to "1" aka "true" by default in case
        // this partition doesn't have a part of the interface!

        for (auto& interface_obj : mInterfaceObjects)
        {
            if (!interface_obj->NeighborOrApproximationFound())
            {
                all_neighbors_found = 0;
            }
        }

        // This is necessary bcs not all partitions would start a new search iteration!
        mrModelPart.GetCommunicator().MinAll(all_neighbors_found);

        return all_neighbors_found;
    }

    void CheckResults()
    {
        for (auto& interface_obj : mInterfaceObjects)
        {
            const int pairing_status = interface_obj->GetPairingStatus();
            if (pairing_status == InterfaceObject::PairingStatus::NoNeighbor)
            {
                std::cout << "MAPPER WARNING, Rank " << mCommRank
                          << "\tPoint [ "
                          << interface_obj->X() << " | "
                          << interface_obj->Y() << " | "
                          << interface_obj->Z() << " ] "
                          << "has not found a neighbor!" << std::endl;
            }
            else if (pairing_status == InterfaceObject::PairingStatus::Approximation)
            {
                std::cout << "MAPPER WARNING, Rank " << mCommRank
                          << "\tPoint [ "
                          << interface_obj->X() << " | "
                          << interface_obj->Y() << " | "
                          << interface_obj->Z() << " ] "
                          << "uses an approximation" << std::endl;
            }
        }
    }

    InterfaceObjectConfigure::ContainerType& GetInterfaceObjects()
    {
        return mInterfaceObjects;
    }

    // ***** InterfaceObjectManagerSerial *****
    virtual void GetInterfaceObjectsSerialSearch(InterfaceObjectConfigure::ContainerType& rCandidateSendObjects)
    {
        KRATOS_ERROR << "Base class function called!" << std::endl;
    }

    virtual void PostProcessReceivedResults(const InterfaceObjectConfigure::ContainerType& rCandidateSendObjects,
                                            const std::vector<double>& rDistances,
                                            const std::vector<int>& rPairingIndices)
    {
        KRATOS_ERROR << "Base class function called!" << std::endl;
    }

    // ***** InterfaceObjectManagerParallel *****
    virtual void ComputeCandidatePartitions(CandidateManager& rCandidateManager, int* pLocalCommList,
                                            int* pLocalMemorySizeArray,
                                            double* pGlobalBoundingBoxes,
                                            const bool LastIteration)
    {
        KRATOS_ERROR << "Base class function called!" << std::endl;
    }

    // Function for debugging
    virtual void PrintCandidatePartitions(const InterfaceObject::Pointer pInterfaceObject,
                                          std::vector<int>& rPartitionList)
    {
        KRATOS_ERROR << "Base class function called!" << std::endl;
    }

    virtual void PrepareMatching(CandidateManager& rCandidateManager, int* pLocalCommList,
                                 int* pLocalMemorySizeArray)
    {
        KRATOS_ERROR << "Base class function called!" << std::endl;
    }

    virtual void FillSendBufferWithMatchInformation(CandidateManager& rCandidateManager,
            int* pSendBuffer, int& rSendBufferSize,
            const int CommPartner)
    {
        KRATOS_ERROR << "Base class function called!" << std::endl;
    }

    virtual void FillBufferLocalSearch(CandidateManager& rCandidateManager,
                                       InterfaceObjectConfigure::ContainerType& rSendObjects,
                                       int& rNumObjects)
    {
        KRATOS_ERROR << "Base class function called!" << std::endl;
    }

    virtual void FillSendBufferRemoteSearch(CandidateManager& rCandidateManager,
                                            double* pSendBuffer, int& rSendBufferSize,
                                            const int CommPartner)
    {
        KRATOS_ERROR << "Base class function called!" << std::endl;
    }

    virtual void PostProcessReceivedResults(CandidateManager& rCandidateManager,
                                            const std::vector<double>& rDistances,
                                            const std::vector<int>& rPairingIndices,
                                            const int CommPartner)
    {
        KRATOS_ERROR << "Base class function called!" << std::endl;
    }

    virtual void PostProcessReceivedResults(CandidateManager& rCandidateManager,
                                            const double* pDistances,
                                            const int* pPairingIndices,
                                            const int CommPartner)
    {
        KRATOS_ERROR << "Base class function called!" << std::endl;
    }

    // **********************************************************************
    // Side where we search neighbors aka origin ****************************
    // **********************************************************************
    // ***** InterfaceObjectManagerSerial *****
    virtual void StoreSearchResults(const std::vector<double>& rDistances,
                                    const std::vector<InterfaceObject::Pointer> TempClosestResults,
                                    const std::vector<std::vector<double>> TempShapeFunctionValues)
    {
        KRATOS_ERROR << "Base class function called!" << std::endl;
    }

    // ***** InterfaceObjectManagerParallel *****
    virtual void ProcessReceiveBuffer(InterfaceObjectConfigure::ContainerType& rRemotePointList,
                                      const double* pCoordinateList, const int CoordinateListSize,
                                      int& rNumObjects)
    {
        KRATOS_ERROR << "Base class function called!" << std::endl;
    }

    virtual void FillSendBufferWithResults(double* pSendBuffer, const int SendBufferSize,
                                           const std::vector<double>& rMinDistances)
    {
        KRATOS_ERROR << "Base class function called!" << std::endl;
    }

    virtual void FillSendBufferWithResults(int* pSendBuffer, const int SendBufferSize,
                                           const std::vector<int>& rPairingIndices)
    {
        KRATOS_ERROR << "Base class function called!" << std::endl;
    }

    virtual void StoreTempSearchResults(CandidateManager& rCandidateManager,
                                        const std::vector<InterfaceObject::Pointer> TempClosestResults,
                                        const std::vector<std::vector<double>> TempShapeFunctionValues,
                                        const int CommPartner)
    {
        KRATOS_ERROR << "Base class function called!" << std::endl;
    }

    virtual void ProcessMatchInformation(CandidateManager& rCandidateManager,
                                         int* pBuffer, const int BufferSize,
                                         const int CommPartner)
    {
        KRATOS_ERROR << "Base class function called!" << std::endl;
    }

    // **********************************************************************
    // Functions for Mapping ************************************************
    // **********************************************************************
    // ** InterfaceObjectManagerSerial and InterfaceObjectManagerParallel **
    template <typename T>
    void FillBufferWithValues(const std::function<T(InterfaceObject::Pointer, const std::vector<double>&)>& FunctionPointer,
                              std::vector< T >& rBuffer)
    {
        int i = 0;

        std::vector<InterfaceObject::Pointer> interface_objects;
        if (mReceiveObjects.count(mCommRank) > 0)
        {
            interface_objects = mReceiveObjects.at(mCommRank);
        }
        // bind shape function values // TODO
        rBuffer.resize(interface_objects.size());

        for (auto interface_obj : interface_objects)
        {
            rBuffer[i] = FunctionPointer(interface_obj, mShapeFunctionValues.at(mCommRank)[i]);
            ++i;
        }
    }

    template <typename T>
    void ProcessValues(const std::function<void(InterfaceObject::Pointer, T)>& FunctionPointer,
                       const std::vector< T >& rBuffer)
    {
        std::vector<InterfaceObject::Pointer> interface_objects;
        if (mSendObjects.count(mCommRank) > 0)
        {
            interface_objects = mSendObjects.at(mCommRank);
        }

        KRATOS_DEBUG_ERROR_IF_NOT(interface_objects.size() == rBuffer.size())
            << "Wrong number of results received!;"
            << " \"interface_objects.size() = "
            << interface_objects.size() << ", rBuffer.size() = "
            << rBuffer.size() << std::endl;

        for (std::size_t i = 0; i < interface_objects.size(); ++i)
        {
            FunctionPointer(interface_objects[i], rBuffer[i]);
        }
    }

    // ***** InterfaceObjectManagerParallel *****
    virtual void ComputeBufferSizesAndCommunicationGraph(int& rMaxSendBufferSize,
            int& rMaxReceiveBufferSize,
            GraphType& rColoredGraph,
            int& rMaxColors)
    {
        KRATOS_ERROR << "Base class function called!" << std::endl;
    }

    virtual void FillBufferWithValues(double* pBuffer, int& rBufferSize, const int CommPartner,
                                      const std::function<double(InterfaceObject::Pointer, const std::vector<double>&)>& FunctionPointer)
    {
        KRATOS_ERROR << "Base class function called!" << std::endl;
    }

    virtual void FillBufferWithValues(double* pBuffer, int& rBufferSize, const int CommPartner,
                                      const std::function<array_1d<double, 3>(InterfaceObject::Pointer, const std::vector<double>&)>& FunctionPointer)
    {
        KRATOS_ERROR << "Base class function called!" << std::endl;
    }

    virtual void ProcessValues(const double* pBuffer, const int BufferSize, const int CommPartner,
                               const std::function<void(InterfaceObject::Pointer, double)>& FunctionPointer)
    {
        KRATOS_ERROR << "Base class function called!" << std::endl;
    }

    virtual void ProcessValues(const double* pBuffer, const int BufferSize, const int CommPartner,
                               const std::function<void(InterfaceObject::Pointer, array_1d<double, 3>)>& FunctionPointer)
    {
        KRATOS_ERROR << "Base class function called!" << std::endl;
    }

    // Functions used for Debugging
    void WriteNeighborRankAndCoordinates()
    {
        for (auto& interface_obj : mInterfaceObjects)
        {
            interface_obj->WriteRankAndCoordinatesToVariable(mCommRank);
        }
    }

    void PrintNeighbors()
    {
        for (auto& interface_obj : mInterfaceObjects)
        {
            interface_obj->PrintNeighbors(mCommRank);
        }
    }

    void PrintInterfaceObjects(const std::string& rInterfaceSide)
    {
        for (auto& r_interface_obj : mInterfaceObjects)
        {
            std::cout << rInterfaceSide << " , Rank " << mCommRank << " , "
                      << r_interface_obj->Info() << " [ "
                      << r_interface_obj->X() << " "
                      << r_interface_obj->Y() << " "
                      << r_interface_obj->Z() << "]" << std::endl;
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
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "InterfaceObjectManagerBase" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "InterfaceObjectManagerBase";
    }

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

    InterfaceObjectManagerBase(ModelPart& rModelPart, int CommRank, int CommSize,
                               MapperUtilities::InterfaceObjectConstructionType InterfaceObjectType,
                               GeometryData::IntegrationMethod IntegrationMethod, const int EchoLevel,
                               const double ApproximationTolerance) :
        mrModelPart(rModelPart)
    {

        mCommRank = CommRank;
        mCommSize = CommSize;

        mEchoLevel = EchoLevel;

        if (InterfaceObjectType == MapperUtilities::Node_Coords)
        {
            InitializeInterfaceNodeManager(rModelPart);
        }
        else if (InterfaceObjectType == MapperUtilities::Condition_Center ||
                 InterfaceObjectType == MapperUtilities::Condition_Gauss_Point)
        {
            InitializeInterfaceGeometryObjectManager(rModelPart, IntegrationMethod, ApproximationTolerance);
        }
        else
        {
            KRATOS_ERROR << "Type of interface object construction not implemented" << std::endl;
        }

        int num_interface_objects = mInterfaceObjects.size();
        mrModelPart.GetCommunicator().SumAll(num_interface_objects);

        KRATOS_ERROR_IF_NOT(num_interface_objects > 0) 
            << "No interface objects were created in ModelPart \""
            << mrModelPart.Name() << "\"!" << std::endl;
    }

    ModelPart& mrModelPart;

    InterfaceObjectConfigure::ContainerType mInterfaceObjects;

    int mCommRank = 0;
    int mCommSize = 0;
    int mEchoLevel = 0;

    // point-sending interface (destination)
    std::unordered_map<int, std::vector<InterfaceObject::Pointer> > mSendObjects;

    // point-receiving interface
    std::unordered_map<int, std::vector<InterfaceObject::Pointer> > mReceiveObjects;
    std::unordered_map<int, std::vector<std::vector<double> > > mShapeFunctionValues;


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


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void InitializeInterfaceNodeManager(ModelPart& rModelPart)
    {
        mInterfaceObjects.resize(rModelPart.GetCommunicator().LocalMesh().NumberOfNodes());

        int i = 0;
        for (auto &local_node : rModelPart.GetCommunicator().LocalMesh().Nodes())
        {
            mInterfaceObjects[i] = InterfaceObject::Pointer( new InterfaceNode(local_node, mEchoLevel) );
            ++i;
        }
    }

    void InitializeInterfaceGeometryObjectManager(ModelPart& rModelPart,
            GeometryData::IntegrationMethod IntegrationMethod,
            const double ApproximationTolerance)
    {
        bool construct_with_center;
        int size_factor = 1;
        if (IntegrationMethod == GeometryData::NumberOfIntegrationMethods)
        {
            construct_with_center = true;
            size_factor = 1;
        }
        else
        {
            construct_with_center = false;
            if (rModelPart.GetCommunicator().LocalMesh().NumberOfConditions() > 0)
            {
                size_factor = rModelPart.GetCommunicator().LocalMesh().ConditionsBegin()->GetGeometry().IntegrationPointsNumber(IntegrationMethod);
            }
            else if (rModelPart.GetCommunicator().LocalMesh().NumberOfElements() > 0)
            {
                size_factor = rModelPart.GetCommunicator().LocalMesh().ElementsBegin()->GetGeometry().IntegrationPointsNumber(IntegrationMethod);
            }
        }

        mInterfaceObjects.reserve(size_factor * rModelPart.GetCommunicator().LocalMesh().NumberOfConditions());

        if (construct_with_center)   // construct with condition center point
        {
            for (auto& condition : rModelPart.GetCommunicator().LocalMesh().Conditions())
            {
                mInterfaceObjects.push_back(InterfaceObject::Pointer( new InterfaceGeometryObject(condition.GetGeometry(),
                                            ApproximationTolerance,
                                            mEchoLevel, 
                                            0) ));
            }
            for (auto& element : rModelPart.GetCommunicator().LocalMesh().Elements())
            {
                mInterfaceObjects.push_back(InterfaceObject::Pointer( new InterfaceGeometryObject(element.GetGeometry(),
                                            ApproximationTolerance,
                                            mEchoLevel,
                                            0) ));
            }
        }
        else     // construct with condition gauss points
        {
            KRATOS_ERROR << "This is not implemented at the moment" << std::endl;
            for (auto& condition : rModelPart.GetCommunicator().LocalMesh().Conditions())
            {
                for (int g = 0; g < 111111; ++g) // TODO fix this, should be number of GPs
                {
                    mInterfaceObjects.push_back(InterfaceObject::Pointer( new InterfaceGeometryObject(condition.GetGeometry(),
                                                ApproximationTolerance,
                                                mEchoLevel,
                                                g, IntegrationMethod) ));
                }
            }
            for (auto& element : rModelPart.GetCommunicator().LocalMesh().Elements())
            {
                for (int g = 0; g < 111111; ++g) // TODO fix this, should be number of GPs
                {
                    mInterfaceObjects.push_back(InterfaceObject::Pointer( new InterfaceGeometryObject(element.GetGeometry(),
                                                ApproximationTolerance,
                                                mEchoLevel,
                                                g, IntegrationMethod) ));
                }
            }

        }
    }

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
    InterfaceObjectManagerBase& operator=(InterfaceObjectManagerBase const& rOther);

    //   /// Copy constructor.
    //   InterfaceObjectManagerBase(InterfaceObjectManagerBase const& rOther){}


    ///@}

}; // Class InterfaceObjectManagerBase

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  InterfaceObjectManagerBase& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const InterfaceObjectManagerBase& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block


}  // namespace Kratos.

#endif // KRATOS_INTERFACE_OBJECT_MANAGER_BASE_H_INCLUDED  defined
