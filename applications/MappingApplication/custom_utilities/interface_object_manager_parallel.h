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

#if !defined(KRATOS_INTERFACE_OBJECT_MANAGER_PARALLEL_H_INCLUDED )
#define  KRATOS_INTERFACE_OBJECT_MANAGER_PARALLEL_H_INCLUDED

// System includes

// External includes

// Project includes
#include "interface_object_manager_base.h"
#include "mapper_utilities_mpi.h"


namespace Kratos
{
///@addtogroup ApplicationNameApplication
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

/// MPI-Parallel Verison of the Entity that manages the InterfaceObjects
/** It implements the functions that are only needed if the mapper is a mpi-parallel mapper. 
* These functions are implemented as virtual functions in the BaseClass. Besides handeling 
* buffers it also computes the communication graph and the buffer sizes.
* Look into the class description of the MapperCommunicator to see how this Object is used in the application
*/
class InterfaceObjectManagerParallel : public InterfaceObjectManagerBase
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of InterfaceObjectManagerParallel
    KRATOS_CLASS_POINTER_DEFINITION(InterfaceObjectManagerParallel);

    ///@}
    ///@name Life Cycle
    ///@{

    InterfaceObjectManagerParallel(ModelPart& rModelPart, int CommRank, int CommSize,
                                   MapperUtilities::InterfaceObjectConstructionType InterfaceObjectType,
                                   GeometryData::IntegrationMethod IntegrationMethod, const int EchoLevel,
                                   const double ApproximationTolerance) :
        InterfaceObjectManagerBase(rModelPart, CommRank, CommSize,
                                   InterfaceObjectType, IntegrationMethod, EchoLevel, ApproximationTolerance) {}

    /// Destructor.
    virtual ~InterfaceObjectManagerParallel() { }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    // **********************************************************************
    // Side we want to find neighbors for aka destination *******************
    // **********************************************************************
    void ComputeCandidatePartitions(CandidateManager& rCandidateManager, int* pLocalCommList,
                                    int* pLocalMemorySizeArray,
                                    double* pGlobalBoundingBoxes,
                                    const bool LastIteration) override
    {
        double* bounding_box[6];
        std::vector<int> partition_list; // For debugging
        for (auto& interface_obj : mInterfaceObjects)
        {
            for (int partition_index = 0; partition_index < mCommSize; ++partition_index)   // loop over partitions
            {
                for (int j = 0; j < 6; ++j)   // retrieve bounding box of partition
                {
                    bounding_box[j] = &pGlobalBoundingBoxes[(partition_index * 6) + j];
                }
                if (!interface_obj->NeighborOrApproximationFound())   // check if the interface object already found a neighbor
                {
                    if (interface_obj->IsInBoundingBox(bounding_box))
                    {
                        rCandidateManager.mCandidateSendObjects[partition_index].push_back(interface_obj);
                        pLocalCommList[partition_index] = 1;
                        ++pLocalMemorySizeArray[partition_index];
                        interface_obj->SetIsBeingSent();
                        partition_list.push_back(partition_index); // For debugging
                    }
                }
            }

            if (mEchoLevel >= 4)
            {
                PrintCandidatePartitions(interface_obj, partition_list); // For debugging
            }

            if (LastIteration)
            {
                // Robustness check, if interface_obj has not found a neighbor, it is sent to every partition
                if (!interface_obj->GetIsBeingSent())
                {
                    // Send interface_obj to all Partitions
                    if (mEchoLevel >= 2)
                    {
                        std::cout << "MAPPER WARNING, Rank " << mCommRank
                                  << ", interface_obj [ "
                                  << interface_obj->X() << " "
                                  << interface_obj->Y() << " "
                                  << interface_obj->Z() << " ] has not found "
                                  << "a neighbor/approximation yet and "
                                  << "is sent to all partitions!" << std::endl;
                    }

                    for (int partition_index = 0; partition_index < mCommSize; ++partition_index)
                    {
                        rCandidateManager.mCandidateSendObjects[partition_index].push_back(interface_obj);
                        pLocalCommList[partition_index] = 1;
                        ++pLocalMemorySizeArray[partition_index];
                    }
                }
            }
        }
    }

    // Function for debugging
    void PrintCandidatePartitions(const InterfaceObject::Pointer pInterfaceObject,
                                  std::vector<int>& rPartitionList) override
    {
        std::cout << "Rank " << mCommRank << ", interface_obj ["
                  << pInterfaceObject->X() << " "
                  << pInterfaceObject->Y() << " "
                  << pInterfaceObject->Z() << "] "
                  << " sent to ranks: [";
        for (size_t i = 0; i < rPartitionList.size(); ++i)
        {
            std::cout << rPartitionList[i] << " ";
        }
        std::cout << "]" << std::endl;
        rPartitionList.clear();
    }

    void PrepareMatching(CandidateManager& rCandidateManager, int* pLocalCommList,
                         int* pLocalMemorySizeArray) override
    {
        for (int partition_index = 0; partition_index < mCommSize; ++partition_index)
        {
            if (rCandidateManager.mCandidateSendObjects.count(partition_index) > 0)
            {
                rCandidateManager.mMatchingInformation.reserve(rCandidateManager.mCandidateSendObjects.at(partition_index).size());
                for (auto interface_obj : rCandidateManager.mCandidateSendObjects.at(partition_index))
                {
                    if (interface_obj->HasNeighborOrApproximationInPartition(partition_index))
                    {
                        mSendObjects[partition_index].push_back(interface_obj);
                        rCandidateManager.mMatchingInformation[partition_index].push_back(1);
                        pLocalCommList[partition_index] = 1;
                        ++pLocalMemorySizeArray[partition_index];
                    }
                    else
                    {
                        rCandidateManager.mMatchingInformation[partition_index].push_back(0);
                    }
                }
            }
        }
    }

    void FillSendBufferWithMatchInformation(CandidateManager& rCandidateManager,
                                            int* pSendBuffer, int& rSendBufferSize,
                                            const int CommPartner) override
    {
        int i = 0;
        if (rCandidateManager.mMatchingInformation.count(CommPartner) > 0)
        {
            for (auto info : rCandidateManager.mMatchingInformation.at(CommPartner))
            {
                pSendBuffer[i] = info;
                ++i;
            }
        }
        rSendBufferSize = i;
    }

    void FillBufferLocalSearch(CandidateManager& rCandidateManager,
                               InterfaceObjectConfigure::ContainerType& rSendObjects,
                               int& rNumObjects) override
    {
        int i = 0;
        if (rCandidateManager.mCandidateSendObjects.count(mCommRank) > 0)
        {
            for (auto interface_obj : rCandidateManager.mCandidateSendObjects.at(mCommRank))
            {
                rSendObjects[i] = interface_obj;
                ++i;
            }
        }
        rNumObjects = i;
    }

    void FillSendBufferRemoteSearch(CandidateManager& rCandidateManager,
                                    double* pSendBuffer, int& rSendBufferSize,
                                    const int CommPartner) override
    {
        int i = 0;
        if (rCandidateManager.mCandidateSendObjects.count(CommPartner) > 0)
        {
            for (auto interface_obj : rCandidateManager.mCandidateSendObjects.at(CommPartner))
            {
                pSendBuffer[(i * 3) + 0] = interface_obj->X();
                pSendBuffer[(i * 3) + 1] = interface_obj->Y();
                pSendBuffer[(i * 3) + 2] = interface_obj->Z();
                ++i;
            }
        }
        rSendBufferSize = 3 * i;
    }

    void PostProcessReceivedResults(CandidateManager& rCandidateManager,
                                    const std::vector<double>& rDistances,
                                    const std::vector<int>& rPairingIndices,
                                    const int CommPartner) override
    {
        int i = 0;
        if (rCandidateManager.mCandidateSendObjects.count(CommPartner) > 0)
        {
            for (auto interface_obj : rCandidateManager.mCandidateSendObjects.at(CommPartner))
            {
                if (rDistances[i] > -0.5f) // failed search has value "-1"
                    interface_obj->ProcessSearchResult(rDistances[i], rPairingIndices[i], CommPartner);
                ++i;
            }
        }
    }

    void PostProcessReceivedResults(CandidateManager& rCandidateManager,
                                    const double* pDistances,
                                    const int* pPairingIndices,
                                    const int CommPartner) override
    {
        int i = 0;
        if (rCandidateManager.mCandidateSendObjects.count(CommPartner) > 0)
        {
            for (auto interface_obj : rCandidateManager.mCandidateSendObjects.at(CommPartner))
            {
                if (pDistances[i] > -0.5f) // failed search has value "-1"
                    interface_obj->ProcessSearchResult(pDistances[i], pPairingIndices[i], CommPartner);
                ++i;
            }
        }
    }

    // **********************************************************************
    // Side where we search neighbors aka origin ****************************
    // **********************************************************************
    void ProcessReceiveBuffer(InterfaceObjectConfigure::ContainerType& rRemotePointList,
                              const double* pCoordinateList, const int CoordinateListSize,
                              int& rNumObjects) override
    {
        rNumObjects = CoordinateListSize / 3;

        for (int i = 0; i < rNumObjects; ++i)   // create InterfaceObjects
        {
            rRemotePointList[i] = InterfaceObject::Pointer(new InterfaceObject(
                                      pCoordinateList[(i * 3) + 0], pCoordinateList[(i * 3) + 1], pCoordinateList[(i * 3) + 2]));
        }
    }

    void FillSendBufferWithResults(double* pSendBuffer, const int SendBufferSize,
                                   const std::vector<double>& rMinDistances) override
    {
        for (int i = 0; i < SendBufferSize; ++i)
            pSendBuffer[i] = rMinDistances[i];
    }

    void FillSendBufferWithResults(int* pSendBuffer, const int SendBufferSize,
                                   const std::vector<int>& rPairingIndices) override
    {
        for (int i = 0; i < SendBufferSize; ++i)
            pSendBuffer[i] = rPairingIndices[i];
    }

    void StoreTempSearchResults(CandidateManager& rCandidateManager,
                                std::vector<InterfaceObject::Pointer> TempClosestResults,
                                std::vector<std::vector<double>> TempShapeFunctionValues,
                                const int CommPartner) override
    {
        MapInsertElement(rCandidateManager.mCandidateReceiveObjects, CommPartner, TempClosestResults);
        MapInsertElement(rCandidateManager.mCandidateShapeFunctionValues, CommPartner, TempShapeFunctionValues);
    }

    void ProcessMatchInformation(CandidateManager& rCandidateManager,
                                 int* pBuffer, const int BufferSize,
                                 const int CommPartner) override
    {
        for (int i = 0; i < BufferSize; ++i)
        {
            if (pBuffer[i] == 1)   // Match
            {
                KRATOS_DEBUG_ERROR_IF_NOT(rCandidateManager.mCandidateReceiveObjects.at(CommPartner)[i]) 
                    << "interface_obj pointer mismatch"
                    << std::endl;

                mReceiveObjects[CommPartner].push_back(rCandidateManager.mCandidateReceiveObjects.at(CommPartner)[i]);
                mShapeFunctionValues[CommPartner].push_back(rCandidateManager.mCandidateShapeFunctionValues.at(CommPartner)[i]);
            }
        }
    }

    // **********************************************************************
    // Functions for Mapping ************************************************
    // **********************************************************************
    void ComputeBufferSizesAndCommunicationGraph(int& rMaxSendBufferSize,
            int& rMaxReceiveBufferSize,
            GraphType& rColoredGraph,
            int& rMaxColors) override
    {
        int* local_comm_list = new int[mCommSize]();
        int* local_memory_size_array = new int[mCommSize]();

        for (auto& interface_obj : mInterfaceObjects)
        {
            if (interface_obj->NeighborOrApproximationFound())   // check if the interface object already found a neighbor
            {
                int neighbor_rank = interface_obj->GetNeighborRank();
                local_comm_list[neighbor_rank] = 1;
                ++local_memory_size_array[neighbor_rank];
            }
        }

        // sizes are switched bcs transfer direction is inverted from searching to mapping
        MapperUtilitiesMPI::ComputeMaxBufferSizes(local_memory_size_array,
                rMaxReceiveBufferSize,
                rMaxSendBufferSize,
                mCommRank,
                mCommSize);

        MapperUtilitiesMPI::ComputeColoringGraph(local_comm_list, mCommSize,
                rColoredGraph, rMaxColors);

        delete [] local_comm_list;
        delete [] local_memory_size_array;
    }

    void FillBufferWithValues(double* pBuffer, int& rBufferSize, const int CommPartner,
                              const std::function<double(InterfaceObject::Pointer, const std::vector<double>&)>& FunctionPointer) override
    {
        int i = 0;
        std::vector<InterfaceObject::Pointer> interface_objects;
        if (mReceiveObjects.count(CommPartner) > 0)
        {
            interface_objects = mReceiveObjects.at(CommPartner);
        }

        for (auto interface_obj : interface_objects)
        {
            pBuffer[i] = FunctionPointer(interface_obj, mShapeFunctionValues.at(CommPartner)[i]);
            ++i;
        }

        rBufferSize = static_cast<int>(interface_objects.size());

        KRATOS_DEBUG_ERROR_IF_NOT(rBufferSize == i) << "size mismatch" << std::endl;
    }

    void FillBufferWithValues(double* pBuffer, int& rBufferSize, const int CommPartner,
                              const std::function<array_1d<double, 3>(InterfaceObject::Pointer, const std::vector<double>&)>& FunctionPointer) override
    {
        int i = 0;
        std::vector<InterfaceObject::Pointer> interface_objects;
        if (mReceiveObjects.count(CommPartner) > 0)
        {
            interface_objects = mReceiveObjects.at(CommPartner);
        }

        array_1d<double, 3> value;

        for (auto interface_obj : interface_objects)
        {
            value = FunctionPointer(interface_obj, mShapeFunctionValues.at(CommPartner)[i]);

            pBuffer[(i * 3) + 0] = value[0];
            pBuffer[(i * 3) + 1] = value[1];
            pBuffer[(i * 3) + 2] = value[2];

            ++i;
        }

        rBufferSize = static_cast<int>(interface_objects.size()) * 3;

        KRATOS_DEBUG_ERROR_IF_NOT(rBufferSize == i * 3) << "size mismatch" << std::endl;
    }

    void ProcessValues(const double* pBuffer, const int BufferSize, const int CommPartner,
                       const std::function<void(InterfaceObject::Pointer, double)>& FunctionPointer) override
    {
        std::vector<InterfaceObject::Pointer> interface_objects;
        if (mSendObjects.count(CommPartner) > 0)
        {
            interface_objects = mSendObjects.at(CommPartner);
        }

        KRATOS_DEBUG_ERROR_IF_NOT(static_cast<int>(interface_objects.size()) == BufferSize)
            << "Wrong number of results received!; "
            << "interface_objects.size() = " << interface_objects.size()
            << ", BufferSize = " << BufferSize << std::endl;

        for (int i = 0; i < BufferSize; ++i)
        {
            FunctionPointer(interface_objects[i], pBuffer[i]);
        }
    }

    void ProcessValues(const double* pBuffer, const int BufferSize, const int CommPartner,
                       const std::function<void(InterfaceObject::Pointer, array_1d<double, 3>)>& FunctionPointer) override
    {
        KRATOS_DEBUG_ERROR_IF_NOT(BufferSize % 3 == 0)
            << "Uneven number of results "
            << "received!; BufferSize modulo 3 = "
            << BufferSize % 3 << std::endl;

        const int num_values = BufferSize / 3;

        std::vector<InterfaceObject::Pointer> interface_objects;
        if (mSendObjects.count(CommPartner) > 0)
        {
            interface_objects = mSendObjects.at(CommPartner);
        }

        KRATOS_DEBUG_ERROR_IF_NOT(static_cast<int>(interface_objects.size()) == num_values)
            << "Wrong number of results received!; "
            << "interface_objects.size() = "
            << interface_objects.size() << ", num_values = "
            << num_values << std::endl;

        array_1d<double, 3> value;

        for (int i = 0; i < num_values; ++i)
        {
            value[0] = pBuffer[(i * 3) + 0];
            value[1] = pBuffer[(i * 3) + 1];
            value[2] = pBuffer[(i * 3) + 2];

            FunctionPointer(interface_objects[i], value);
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
    virtual std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "InterfaceObjectManagerParallel" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "InterfaceObjectManagerParallel";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override {}


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
    InterfaceObjectManagerParallel& operator=(InterfaceObjectManagerParallel const& rOther);

    //   /// Copy constructor.
    //   InterfaceObjectManagerParallel(InterfaceObjectManagerParallel const& rOther){}


    ///@}

}; // Class InterfaceObjectManagerParallel

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  InterfaceObjectManagerParallel& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const InterfaceObjectManagerParallel& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_INTERFACE_OBJECT_MANAGER_PARALLEL_H_INCLUDED  defined
