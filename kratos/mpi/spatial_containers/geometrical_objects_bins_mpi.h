//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes

// External includes

// Project includes
#include "spatial_containers/geometrical_objects_bins.h"

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

class GeometricalObject; // forward declaration, to be included in the cpp. This is needed to reduce the compilation time. Can be done as we consider the GeometricalObject as a pointer
class DataCommunicator;  // forward declaration, to be included in the cpp. This is needed to reduce the compilation time. Can be done as we consider the GeometricalObject as a pointer

/**
 * @class GeometricalObjectsBinsMPI
 * @ingroup KratosCore
 * @brief A bins container for 3 dimensional GeometricalObject entities (MPI version)
 * @details This is the MPI version of the GeometricalObjectsBins, which is a container for geometrical objects. It is used to perform fast search of geometrical objects in a given space.
 * @author Vicente Mataix Ferrandiz
*/
class KRATOS_API(KRATOS_CORE) GeometricalObjectsBinsMPI
    : public GeometricalObjectsBins
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of GeometricalObjectsBinsMPI
    KRATOS_CLASS_POINTER_DEFINITION(GeometricalObjectsBinsMPI);

    /// The base type
    using BaseType = GeometricalObjectsBins;

    /// The buffer type definition
    using BufferTypeDouble = std::vector<std::vector<double>>;
    using BufferTypeChar = std::vector<std::vector<char>>;

    /// The type of geometrical object to be stored in the bins
    using BaseType::CellType;
    using BaseType::ResultType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor deleted.
    GeometricalObjectsBinsMPI() = delete;

    /**
     * @brief The constructor with all geometries to be stored. Please note that all of them should be available at construction time and cannot be modified after.
     * @param GeometricalObjectsBegin The begin iterator of the geometries to be stored
     * @param GeometricalObjectsEnd The end iterator of the geometries to be stored
     * @param rDataCommunicator The data communicator
     * @tparam TIteratorType The type of the iterator
     */
    template<typename TIteratorType>
    GeometricalObjectsBinsMPI(
        TIteratorType GeometricalObjectsBegin,
        TIteratorType GeometricalObjectsEnd,
        const DataCommunicator& rDataCommunicator
        ) : mLocalGeometricalObjectsBins(GeometricalObjectsBegin, GeometricalObjectsEnd),
            mrDataCommunicator(rDataCommunicator)
    {
        // We get the world size
        const int world_size = GetWorldSize();

        // Set up the buffers
        mSendSizes.resize(world_size);
        mRecvSizes.resize(world_size);

        mSendBufferDouble.resize(world_size);
        mRecvBufferDouble.resize(world_size);

        mSendBufferChar.resize(world_size);
        mRecvBufferChar.resize(world_size);

        // Set up the global bounding boxes
        InitializeGlobalBoundingBoxes();
    }

    /**
     * @brief The constructor with all geometries to be stored. Please note that all of them should be available at construction time and cannot be modified after.
     * @param rGeometricalObjectsVector The geometries to be stored
     * @param rDataCommunicator The data communicator
     * @tparam TContainer The container type
     */
    template<typename TContainer>
    GeometricalObjectsBinsMPI(
        TContainer& rGeometricalObjectsVector,
        const DataCommunicator& rDataCommunicator
        ) : GeometricalObjectsBinsMPI(rGeometricalObjectsVector.begin(), rGeometricalObjectsVector.end(), rDataCommunicator)
    {
    }

    /// Destructor.
    ~GeometricalObjectsBinsMPI() override {}

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
        buffer << "GeometricalObjectsBinsMPI" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {rOStream << "GeometricalObjectsBinsMPI";}

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {}

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

    std::vector<double> mGlobalBoundingBoxes; /// All the global BB, data is xmax, xmin,  ymax, ymin,  zmax, zmin

    GeometricalObjectsBins mLocalGeometricalObjectsBins; /// The local bins

    const DataCommunicator& mrDataCommunicator; /// The data communicator

    std::vector<int> mSendSizes; /// The sizes of the send buffers
    std::vector<int> mRecvSizes; /// The sizes of the recv buffers

    BufferTypeDouble mSendBufferDouble; /// The send buffer (double)
    BufferTypeDouble mRecvBufferDouble; /// The recv buffer (double)

    BufferTypeChar mSendBufferChar; /// The send buffer (char)
    BufferTypeChar mRecvBufferChar; /// The recv buffer (char)

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Returns the current rank
     * @return The current rank
     */
    int GetRank() const;

    /**
     * @brief Returns the world size
     * @return The world size
     */
    int GetWorldSize() const;

    /**
     * @brief Initializes the global bounding boxes
     */
    void InitializeGlobalBoundingBoxes();

    /**
     * @brief This method checks if a point is inside any bounding box of the global bounding boxes
     * @param rCoords The coordinates of the point
     * @return The index of the ranks inside the bounding box
     */
    std::vector<int> RansksPointIsInsideBoundingBox(const array_1d<double, 3>& rCoords);

    // /**
    //  * @brief This method synchronizes the search in radius
    //  * @param rResults The results of the search (local version)
    //  */
    // void SynchronizeSearchInRadius(std::vector<ResultType>& rLocalResults);

    // /**
    //  * @brief This method synchronizes the search nearest in radius
    //  * @param rLocalResult The results of the search (local version)
    //  * @return The results of the search (global version)
    //  */
    // ResultType SynchronizeSearchNearestInRadius(ResultType& rLocalResult);

    // /**
    //  * @brief This method synchronizes the search nearest
    //  * @param rLocalResult The results of the search (local version)
    //  * @return The results of the search (global version)
    //  */
    // ResultType SynchronizeSearchNearest(ResultType& rLocalResult);

    // /**
    //  * @brief This method synchronizes the search is inside
    //  * @param rLocalResult The results of the search (local version)
    //  * @return The results of the search (global version)
    //  */
    // ResultType SynchronizeSearchIsInside(ResultType& rLocalResult);

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator deleted.
    GeometricalObjectsBinsMPI& operator=(GeometricalObjectsBinsMPI const& rOther) = delete;

    /// Copy constructor deleted.
    GeometricalObjectsBinsMPI(GeometricalObjectsBinsMPI const& rOther) = delete;

    ///@}

}; // Class GeometricalObjectsBinsMPI

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@} addtogroup block

}  // namespace Kratos.