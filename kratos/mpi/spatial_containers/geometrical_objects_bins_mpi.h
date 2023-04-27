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
#include "mpi.h"

// Project includes
#include "mpi/utilities/mpi_search_utilities.h"
#include "spatial_containers/geometrical_objects_bins.h"

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

class GeometricalObject; // forward declaration, to be included in the cpp. This is needed to reduce the compilation time. Can be done as we consider the GeometricalObject as a pointer

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
     * @tparam TIteratorType The type of the iterator
     */
    template<typename TIteratorType>
    GeometricalObjectsBinsMPI(
        TIteratorType GeometricalObjectsBegin,
        TIteratorType GeometricalObjectsEnd
        )
    {
        // Set up the buffers
        MPI_Comm_rank(MPI_COMM_WORLD, &mCommRank);
        MPI_Comm_size(MPI_COMM_WORLD, &mCommSize);

        mSendSizes.resize(mCommSize);
        mRecvSizes.resize(mCommSize);

        mSendBufferDouble.resize(mCommSize);
        mRecvBufferDouble.resize(mCommSize);

        mSendBufferChar.resize(mCommSize);
        mRecvBufferChar.resize(mCommSize);

        // TODO: Do equivalent
        // mMapperInterfaceInfosContainer.resize(mCommSize);

        // We compute the local bounding box
        const std::size_t local_number_of_objects = std::distance(GeometricalObjectsBegin, GeometricalObjectsEnd);
        if (local_number_of_objects > 0){
            mBoundingBox.Set(GeometricalObjectsBegin->GetGeometry().begin(), GeometricalObjectsBegin->GetGeometry().end());
            for (TIteratorType i_object = GeometricalObjectsBegin ; i_object != GeometricalObjectsEnd ; i_object++){
                mBoundingBox.Extend(i_object->GetGeometry().begin() , i_object->GetGeometry().end());
            }
        }
        mBoundingBox.Extend(Tolerance);

        // Now we compute the global bounding boxes
        ComputeGlobalBoundingBoxes();

        // We compute the cell size (WIP)
        // TODO: Update for global bounding boxes
        const std::size_t global_number_of_objects = local_number_of_objects;
        CalculateCellSize(global_number_of_objects);
        mCells.resize(GetTotalNumberOfCells());
        AddObjectsToCells(GeometricalObjectsBegin, GeometricalObjectsEnd);
    }

    /**
     * @brief The constructor with all geometries to be stored. Please note that all of them should be available at construction time and cannot be modified after.
     * @param rGeometricalObjectsVector The geometries to be stored
     * @tparam TContainer The container type
     */
    template<typename TContainer>
    GeometricalObjectsBinsMPI(TContainer& rGeometricalObjectsVector)
        : GeometricalObjectsBinsMPI(rGeometricalObjectsVector.begin(), rGeometricalObjectsVector.end())
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

    // TODO: Replace with the BB class instead
    std::vector<double> mGlobalBoundingBoxes; /// The global bounding boxes of the model part, in the form: xmax, xmin,  ymax, ymin,  zmax, zmin

    double mRadius = 0.0;                     /// The radius of the search

    int mCommRank;
    int mCommSize;

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
     * @brief This method initializes the search
     * @details This method initializes the search
     */
    void InitializeSearch() override;

    /**
     * @brief This method finalizes the search
     * @details This method finalizes the search
     */
    void FinalizeSearch() override;

    /**
     * @brief This method computes the global bounding boxes
     */
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