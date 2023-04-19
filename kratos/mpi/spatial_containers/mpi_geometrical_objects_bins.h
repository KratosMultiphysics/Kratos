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

/**
 * @class MPIGeometricalObjectsBins
 * @ingroup KratosCore
 * @brief A bins container for 3 dimensional GeometricalObject entities (MPI version)
 * @details This is the MPI version of the GeometricalObjectsBins, which is a container for geometrical objects. It is used to perform fast search of geometrical objects in a given space.
 * @author Vicente Mataix Ferrandiz
*/
class KRATOS_API(KRATOS_CORE) MPIGeometricalObjectsBins
    : public GeometricalObjectsBins
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MPIGeometricalObjectsBins
    KRATOS_CLASS_POINTER_DEFINITION(MPIGeometricalObjectsBins);

    /// The base type
    using BaseType = GeometricalObjectsBins;

    /// The type of geometrical object to be stored in the bins
    using BaseType::CellType;
    using BaseType::ResultType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor deleted.
    MPIGeometricalObjectsBins() = delete;

    /**
     * @brief The constructor with all geometries to be stored. Please note that all of them should be available at construction time and cannot be modified after.
     * @param GeometricalObjectsBegin The begin iterator of the geometries to be stored
     * @param GeometricalObjectsEnd The end iterator of the geometries to be stored
     * @tparam TIteratorType The type of the iterator
     */
    template<typename TIteratorType>
    MPIGeometricalObjectsBins(
        TIteratorType GeometricalObjectsBegin,
        TIteratorType GeometricalObjectsEnd
        ) : BaseType(GeometricalObjectsBegin, GeometricalObjectsEnd)
    {
        // DOING NOTHING
    }

    /**
     * @brief The constructor with all geometries to be stored. Please note that all of them should be available at construction time and cannot be modified after.
     * @param rGeometricalObjectsVector The geometries to be stored
     * @tparam TContainer The container type
     */
    template<typename TContainer>
    MPIGeometricalObjectsBins(TContainer& rGeometricalObjectsVector)
        : MPIGeometricalObjectsBins(rGeometricalObjectsVector.begin(), rGeometricalObjectsVector.end())
    {
    }

    /// Destructor.
    ~MPIGeometricalObjectsBins() override {}

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
        buffer << "MPIGeometricalObjectsBins" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {rOStream << "MPIGeometricalObjectsBins";}

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

    /// Assignment operator deleted.
    MPIGeometricalObjectsBins& operator=(MPIGeometricalObjectsBins const& rOther) = delete;

    /// Copy constructor deleted.
    MPIGeometricalObjectsBins(MPIGeometricalObjectsBins const& rOther) = delete;

    ///@}

}; // Class MPIGeometricalObjectsBins

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@} addtogroup block

}  // namespace Kratos.