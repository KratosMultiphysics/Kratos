//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nelson Lafontaine
//                   Vicente Mataix Ferrandiz
//

#pragma once

// System includes
#include <string>
#include <iostream>
#include <cmath>

// External includes

// Project includes

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

using SizeType = std::size_t;

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class SpatialContainersConfigure
 * @ingroup KratosCore
 * @brief Thhis class is a container for spatial search
 * @details It is used in binbased locator among other classes and utilities
 * @tparam TDimension The working dimension
 * @tparam TEntity The entity considered
 * @author Nelson Lafontaine
*/
template <SizeType TDimension, class TEntity = Element>
class SpatialContainersConfigure
{
public:
    ///@name Type Definitions
    ///@{

    /// Point definition
    using PointType = Point;

    /// The node definition
    using NodeType = Node;

    /// The geometry definition
    using GeometryType = Geometry<NodeType>;

    /// Dditance iterator
    using DistanceIteratorType = std::vector<double>::iterator;

    /// The entity definition
    using EntityType = TEntity;

    /// Container definition
    using ContainerType = typename PointerVectorSet<TEntity, IndexedObject>::ContainerType;
    using PointerType = typename ContainerType::value_type;
    using IteratorType = typename ContainerType::iterator;
    using ResultContainerType = typename PointerVectorSet<TEntity, IndexedObject>::ContainerType;
    using ResultPointerType = typename ResultContainerType::value_type;
    using ResultIteratorType = typename ResultContainerType::iterator;

    /// Definition of the Dimension (it is the template argument, but it is needed to be defined as a static member for legacy reasons as it is used that way in other places)
    static constexpr std::size_t Dimension = TDimension;

    /// Definition of the DIMENSION (it is the template argument, but it is needed to be defined as a static member for legacy reasons as it is used that way in other places)
    static constexpr std::size_t DIMENSION = TDimension;

    /// Definition of the maximum level
    static constexpr std::size_t MAX_LEVEL = 16;

    /// Definition of the minimum level
    static constexpr std::size_t MIN_LEVEL = 2;

    /// Pointer definition of SpatialContainersConfigure
    KRATOS_CLASS_POINTER_DEFINITION(SpatialContainersConfigure);

    ///@}
    ///@name  Enum's
    ///@{

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SpatialContainersConfigure() {}

    /// Destructor.
    virtual ~SpatialContainersConfigure() {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method computes the bounding box
     * @param rObject The object to compute the bounding box
     * @param rLowPoint The lower point of the bounding box
     * @param rHighPoint The higher point of the bounding box
     */
    static inline void CalculateBoundingBox(
        const PointerType& rObject,
        PointType& rLowPoint,
        PointType& rHighPoint
        )
    {
        rHighPoint = rObject->GetGeometry().GetPoint(0);
        rLowPoint  = rObject->GetGeometry().GetPoint(0);
        for (unsigned int point = 0; point<rObject->GetGeometry().PointsNumber(); point++) {
            for(unsigned int i = 0; i<TDimension; i++) {
                rLowPoint[i]  =  (rLowPoint[i]  >  rObject->GetGeometry().GetPoint(point)[i] ) ?  rObject->GetGeometry().GetPoint(point)[i] : rLowPoint[i];
                rHighPoint[i] =  (rHighPoint[i] <  rObject->GetGeometry().GetPoint(point)[i] ) ?  rObject->GetGeometry().GetPoint(point)[i] : rHighPoint[i];
            }
        }
    }

    /**
     * @brief This method computes the bounding box
     * @param rObject The object to compute the bounding box
     * @param rLowPoint The lower point of the bounding box
     * @param rHighPoint The higher point of the bounding box
     * @param Radius The radius
     */
    static inline void CalculateBoundingBox(
        const PointerType& rObject,
        PointType& rLowPoint,
        PointType& rHighPoint,
        const double Radius
        )
    {
        (void)Radius;
        CalculateBoundingBox(rObject, rLowPoint, rHighPoint);
    }

    /**
     * @brief This method computes the intersection
     * @param rObj_1 The first object
     * @param rObj_2 The second object
     * @return True if there is an intersection, false otherwise
     */
    static inline bool Intersection(
        const PointerType& rObj_1,
        const PointerType& rObj_2
        )
    {
        GeometryType& r_geom_1 = rObj_1->GetGeometry();
        GeometryType& r_geom_2 = rObj_2->GetGeometry();
        return r_geom_1.HasIntersection(r_geom_2);
    }

    /**
     * @brief This method computes the intersection
     * @param rObj_1 The first object
     * @param rObj_2 The second object
     * @param Radius The radius
     */
    static inline bool Intersection(
        const PointerType& rObj_1,
        const PointerType& rObj_2,
        const double Radius
        )
    {
        (void)Radius;
        return Intersection(rObj_1, rObj_2);
    }

    /**
     * @brief This method computes the intersection box
     * @param rObject The object considered
     * @param rLowPoint The low point of the box
     * @param rHighPoint The high point of the box
     */
    static inline bool IntersectionBox(
        const PointerType& rObject,
        const PointType& rLowPoint,
        const PointType& rHighPoint
        )
    {
        return rObject->GetGeometry().HasIntersection(rLowPoint, rHighPoint);
    }

    /**
     * @brief This method computes the intersection box
     * @param rObject The object considered
     * @param rLowPoint The low point of the box
     * @param rHighPoint The high point of the box
     * @param Radius The radius
     */
    static inline bool IntersectionBox(
        const PointerType& rObject,
        const PointType& rLowPoint,
        const PointType& rHighPoint,
        const double Radius
        )
    {
        (void)Radius;
        return IntersectionBox(rObject, rLowPoint, rHighPoint);
    }

    /**
     * @brief This method computes the distance
     * @param rObj_1 The first object
     * @param rObj_2 The second object
     * @param rDistance The distance
     */
    static inline void Distance(
        const PointerType& rObj_1,
        const PointerType& rObj_2,
        double& rDistance
        )
    {

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
        return " Spatial Containers Configure";
    }

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
    SpatialContainersConfigure& operator=(SpatialContainersConfigure const& rOther);

    /// Copy constructor.
    SpatialContainersConfigure(SpatialContainersConfigure const& rOther);

    ///@}

}; // Class SpatialContainersConfigure

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template <std::size_t TDimension>
inline std::istream& operator >> (std::istream& rIStream,
                                  SpatialContainersConfigure<TDimension> & rThis)
{
    return rIStream;
}

/// output stream function
template <std::size_t TDimension>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const SpatialContainersConfigure<TDimension>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.
