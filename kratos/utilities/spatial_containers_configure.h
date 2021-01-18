//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Nelson Lafontaine
//   Last update:    Vicente Mataix Ferrandiz
//                    
//

#if !defined(KRATOS_SPATIAL_CONTAINERS_CONFIGURE_INCLUDED )
#define  KRATOS_SPATIAL_CONTAINERS_CONFIGURE_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <cmath>

// Project includes
#include "utilities/contact_pair.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    typedef std::size_t SizeType;
    
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
    typedef Point                                                                      PointType; 
    
    /// The node definition
    typedef Node<3>                                                                     NodeType;
    
    /// The geometry definition
    typedef Geometry<NodeType>                                                      GeometryType;
    
    /// Dditance iterator
    typedef std::vector<double>::iterator                                   DistanceIteratorType;
    
    /// The entity definition
    typedef TEntity                                                                   EntityType;
    
    /// Container definition
    typedef typename PointerVectorSet<TEntity, IndexedObject>::ContainerType       ContainerType;
    typedef typename ContainerType::value_type                                       PointerType;
    typedef typename ContainerType::iterator                                        IteratorType;
    typedef typename PointerVectorSet<TEntity, IndexedObject>::ContainerType ResultContainerType;
    typedef typename ResultContainerType::value_type                           ResultPointerType;
    typedef typename ResultContainerType::iterator                            ResultIteratorType;

    /// Contact Pairs
    typedef ContactPair<PointerType>                                             ContactPairType;
//     typedef array_1d<PointerType, 2>                                          ContactPairType;
    typedef std::vector<ContactPairType>                                    ContainerContactType;
    typedef typename ContainerContactType::iterator                          IteratorContactType;
    typedef typename ContainerContactType::value_type                         PointerContactType;
    typedef typename std::vector<PointerType>::iterator                      PointerTypeIterator;

    /// Pointer definition of SpatialContainersConfigure
    KRATOS_CLASS_POINTER_DEFINITION(SpatialContainersConfigure);

    ///@}
    ///@name  Enum's
    ///@{

    enum { Dimension = TDimension,
           DIMENSION = TDimension,
           MAX_LEVEL = 16,
           MIN_LEVEL = 2
         };
    
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
     */
    static inline void CalculateBoundingBox(const PointerType& rObject, PointType& rLowPoint, PointType& rHighPoint)
    {
        rHighPoint = rObject->GetGeometry().GetPoint(0);
        rLowPoint  = rObject->GetGeometry().GetPoint(0);
        for (unsigned int point = 0; point<rObject->GetGeometry().PointsNumber(); point++)
        {
            for(std::size_t i = 0; i<TDimension; i++)
            {
                rLowPoint[i]  =  (rLowPoint[i]  >  rObject->GetGeometry().GetPoint(point)[i] ) ?  rObject->GetGeometry().GetPoint(point)[i] : rLowPoint[i];
                rHighPoint[i] =  (rHighPoint[i] <  rObject->GetGeometry().GetPoint(point)[i] ) ?  rObject->GetGeometry().GetPoint(point)[i] : rHighPoint[i];
            }
        }
    }

    /**
     * @brief This method computes the bounding box
     */
    static inline void CalculateBoundingBox(const PointerType& rObject, PointType& rLowPoint, PointType& rHighPoint, const double Radius)
    {
        (void)Radius;
        CalculateBoundingBox(rObject, rLowPoint, rHighPoint);
    }

    /**
     * @brief This method computes the intersection
     */
    static inline bool Intersection(const PointerType& rObj_1, const PointerType& rObj_2)
    {
        GeometryType& geom_1 = rObj_1->GetGeometry();
        GeometryType& geom_2 = rObj_2->GetGeometry();
        return  geom_1.HasIntersection(geom_2);

    }

    /**
     * @brief This method computes the intersection
     */
    static inline bool Intersection(const PointerType& rObj_1, const PointerType& rObj_2, const double Radius)
    {
        (void)Radius;
        return Intersection(rObj_1, rObj_2);
    }

    /**
     * @brief This method computes the intersection box
     */
    static inline bool IntersectionBox(const PointerType& rObject,  const PointType& rLowPoint, const PointType& rHighPoint)
    {
        return rObject->GetGeometry().HasIntersection(rLowPoint, rHighPoint);
    }

    /**
     * @brief This method computes the intersection box
     */
    static inline bool  IntersectionBox(const PointerType& rObject,  const PointType& rLowPoint, const PointType& rHighPoint, const double Radius)
    {
        (void)Radius;
        return IntersectionBox(rObject, rLowPoint, rHighPoint);
    }

    /**
     * @brief This method computes the distance
     */
    static inline void Distance(const PointerType& rObj_1, const PointerType& rObj_2, double& distance) {}


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

#endif // KRATOS_SPATIAL_CONTAINERS_CONFIGURE_INCLUDED  defined
