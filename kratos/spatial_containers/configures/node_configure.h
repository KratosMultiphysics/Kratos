//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos A. Roig
//

#if !defined(KRATOS_NODE_CONFIGURE_INCLUDED)
#define  KRATOS_NODE_CONFIGURE_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <limits>
#include <cmath>


#include "spatial_containers/tree.h"
#include "spatial_containers/cell.h"

// Kratos includes
#include "includes/define.h"
#include "geometries/point.h"
#include "containers/pointer_vector_set.h"
#include "utilities/indexed_object.h"
#include "utilities/contact_pair.h"

namespace Kratos {

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

/** @brief Configuration file for Nodes.
 * @details This class provides a configuration file to calculate a 'Bins'
 * using nodes.
 */
class NodeConfigure {
public:
    /// Pointer definition of NodeConfigure
    KRATOS_CLASS_POINTER_DEFINITION(NodeConfigure);

    /** @brief Compile time definitions
    * @param Epsilon Error tolerance for cmparison operations with doubles
    * @param Dimension Dimension of the problem. Fixed to 3.
    */
    static constexpr auto Epsilon   = std::numeric_limits<double>::epsilon();
    static constexpr auto Dimension = 3;

    /** @brief Point and Pointer Types
    * @param PointType Point of doubles with 3 coordinates (Dimension = 3)
    * @param PointerType Pointer to nodes of doubles with 3 coordinates (Dimension = 3)
    */
    typedef Point                   PointType;
    typedef Node<3>                 ObjectType;
    typedef ObjectType::Pointer     PointerType;

    /** Additional types needed by the bins.
    * @param PointContainerType    Point Container.
    * @param ContainerType         Base container Type.
    * @param ResultContainerType   Result Container. For this configure should be the same as ContainerType.
    * @param ContactPairType       Contact pair for points.
    * @param IteratorType          Iterator of points.
    * @param ResultIteratorType    Iterator of results. For this configure should be the same as PointIteratorType.
    * @param DistanceIteratorType  Iterato of distances (doubles)
    * @param ContainerContactType  Container type for contacts
    * @param IteratorContactType   Iterator type for contacts
    */

    typedef PointerVectorSet<ObjectType,IndexedObject>  ObjectContainerType;

    typedef ObjectContainerType::ContainerType          ContainerType;
    typedef ObjectContainerType::ContainerType          ResultContainerType;
    typedef ContactPair<PointerType>                    ContactPairType;

    typedef ContainerType::iterator                     IteratorType;
    typedef ResultContainerType::iterator               ResultIteratorType;
    typedef std::vector<double>::iterator               DistanceIteratorType;

    typedef std::vector<ContactPairType>                ContainerContactType;
    typedef ContainerContactType::iterator              IteratorContactType;

    typedef double                                      CoordinateType;
    typedef Tvector<CoordinateType,Dimension>           CoordinateArray;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default consturctor
    NodeConfigure(){};

    /// Default destructor
    virtual ~NodeConfigure(){}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /** @brief Calculates the bounding box for the given object.
    * @details For this configuation file, the bounding box is the equal to the point given in 'rObject'.
    * @param rObject Point for which the bounding box will be calculated.
    * @param rLowPoint Lower point of the boundingbox.
    * @param rHighPoint Higher point of the boundingbox.
    */
    static inline void CalculateBoundingBox(
        const PointerType& rObject, 
        PointType& rLowPoint, 
        PointType& rHighPoint
        ) 
    {
        rHighPoint = rLowPoint = *rObject;
    }

    /** @brief Calculates the bounding box for the given object extended with a Radius.
    * @details For this configuation file, the bounding box is the equal to the point given in 'rObject' + - a radius.
    * @param rObject Point for which the bounding box will be calculated.
    * @param rLowPoint Lower point of the boundingbox.
    * @param rHighPoint Higher point of the boundingbox.
    * @param Radius The extension radius to be applied to the boundingbox.
    */
    static inline void CalculateBoundingBox(
        const PointerType& rObject, 
        PointType& rLowPoint, 
        PointType& rHighPoint, 
        const double Radius
        ) 
    {
        auto radiusExtension = PointType(Radius, Radius, Radius);

        rLowPoint  = *rObject - radiusExtension;
        rHighPoint = *rObject + radiusExtension;
    }

    /** @brief Calculates the Center of the object.
    * @param rObject Point for which the bounding box will be calculated.
    * @param rCentralPoint The center point of the object.
    */
    static inline void CalculateCenter(
        const PointerType& rObject, 
        PointType& rCentralPoint
        ) 
    {
        rCentralPoint = *rObject;
    }

    /** @brief Tests the intersection of two objects
    * @details For this configuation file, tests if the two points are the same within a Epsilon tolerance range.
    * @param rObj_1 First point of the tests
    * @param rObj_2 Second point of the tests
    * @return Boolean indicating the result of the intersection test described.
    */
    static inline bool Intersection(
        const PointerType& rObj_1, 
        const PointerType& rObj_2
        ) 
    {
        for(std::size_t i = 0; i < Dimension; i++) {
            if(std::fabs((*rObj_1)[i] - (*rObj_2)[i]) > Epsilon) {
                return false;
            }
        }

        return true;
    }

    /** @brief Tests the intersection of two objects extended with a given radius.
    * @details For this configuation file, tests if the two points extended with a radius
    * are the same within a Epsilon tolerance range.
    * @param rObj_1 First point of the tests
    * @param rObj_2 Second point of the tests
    * @param Radius The extension radius to be applied in the intersection.
    * @return Boolean indicating the result of the intersection test described.
    */
    static inline bool Intersection(
        const PointerType& rObj_1, 
        const PointerType& rObj_2, 
        double Radius
    ) 
    {
        for(std::size_t i = 0; i < Dimension; i++) {
            if(std::fabs((*rObj_1)[i] - (*rObj_2)[i]) > Epsilon + Radius) {
                return false;
            }
        }

        return true;
    }

    /** @brief Tests the intersection of one object with a boundingbox descrived by 'rLowPoint' and 'rHighPoint'.
    * @details For this configuation file, tests if one point is inside the boundingbox
    * described by 'rLowPoint' and 'rHighPoint' within a Epsilon tolerance range.
    * @param rObject Point of the tests.
    * @param rLowPoint Lower point of the boundingbox.
    * @param rHighPoint Higher point of the boundingbox.
    * @return Boolean indicating the result of the intersection test described.
    */
    static inline bool IntersectionBox(
        const PointerType& rObject, 
        const PointType& rLowPoint, 
        const PointType& rHighPoint
        ) 
    {
        for(std::size_t i = 0; i < Dimension; i++) {
            if( (*rObject)[i] < rLowPoint[i] - Epsilon || (*rObject)[i] > rHighPoint[i] + Epsilon) {
                return false;
            }
        }

        return true;
    }

    /** @brief Tests the intersection of one object with a boundingbox descrived by 'rLowPoint' and 'rHighPoint'.
    * @details For this configuation file, tests if one point extended by radius is inside the boundingbox
    * described by 'rLowPoint' and 'rHighPoint' within a Epsilon tolerance range.
    * @param rObject Point of the tests.
    * @param rLowPoint Lower point of the boundingbox.
    * @param rHighPoint Higher point of the boundingbox.
    * @param Radius The extension radius to be applied in the intersection.
    * @return Boolean indicating the result of the intersection test described.
    */
    static inline bool IntersectionBox(
        const PointerType& rObject, 
        const PointType& rLowPoint, 
        const PointType& rHighPoint, 
        const double Radius
        ) 
    {
        for(std::size_t i = 0; i < Dimension; i++) {
            if( ((*rObject)[i] + Radius) < rLowPoint[i] - Epsilon || ((*rObject)[i] - Radius) > rHighPoint[i] + Epsilon) {
                return false;
            }
        }

        return true;
    }

    /** @brief Calculates the distance betwen two objects.
    * @details For this configuation file, calculates the euclidean distance between 'rObj_1' and 'rObj_2'.
    * -# Performance
    *   - In C++11 'std::pow(T, int)' provides the optimal solution in terms of speed.
    * -# References
    *  - <a href="http://en.cppreference.com/w/cpp/numeric/math/pow">CPP References. 'std::pow' method </a>
    *  - <a href="http://stackoverflow.com/questions/2940367">Stackoverflow</a>
    * @param rObj_1 First point.
    * @param rObj_2 Second point
    * @param distance The euclidean distance between 'rObj_1' and 'rObj_2'.
    */
    static inline void Distance(
        const PointerType& rObj_1, 
        const PointerType& rObj_2, 
        double& distance
        ) 
    {
        double pwdDistance = 0.0f;

        for(std::size_t i = 0; i < Dimension; i++) {
            pwdDistance += std::pow((*rObj_1)[i] - (*rObj_2)[i], 2);
        }

        distance = std::sqrt(pwdDistance);
    }

    /** @brief Returns a radius associated to the object
    * @details Returns a radius associated to the object
    * @param rObject The object
    * @param Radius An extension factor.
    * @return 0.0f always.
    */
    static inline double GetObjectRadius(
        const PointerType& rObject, 
        const double Radius
        ) 
    {
        return 0.0f;
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

    /// Turns back information as a string.
    virtual std::string Info() const {
        return "Spatial Containers Configure for 'Points'";
    }

    /// Turns back data as a string.
    virtual std::string Data() const {
        return "Dimension: " + std::to_string(Dimension);
    }

    /// Prints object's information.
    virtual void PrintInfo(std::ostream& rOStream) const {
        rOStream << Info() << std::endl;
    }

    /// Prints object's data.
    virtual void PrintData(std::ostream& rOStream) const {
        rOStream << Data() << Dimension << std::endl;
    }

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
    NodeConfigure& operator=(NodeConfigure const& rOther);

    /// Copy constructor.
    NodeConfigure(NodeConfigure const& rOther);

    ///@}

}; // Class NodeConfigure

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream, NodeConfigure& rThis){
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const NodeConfigure& rThis){
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

  return rOStream;
}

///@}

} // namespace Kratos.
#endif /* POINT_CONFIGURE */
