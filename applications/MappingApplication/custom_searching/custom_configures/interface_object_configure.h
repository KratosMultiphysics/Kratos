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

#if !defined(KRATOS_INTERFACE_OBJECT_CONFIGURE_H_INCLUDED)
#define  KRATOS_INTERFACE_OBJECT_CONFIGURE_H_INCLUDED

// System includes

// Kratos includes
#include "utilities/contact_pair.h"
#include "custom_searching/interface_object.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/** Configuration file for InterfaceObject.
 * This class provides a configuration file to calculate a 'Bins'
 * using InterfaceObject.
 * It is an adapted copy of "kratos/spatial_containers/configures/point_configure.h"
 */
class InterfaceObjectConfigure
{
public:

    /// Pointer definition of InterfaceObjectConfigure
    KRATOS_CLASS_POINTER_DEFINITION(InterfaceObjectConfigure);

    /** Compile time definitions
     * @param Epsilon   Error tolerance for cmparison operations with doubles
     * @param Dimension Dimension of the problem. Fixed to 3.
     */
    static constexpr auto Epsilon   = std::numeric_limits<double>::epsilon();
    static constexpr auto Dimension = 3;

    /** Point and Pointer Types
     * @param PointType   Point of doubles with 3 coordinates (Dimension = 3)
     * @param PointerType Pointer to Point of doubles with 3 coordinates (Dimension = 3)
     */
    typedef Point          PointType;

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
    typedef InterfaceObject                     ObjectType;
    typedef PointerVectorSet <
    ObjectType,
    IndexedObject
    >                                           ObjectContainerType;
    typedef ObjectType::Pointer                 PointerType;

    typedef ObjectContainerType::ContainerType  ContainerType;
    typedef ObjectContainerType::ContainerType  ResultContainerType;
    typedef ContactPair<PointerType>            ContactPairType;

    typedef ContainerType::iterator             IteratorType;
    typedef ResultContainerType::iterator       ResultIteratorType;
    typedef std::vector<double>::iterator       DistanceIteratorType;

    typedef std::vector<ContactPairType>        ContainerContactType;
    typedef ContainerContactType::iterator      IteratorContactType;

    typedef double                              CoordinateType;
    typedef Tvector<CoordinateType, Dimension>   CoordinateArray;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default consturctor
    InterfaceObjectConfigure() {};

    /// Default destructor
    virtual ~InterfaceObjectConfigure() {}

    ///@}
    ///@name Operations
    ///@{

    /** Calculates the bounding box for the given object.
     * For this configuation file, the bounding box is the equal to the point given in 'rObject'.
     * @param rObject    Point for which the bounding box will be calculated.
     * @param rLowPoint  Lower point of the boundingbox.
     * @param rHighPoint Higher point of the boundingbox.
     */
    static inline void CalculateBoundingBox(const PointerType& rObject, PointType& rLowPoint, PointType& rHighPoint)
    {
        rHighPoint = rLowPoint = *rObject;
    }

    /** Calculates the bounding box for the given object extended with a Radius.
     * For this configuation file, the bounding box is the equal to the point given in 'rObject' + - a radius.
     * @param rObject    Point for which the bounding box will be calculated.
     * @param rLowPoint  Lower point of the boundingbox.
     * @param rHighPoint Higher point of the boundingbox.
     * @param Radius     The extension radius to be applied to the boundingbox.
     */
    static inline void CalculateBoundingBox(const PointerType& rObject, PointType& rLowPoint, PointType& rHighPoint, const double& Radius)
    {
        auto radiusExtension = PointType(Radius, Radius, Radius);

        rLowPoint  = *rObject - radiusExtension;
        rHighPoint = *rObject + radiusExtension;
    }

    /** Calculates the Center of the object.
     * @param rObject        Point for which the bounding box will be calculated.
     * @param rCentralPoint  The center point of the object.
     */
    static inline void CalculateCenter(const PointerType& rObject, PointType& rCentralPoint)
    {
        rCentralPoint = *rObject;
    }

    /** Tests the intersection of two objects
     * For this configuation file, tests if the two points are the same within a Epsilon tolerance range.
     * @param  rObj_1 First point of the tests
     * @param  rObj_2 Second point of the tests
     * @return        Boolean indicating the result of the intersection test described.
     */
    static inline bool Intersection(const PointerType& rObj_1, const PointerType& rObj_2)
    {
        for(std::size_t i = 0; i < Dimension; i++)
        {
            if(std::fabs((*rObj_1)[i] - (*rObj_2)[i]) > Epsilon)
            {
                return false;
            }
        }

        return true;
    }

    /** Tests the intersection of two objects extended with a given radius.
     * For this configuation file, tests if the two points extended with a radius
     * are the same within a Epsilon tolerance range.
     * @param  rObj_1 First point of the tests
     * @param  rObj_2 Second point of the tests
     * @param  Radius The extension radius to be applied in the intersection.
     * @return        Boolean indicating the result of the intersection test described.
     */
    static inline bool Intersection(const PointerType& rObj_1, const PointerType& rObj_2, double Radius)
    {
        // TODO change to squared distance function
        // for(std::size_t i = 0; i < Dimension; i++) {
        //   if(std::fabs((*rObj_1)[i] - (*rObj_2)[i]) > Epsilon + Radius) {
        //     return false;
        //   }
        // }
        //
        // return true;

        double pwdDistance = 0.0f;

        for(std::size_t i = 0; i < Dimension; i++)
        {
            pwdDistance += std::pow((*rObj_1)[i] - (*rObj_2)[i], 2);
        }

        if (std::sqrt(pwdDistance) > Epsilon + Radius)
        {
            return false;
        }
        else
        {
            return true;
        }
    }

    /** Tests the intersection of one object with a boundingbox descrived by 'rLowPoint' and 'rHighPoint'.
     * For this configuation file, tests if one point is inside the boundingbox
     * described by 'rLowPoint' and 'rHighPoint' within a Epsilon tolerance range.
     * @param  rObject    Point of the tests.
     * @param  rLowPoint  Lower point of the boundingbox.
     * @param  rHighPoint Higher point of the boundingbox.
     * @return            Boolean indicating the result of the intersection test described.
     */
    static inline bool IntersectionBox(const PointerType& rObject, const PointType& rLowPoint, const PointType& rHighPoint)
    {
        for(std::size_t i = 0; i < Dimension; i++)
        {
            if( (*rObject)[i] < rLowPoint[i] - Epsilon || (*rObject)[i] > rHighPoint[i] + Epsilon)
            {
                return false;
            }
        }

        return true;
    }

    /** Tests the intersection of one object with a boundingbox descrived by 'rLowPoint' and 'rHighPoint'.
     * For this configuation file, tests if one point extended by radius is inside the boundingbox
     * described by 'rLowPoint' and 'rHighPoint' within a Epsilon tolerance range.
     * @param  rObject    Point of the tests.
     * @param  rLowPoint  Lower point of the boundingbox.
     * @param  rHighPoint Higher point of the boundingbox.
     * @param  Radius     The extension radius to be applied in the intersection.
     * @return            Boolean indicating the result of the intersection test described.
     */
    static inline bool IntersectionBox(const PointerType& rObject, const PointType& rLowPoint, const PointType& rHighPoint, const double& Radius)
    {
        for(std::size_t i = 0; i < Dimension; i++)
        {
            if( ((*rObject)[i] + Radius) < rLowPoint[i] - Epsilon || ((*rObject)[i] - Radius) > rHighPoint[i] + Epsilon)
            {
                return false;
            }
        }

        return true;
    }

    /** Calculates the distance betwen two objects.
     * For this configuation file, calculates the euclidean distance between 'rObj_1' and 'rObj_2'.
     * # Performance
     * In C++11 'std::pow(T, int)' provides the optimal solution in terms of speed.
     * # References
     * (http://en.cppreference.com/w/cpp/numeric/math/pow)
     * (http://stackoverflow.com/questions/2940367)
     * @param rObj_1      First point.
     * @param rLowPoint   Lower point.
     * @param rHighPoint  Higher point of the boundingbox.
     * @param distance    The euclidean distance between 'rObj_1' and 'rObj_2'.
     */
    static inline void Distance(const PointerType& rObj_1, const PointerType& rObj_2, double& distance)
    {
        double pwdDistance = 0.0f;

        for(std::size_t i = 0; i < Dimension; i++)
        {
            pwdDistance += std::pow((*rObj_1)[i] - (*rObj_2)[i], 2);
        }

        distance = std::sqrt(pwdDistance);
    }

    /** Returns a radius associated to the object
     * Returns a radius associated to the object
     * @param  rObject the object
     * @param  Radius  an extension factor.
     * @return         0.0f always.
     */
    static inline double GetObjectRadius(const PointerType& rObject, const double& Radius)
    {
        return 0.0f;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turns back information as a string.
    virtual std::string Info() const
    {
        return "Spatial Containers Configure for 'Points'";
    }

    /// Turns back data as a string.
    virtual std::string Data() const
    {
        return "Dimension: " + std::to_string(Dimension);
    }

    /// Prints object's information.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info() << std::endl;
    }

    /// Prints object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << Data() << Dimension << std::endl;
    }

    ///@}

private:
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    InterfaceObjectConfigure& operator=(InterfaceObjectConfigure const& rOther);

    /// Copy constructor.
    InterfaceObjectConfigure(InterfaceObjectConfigure const& rOther);

    ///@}

}; // Class InterfaceObjectConfigure

///@}

///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream, InterfaceObjectConfigure& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const InterfaceObjectConfigure& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

} // namespace Kratos.
#endif /* KRATOS_INTERFACE_OBJECT_CONFIGURE_H_INCLUDED */
