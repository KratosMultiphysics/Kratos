//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:       Ilaria Iaconeta

//A fast algorithm to check intersection between two axis-aligned quadrilateral is implemented.
//This axis-aligned rect test is a special case of the separating axis theorem.

#if !defined(KRATOS_FAST_QUAD_SPATIAL_CONTAINERS_CONFIGURE_INCLUDED)
#define  KRATOS_FAST_QUAD_SPATIAL_CONTAINERS_CONFIGURE_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <cmath>

// kratos utils
#include "utilities/spatial_containers_configure.h"


namespace Kratos
{

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


template <std::size_t TDimension>
class FastQuadSpatialContainersConfigure
{

public:

    enum
    {
        Dimension = TDimension,
        DIMENSION = TDimension,
        MAX_LEVEL = 16,
        MIN_LEVEL = 2
    };

    /// Pointer definition of FastQuadSpatialContainersConfigure
    KRATOS_CLASS_POINTER_DEFINITION(FastQuadSpatialContainersConfigure);

    typedef Point                                PointType;  /// always the point 3D
    typedef std::vector<double>::iterator                   DistanceIteratorType;
    typedef ModelPart::ElementsContainerType::ContainerType ContainerType;
    typedef ContainerType::value_type                       PointerType;
    typedef ContainerType::iterator                         IteratorType;
    typedef ModelPart::ElementsContainerType::ContainerType ResultContainerType;
    typedef ResultContainerType::value_type                 ResultPointerType;
    typedef ResultContainerType::iterator                   ResultIteratorType;

    /// Contact Pairs
    // typedef std::pair<PointerType, PointerType>            ContactPairType;
    // typedef std::vector<ContactPairType>                   ContainerContactType;
    // typedef std::vector<ContactPairType>::iterator         IteraratorContactType;

    /// Contact Pairs
    typedef ContactPair<PointerType>                        ContactPairType;
    // typedef array_1d<PointerType, 2>                       ContactPairType;
    typedef  std::vector<ContactPairType>                   ContainerContactType;
    typedef  ContainerContactType::iterator                 IteratorContactType;
    typedef  ContainerContactType::value_type               PointerContactType;
    typedef  std::vector<PointerType>::iterator             PointerTypeIterator;


    ///@}
    ///@name Life Cycle
    ///@{

    FastQuadSpatialContainersConfigure() {};
    virtual ~FastQuadSpatialContainersConfigure() {}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

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

    static inline void CalculateBoundingBox(const PointerType& rObject, PointType& rLowPoint, PointType& rHighPoint, const double & Radius)
    {
        KRATOS_ERROR << "Fast Quad spatial does not define this function" << std::endl;
    }

    //******************************************************************************************************************

    static inline bool Intersection(const PointerType& rObj_1, const PointerType& rObj_2)
    {
        Element::GeometryType& geom_1 = rObj_1->GetGeometry();
        Element::GeometryType& geom_2 = rObj_2->GetGeometry();

        PointType rHighPoint1 = geom_1.GetPoint(0);
        PointType rLowPoint1  = geom_1.GetPoint(0);

        PointType rHighPoint2 = geom_2.GetPoint(0);
        PointType rLowPoint2  = geom_2.GetPoint(0);

        //Firstly the Lowpoint and Highpoint are defined for both quadrilaterals

        for (unsigned int point = 1; point<geom_1.PointsNumber(); point++)
        {
            for(std::size_t i = 0; i<TDimension; i++)
            {
                rLowPoint1[i]  =  (rLowPoint1[i]  >  geom_1.GetPoint(point)[i] ) ?  geom_1.GetPoint(point)[i] : rLowPoint1[i];
                rHighPoint1[i] =  (rHighPoint1[i] <  geom_1.GetPoint(point)[i] ) ?  geom_1.GetPoint(point)[i] : rHighPoint1[i];
            }
        }

        for (unsigned int point = 1; point<geom_2.PointsNumber(); point++)
        {
            for(std::size_t i = 0; i<TDimension; i++)
            {
                rLowPoint2[i]  =  (rLowPoint2[i]  >  geom_2.GetPoint(point)[i] ) ?  geom_2.GetPoint(point)[i] : rLowPoint2[i];
                rHighPoint2[i] =  (rHighPoint2[i] <  geom_2.GetPoint(point)[i] ) ?  geom_2.GetPoint(point)[i] : rHighPoint2[i];
            }
        }

        if (rHighPoint1[0] < rLowPoint2[0])
        {
            return false;
        }
        else if (rLowPoint1[0] > rHighPoint2[0])
        {
            return false;
        }
        else if (rLowPoint1[1] > rHighPoint2[1])
        {
            return false;
        }
        else if (rHighPoint1[1] < rLowPoint2[1])
        {
            return false;
        }

        return true;
    }

    static inline bool Intersection(const PointerType& rObj_1, const PointerType& rObj_2, const double & Radius)
    {
        KRATOS_ERROR << "Fast Quad spatial does not define this function" << std::endl;
    }

    //******************************************************************************************************************

    static inline bool IntersectionBox(const PointerType& rObject,  const PointType& rLowPoint, const PointType& rHighPoint)
    {
        Element::GeometryType& geom_1 = rObject->GetGeometry();

        PointType rHighPoint1 = geom_1.GetPoint(0);
        PointType rLowPoint1  = geom_1.GetPoint(0);

        // std::cout << "POINTINTERBOX" << std::endl;
        // std::cout << "_____________" << std::endl;
        // std::cout << "Point number: " << geom_1.PointsNumber() << std::endl;
        // std::cout << "_____________" << std::endl;
        for (unsigned int point = 1; point<geom_1.PointsNumber(); point++)
        {
            // std::cout << "\t" << geom_1.GetPoint(point)[0] << " " << geom_1.GetPoint(point)[1] << std::endl;
            for(std::size_t i = 0; i<TDimension; i++)
            {
                rLowPoint1[i]  =  (rLowPoint1[i]  >  geom_1.GetPoint(point)[i] ) ?  geom_1.GetPoint(point)[i] : rLowPoint1[i];
                rHighPoint1[i] =  (rHighPoint1[i] <  geom_1.GetPoint(point)[i] ) ?  geom_1.GetPoint(point)[i] : rHighPoint1[i];
            }
        }

        // std::cout << "_____________" << std::endl;
        // std::cout << rLowPoint1[0] << " " << rLowPoint1[1] << std::endl;
        // std::cout << rHighPoint1[0] << " " << rHighPoint1[1] << std::endl;
        // std::cout << rLowPoint[0] << " " << rLowPoint[1] << std::endl;
        // std::cout << rHighPoint[0] << " " << rHighPoint[1] << std::endl;

        if (rHighPoint1[0] < rLowPoint[0])
        {
            return false;
        }
        else if (rLowPoint1[0] > rHighPoint[0])
        {
            return false;
        }
        else if (rLowPoint1[1] > rHighPoint[1])
        {
            return false;
        }
        else if (rHighPoint1[1] < rLowPoint[1])
        {
            return false;
        }

        return true;
    }

    static inline bool IntersectionBox(const PointerType& rObject,  const PointType& rLowPoint, const PointType& rHighPoint, const double & Radius)
    {
        KRATOS_ERROR << "Fast Quad spatial does not define this function" << std::endl;
    }

    static inline void Distance(const PointerType& rObj_1, const PointerType& rObj_2, double& distance)
    {
        KRATOS_ERROR << "Fast Quad spatial does not define this function" << std::endl;
    }

    //static inline void Distance(const PointerType& rObj_1, const PointerType& rObj_2, double& distance)
    //{
    //array_1d<double, 3> center_of_particle1 = rObj_1->GetGeometry()[0];
    //array_1d<double, 3> center_of_particle2 = rObj_2->GetGeometry()[0];

    //distance = sqrt((center_of_particle1[0] - center_of_particle2[0]) * (center_of_particle1[0] - center_of_particle2[0]) +
    //(center_of_particle1[1] - center_of_particle2[1]) * (center_of_particle1[1] - center_of_particle2[1]) +
    //(center_of_particle1[2] - center_of_particle2[2]) * (center_of_particle1[2] - center_of_particle2[2]) );
    //}

    //******************************************************************************************************************

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
        return " Spatial Containers Configure for Quadrilaterals";
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

    //static inline bool floateq(double a, double b) {
    //return std::fabs(a - b) < std::numeric_limits<double>::epsilon();
    //}

    //static inline bool floatle(double a, double b) {
    //return std::fabs(a - b) < std::numeric_limits<double>::epsilon() || a < b;
    //}

    //static inline bool floatge(double a, double b) {
    //return std::fabs(a - b) < std::numeric_limits<double>::epsilon() || a > b;
    //}

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
    FastQuadSpatialContainersConfigure& operator=(FastQuadSpatialContainersConfigure const& rOther);

    /// Copy constructor.
    FastQuadSpatialContainersConfigure(FastQuadSpatialContainersConfigure const& rOther);

    ///@}

}; // Class ParticleConfigure

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template <std::size_t TDimension>
inline std::istream& operator >> (std::istream& rIStream, FastQuadSpatialContainersConfigure<TDimension> & rThis)
{
    return rIStream;
}

/// output stream function
template <std::size_t TDimension>
inline std::ostream& operator << (std::ostream& rOStream, const FastQuadSpatialContainersConfigure<TDimension>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

}   // namespace Kratos.
#endif	/* FAST_QUAD_SPATIAL_CONTAINERS_CONFIGURE_H */
