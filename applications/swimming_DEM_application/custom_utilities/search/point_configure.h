//   Project Name:        Kratos
//   Last Modified by:    $Author: Nelson Lafontaine  $
//   Date:                $Date: 2012-05-24 $
//   Revision:            $Revision: 1.0 $
//


#if !defined(KRATOS_POINT_CONFIGURE)
#define  KRATOS_POINT_CONFIGURE



// System includes
#include <string>
#include <iostream>
#include <cmath>
#include "utilities/spatial_containers_configure.h"

// Kratos includes
#include "includes/variables.h"
#include "spatial_containers/spatial_search.h"
//#include "GeometryFunctions.h"

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
class PointConfigure{
public:

enum {Dimension = TDimension,
      DIMENSION = TDimension,
      MAX_LEVEL = 16,
      MIN_LEVEL = 2
};

typedef SpatialSearch                                                       SearchType;

typedef SearchType::PointType                                               PointType;
typedef PointerVectorSet<Point, IndexedObject>::ContainerType   ContainerType;
typedef PointerVectorSet<Point, IndexedObject>                  PointsContainerType;

typedef SearchType::ElementType                                             ElementType;
typedef ContainerType::value_type                                           PointerType;
typedef ContainerType::iterator                                             IteratorType;

typedef PointerVectorSet<Point, IndexedObject>::ContainerType   ResultContainerType;


typedef ResultContainerType::iterator                           ResultIteratorType;
typedef std::vector<double>::iterator                           DistanceIteratorType;

typedef ContactPair<PointerType>                                ContactPairType;
typedef std::vector<ContactPairType>                            ContainerContactType;
typedef ContainerContactType::iterator                          IteratorContactType;
typedef ContainerContactType::value_type                        PointerContactType;



/// Pointer definition of SpatialContainersConfigure
KRATOS_CLASS_POINTER_DEFINITION(PointConfigure);

///@}
///@name Life Cycle
///@{

PointConfigure(){};
virtual ~PointConfigure(){}

  ///@}
  ///@name Operators
  ///@{


  ///@}
  ///@name Operations
  ///@{


//******************************************************************************************************************
//******************************************************************************************************************

static inline void CalculateBoundingBox(const PointerType& r_p_point, PointType& r_low_point, PointType& r_high_point){
    r_high_point = r_low_point = *r_p_point;
}

//******************************************************************************************************************
//******************************************************************************************************************

static inline void CalculateBoundingBox(const PointerType& r_p_point, PointType& r_low_point, PointType& r_high_point, const double& radius){
    r_high_point = r_low_point  = *r_p_point;

    for (std::size_t i = 0; i < 3; ++i){
        r_low_point[i]  -= radius;
        r_high_point[i] += radius;
    }
}

//******************************************************************************************************************
//******************************************************************************************************************

static inline void CalculateCenter(const PointerType& r_p_point, PointType& rCenter){
    rCenter = *r_p_point;
}

//******************************************************************************************************************
//******************************************************************************************************************

static inline bool Intersection(const PointerType& r_p_point_1, const PointerType& r_p_point_2){
    return false;
}

//******************************************************************************************************************
//******************************************************************************************************************

static inline bool Intersection(const PointerType& r_p_point_1, const PointerType& r_p_point_2, const double& radius){
    array_1d<double, 3> node_2_to_1 = *r_p_point_1 - *r_p_point_2;
    double distance;
    Distance(r_p_point_1, r_p_point_2, distance);

    bool intersect = (distance - radius) <= 0;

    return intersect;
}

//******************************************************************************************************************
//******************************************************************************************************************

static inline bool  IntersectionBox(const PointerType& r_p_point,  const PointType& r_low_point, const PointType& r_high_point)
{
    array_1d<double, 3> center = *r_p_point;

    bool intersect = (r_low_point[0] <= center[0] &&  r_low_point[1] <= center[1] &&  r_low_point[2] <= center[2] &&
                     r_high_point[0] >= center[0] && r_high_point[1] >= center[1] && r_high_point[2] >= center[2]);

    return  intersect;

}

//******************************************************************************************************************
//******************************************************************************************************************

static inline bool  IntersectionBox(const PointerType& r_p_point,  const PointType& r_low_point, const PointType& r_high_point, const double& radius)
{
    array_1d<double, 3> center = *r_p_point;

    bool intersect = (r_low_point[0] - radius <= center[0] &&  r_low_point[1] - radius <= center[1] &&  r_low_point[2] - radius <= center[2] &&
                     r_high_point[0] + radius >= center[0] && r_high_point[1] + radius >= center[1] && r_high_point[2] + radius >= center[2]);

    return  intersect;
}

//******************************************************************************************************************
//******************************************************************************************************************

static inline void Distance(const PointerType& r_p_point_1, const PointerType& r_p_point_2, double& distance)
{
    array_1d<double, 3> center_1 = *r_p_point_1;
    array_1d<double, 3> center_2 = *r_p_point_2;

    distance = sqrt((center_1[0] - center_2[0]) * (center_1[0] - center_2[0]) +
                    (center_1[1] - center_2[1]) * (center_1[1] - center_2[1]) +
                    (center_1[2] - center_2[2]) * (center_1[2] - center_2[2]));
}

//******************************************************************************************************************
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
virtual std::string Info() const {return " Spatial Containers Configure for Particles"; }

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
PointConfigure& operator=(PointConfigure const& rOther);

/// Copy constructor.
PointConfigure(PointConfigure const& rOther);

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
    inline std::istream& operator >> (std::istream& rIStream, PointConfigure<TDimension> & rThis){
        return rIStream;
        }

    /// output stream function
    template <std::size_t TDimension>
    inline std::ostream& operator << (std::ostream& rOStream, const PointConfigure<TDimension>& rThis){
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
        }

    ///@}

}   // namespace Kratos.
#endif	/* KRATOS_POINT_CONFIGURE */
