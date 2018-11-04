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

#if !defined(KRATOS_PARTITION_CONFIGURE_INCLUDED)
#define  KRATOS_PARTITION_CONFIGURE_INCLUDED

// System includes
#include <unordered_set>
#include <string>
#include <iostream>
#include <limits>
#include <cmath>

// Kratos includes
#include "includes/define.h"
#include "geometries/point.h"
#include "containers/pointer_vector.h"
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

/** Configuration file for Points.
 * This class provides a configuration file to calculate a 'Bins'
 * using points.
 */
template<class BaseConfigure>
class PartitionConfigure {
public:

  /// Pointer definition of PartitionConfigure
  KRATOS_CLASS_POINTER_DEFINITION(PartitionConfigure);

  /** Compile time definitions
   * @param Dimension Dimension of the problem. Fixed to 3.
   */
  static constexpr auto Dimension = 3;

  /** PartitionObject to store data inside the bins points.
   * This is a current workaround the be able to store data inside the bins without
   * loosing its geometrical information.
   * [] Operator access the point info as if it were a CoordinateArray or ObjectType.
   * () Operator access the extension.
   */
  template<class T>
  class PartitionObject: public BaseConfigure::ObjectType {
    public:
      T & operator()() { return container;};
    private:
      T container;
   };

  /** Point and Pointer Types
   * @param PointType   Set of integers
   * @param PointerType Pointer to Set of integers
   */
  typedef typename BaseConfigure::PointType             PointType;

  typedef PartitionObject<std::unordered_set<int>>      ObjectType;
  typedef Kratos::shared_ptr<ObjectType>                PointerType;
  typedef typename BaseConfigure::PointerType           BasePointerType;

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
  typedef PointerVector<PointType>                    PointContainerType;

  typedef typename PointContainerType::ContainerType  ContainerType;
  typedef typename PointContainerType::ContainerType  ResultContainerType;
  typedef ContactPair<PointerType>                    ContactPairType;

  typedef typename ContainerType::iterator            IteratorType;
  typedef typename ResultContainerType::iterator      ResultIteratorType;
  typedef typename std::vector<double>::iterator      DistanceIteratorType;

  typedef typename std::vector<ContactPairType>       ContainerContactType;
  typedef typename ContainerContactType::iterator     IteratorContactType;

  ///@}
  ///@name Life Cycle
  ///@{

  /// Default consturctor
  PartitionConfigure(){};

  /// Default destructor
  virtual ~PartitionConfigure(){}

  ///@}
  ///@name Operators
  ///@{

  ///@}
  ///@name Operations
  ///@{

  /** Calculates the bounding box for the given object.
   * @param rObject    Point for which the bounding box will be calculated.
   * @param rLowPoint  Lower point of the boundingbox.
   * @param rHighPoint Higher point of the boundingbox.
   */
  static inline void CalculateBoundingBox(const BasePointerType& rObject, PointType& rLowPoint, PointType& rHighPoint) {
    BaseConfigure::CalculateBoundingBox(rObject, rLowPoint, rHighPoint);
  }

  /** Calculates the bounding box for the given object extended with a Radius.
   * @param rObject    Point for which the bounding box will be calculated.
   * @param rLowPoint  Lower point of the boundingbox.
   * @param rHighPoint Higher point of the boundingbox.
   * @param Radius     The extension radius to be applied to the boundingbox.
   */
  static inline void CalculateBoundingBox(const BasePointerType& rObject, PointType& rLowPoint, PointType& rHighPoint, const double& Radius) {
    BaseConfigure::CalculateBoundingBox(rObject, rLowPoint, rHighPoint, Radius);
  }

  /** Calculates the Center of the object.
   * @param rObject        Point for which the bounding box will be calculated.
   * @param rCentralPoint  The center point of the object.
   */
  static inline void CalculateCenter(const BasePointerType& rObject, PointType& rCentralPoint) {
    BaseConfigure::CalculateCenter(rObject, rCentralPoint);
  }

  /** Tests the intersection of two objects
   * @param  rObj_1 First point of the tests
   * @param  rObj_2 Second point of the tests
   * @return        Boolean indicating the result of the intersection test described.
   */
  static inline bool Intersection(const BasePointerType& rObj_1, const BasePointerType& rObj_2) {
    return BaseConfigure::Intersection(rObj_1, rObj_2);
  }

  /** Tests the intersection of two objects extended with a given radius.
   * @param  rObj_1 First point of the tests
   * @param  rObj_2 Second point of the tests
   * @param  Radius The extension radius to be applied in the intersection.
   * @return        Boolean indicating the result of the intersection test described.
   */
  static inline bool Intersection(const BasePointerType& rObj_1, const BasePointerType& rObj_2, double Radius) {
    return BaseConfigure::Intersection(rObj_1, rObj_2, Radius);
  }

  /** Tests the intersection of one object with a boundingbox descrived by 'rLowPoint' and 'rHighPoint'.
   * @param  rObject    Point of the tests.
   * @param  rLowPoint  Lower point of the boundingbox.
   * @param  rHighPoint Higher point of the boundingbox.
   * @return            Boolean indicating the result of the intersection test described.
   */
  static inline bool IntersectionBox(const BasePointerType& rObject, const PointType& rLowPoint, const PointType& rHighPoint) {
    return BaseConfigure::IntersectionBox(rObject, rLowPoint, rHighPoint);
  }

  /** Tests the intersection of one object with a boundingbox descrived by 'rLowPoint' and 'rHighPoint'.
   * @param  rObject    Point of the tests.
   * @param  rLowPoint  Lower point of the boundingbox.
   * @param  rHighPoint Higher point of the boundingbox.
   * @param  Radius     The extension radius to be applied in the intersection.
   * @return            Boolean indicating the result of the intersection test described.
   */
  static inline bool IntersectionBox(const BasePointerType& rObject, const PointType& rLowPoint, const PointType& rHighPoint, const double& Radius) {
    return BaseConfigure::IntersectionBox(rObject, rLowPoint, rHighPoint, Radius);
  }

  /** Calculates the distance betwen two objects.
   * @param rObj_1      First point.
   * @param rLowPoint   Lower point.
   * @param rHighPoint  Higher point of the boundingbox.
   * @param distance    The euclidean distance between 'rObj_1' and 'rObj_2'.
   */
  static inline void Distance(const BasePointerType& rObj_1, const BasePointerType& rObj_2, double& distance) {
    BaseConfigure::Distance(rObj_1, rObj_2, distance);
  }

  static std::vector<double> mMinPoint;
  static std::vector<double> mMaxPoint;

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
  PartitionConfigure& operator=(PartitionConfigure const& rOther);

  /// Copy constructor.
  PartitionConfigure(PartitionConfigure const& rOther);

  ///@}

}; // Class PartitionConfigure

template<class BaseConfigure>
std::vector<double> PartitionConfigure<BaseConfigure>::mMinPoint = [] {
  std::vector<double> v;
  for(std::size_t i = 0; i < PartitionConfigure<BaseConfigure>::Dimension; i++) {
    v.push_back(0);
  }
  return v;
}();

template<class BaseConfigure>
std::vector<double> PartitionConfigure<BaseConfigure>::mMaxPoint = [] {
  std::vector<double> v;
  for(std::size_t i = 0; i < PartitionConfigure<BaseConfigure>::Dimension; i++) {
    v.push_back(0);
  }
  return v;
}();

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template<class BaseConfigure>
inline std::istream& operator >> (std::istream& rIStream, PartitionConfigure<BaseConfigure>& rThis){
  return rIStream;
}

/// output stream function
template<class BaseConfigure>
inline std::ostream& operator << (std::ostream& rOStream, const PartitionConfigure<BaseConfigure>& rThis){
  rThis.PrintInfo(rOStream);
  rOStream << std::endl;
  rThis.PrintData(rOStream);

  return rOStream;
}

///@}

} // namespace Kratos.
#endif /* POINT_CONFIGURE */
