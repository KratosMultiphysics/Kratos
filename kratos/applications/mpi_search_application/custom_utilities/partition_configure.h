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
#include <string>
#include <iostream>
#include <limits>
#include <cmath>

// Kratos includes
#include "includes/define.h"
#include "geometries/point.h"

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

/** Configuration file for Partitions.
 * This class provides a configuration file to calculate a 'Bins'
 * for the partitions of a model. For this proupose, the objects
 * are Points in the space.
 */
class PartitionConfigure {
public:

  /// Pointer definition of PartitionConfigure
  KRATOS_CLASS_POINTER_DEFINITION(PartitionConfigure);

  /** Compile time definitions
   * @param EPSILON   Error tolerance for cmparison operations with doubles
   * @param DIMENSION Dimension of the problem. Fixed to 3.
   */
  static constexpr auto EPSILON   = std::numeric_limits<double>::epsilon();
  static constexpr auto DIMENSION = 3;

  /** Point and Pointer Types
   * @param PointType   Point of doubles with 3 coordinates (DIMENSION = 3)
   * @param PointerType Pointer to Point of doubles with 3 coordinates (DIMENSION = 3)
   */
  typedef Point<DIMENSION, double>          PointType;
  typedef Point<DIMENSION, double>::Pointer PointerType;

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
   * For this configuation file, the bounding box is the equal to the point given in 'rObject'.
   * @param rObject    Point for which the bounding box will be calculated.
   * @param rLowPoint  Lower point of the boundingbox.
   * @param rHighPoint Higher point of the boundingbox.
   */
  static inline void CalculateBoundingBox(const PointerType& rObject, PointType& rLowPoint, PointType& rHighPoint) {
    rHighPoint = rLowPoint = *rObject;
  }

  /** Calculates the bounding box for the given object extended with a Radius.
   * For this configuation file, the bounding box is the equal to the point given in 'rObject' + - a radius.
   * @param rObject    Point for which the bounding box will be calculated.
   * @param rLowPoint  Lower point of the boundingbox.
   * @param rHighPoint Higher point of the boundingbox.
   * @param Radius     The extension radius to be applied to the boundingbox.
   */
  static inline void CalculateBoundingBox(const PointerType& rObject, PointType& rLowPoint, PointType& rHighPoint, const double& Radius) {
    auto radiusExtension = PointType(Radius, Radius, Radius);

    rLowPoint = *rObject - radiusExtension;
    rLowPoint = *rObject + radiusExtension;
  }

  /** Tests the intersection of two objects
   * For this configuation file, tests if the two points are the same within a EPSILON tolerance range.
   * @param  rObj_1 First point of the tests
   * @param  rObj_2 Second point of the tests
   * @return        Boolean indicating the result of the intersection test described.
   */
  static inline bool Intersection(const PointerType& rObj_1, const PointerType& rObj_2) {
    for(std::size_t i = 0; i < DIMENSION; i++) {
      if(std::fabs((*rObj_1)[i] - (*rObj_2)[i]) > EPSILON) {
        return false;
      }
    }

    return true;
  }

  /** Tests the intersection of two objects extended with a given radius.
   * For this configuation file, tests if the two points extended with a radius
   * are the same within a EPSILON tolerance range.
   * @param  rObj_1 First point of the tests
   * @param  rObj_2 Second point of the tests
   * @param  Radius The extension radius to be applied in the intersection.
   * @return        Boolean indicating the result of the intersection test described.
   */
  static inline bool Intersection(const PointerType& rObj_1, const PointerType& rObj_2, double Radius) {
    for(std::size_t i = 0; i < DIMENSION; i++) {
      if(std::fabs((*rObj_1)[i] - (*rObj_2)[i]) > EPSILON + 2.0f * Radius) {
        return false;
      }
    }

    return true;
  }

  /** Tests the intersection of one object with a boundingbox descrived by 'rLowPoint' and 'rHighPoint'.
   * For this configuation file, tests if one point is inside the boundingbox
   * described by 'rLowPoint' and 'rHighPoint' within a EPSILON tolerance range.
   * @param  rObject    Point of the tests.
   * @param  rLowPoint  Lower point of the boundingbox.
   * @param  rHighPoint Higher point of the boundingbox.
   * @return            Boolean indicating the result of the intersection test described.
   */
  static inline bool IntersectionBox(const PointerType& rObject, const PointType& rLowPoint, const PointType& rHighPoint) {
    for(std::size_t i = 0; i < DIMENSION; i++) {
      if( (*rObject)[i] < rLowPoint[i] - EPSILON || (*rObject)[i] > rHighPoint[i] + EPSILON) {
        return false;
      }
    }

    return true;
  }

  /** Tests the intersection of one object with a boundingbox descrived by 'rLowPoint' and 'rHighPoint'.
   * For this configuation file, tests if one point extended by radius is inside the boundingbox
   * described by 'rLowPoint' and 'rHighPoint' within a EPSILON tolerance range.
   * @param  rObject    Point of the tests.
   * @param  rLowPoint  Lower point of the boundingbox.
   * @param  rHighPoint Higher point of the boundingbox.
   * @param  Radius     The extension radius to be applied in the intersection.
   * @return            Boolean indicating the result of the intersection test described.
   */
  static inline bool IntersectionBox(const PointerType& rObject, const PointType& rLowPoint, const PointType& rHighPoint, const double& Radius) {
    for(std::size_t i = 0; i < DIMENSION; i++) {
      if( ((*rObject)[i] + Radius) < rLowPoint[i] - EPSILON || ((*rObject)[i] - Radius) > rHighPoint[i] + EPSILON) {
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
  static inline void Distance(const PointerType& rObj_1, const PointerType& rObj_2, double& distance) {
    auto pwdDistance = 0;

    for(std::size_t i = 0; i < DIMENSION; i++) {
      pwdDistance += std::pow((*rObj_2)[i] - (*rObj_2)[i], 2);
    }

    distance = std::sqrt(pwdDistance);
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
    return "Spatial Containers Configure for Partitions";
  }

  /// Turns back data as a string.
  virtual std::string Data() const {
    return "DIMENSION: " + std::to_string(DIMENSION);
  }

  /// Prints object's information.
  virtual void PrintInfo(std::ostream& rOStream) const {
    rOStream << Info() << std::endl;
  }

  /// Prints object's data.
  virtual void PrintData(std::ostream& rOStream) const {
    rOStream << Data() << DIMENSION << std::endl;
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

}; // Class ParticleConfigure

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream, PartitionConfigure& rThis){
  return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const PartitionConfigure& rThis){
  rThis.PrintInfo(rOStream);
  rOStream << std::endl;
  rThis.PrintData(rOStream);

  return rOStream;
}

///@}

} // namespace Kratos.
#endif /* MPI_DISCRETE_PARTICLE_CONFIGURE_H */
