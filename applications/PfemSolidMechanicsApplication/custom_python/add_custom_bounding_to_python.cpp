//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                February 2016 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes

// External includes
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_bounding_to_python.h"

// Meshers

// Bounding Boxes
#include "custom_bounding/spatial_bounding_box.hpp"
#include "custom_bounding/rigid_nose_wall_bounding_box.hpp"
#include "custom_bounding/rigid_circle_wall_bounding_box.hpp"
#include "custom_bounding/rigid_plane_wall_bounding_box.hpp"

namespace Kratos
{

namespace Python
{

  typedef SpatialBoundingBox                        BoundingBoxBaseType;
  typedef SpatialBoundingBox::Pointer                BoundingBoxPointer;
  typedef std::vector<SpatialBoundingBox::Pointer> BoundingBoxContainer;


  void  AddCustomBoundingToPython()
  {

    using namespace boost::python;

    //nose-wall
    class_<RigidNoseWallBoundingBox, bases<BoundingBoxBaseType>, boost::noncopyable > 
      ( "RigidNoseWallBoundingBox", 
	init<int, Vector, Vector, Vector, Vector, Matrix, Vector, Vector, Vector>() )
      .def("SetAxisymmetric",&RigidNoseWallBoundingBox::SetAxisymmetric)
      .def("SetDimension",&RigidNoseWallBoundingBox::SetDimension)
      ;
    
    //circle-wall
    class_<RigidCircleWallBoundingBox, bases<BoundingBoxBaseType>, boost::noncopyable > 
      ( "RigidCircleWallBoundingBox", 
	init<int, int, double, Vector, Vector, Vector, Vector>() )
      .def("SetAxisymmetric",&RigidCircleWallBoundingBox::SetAxisymmetric)
      .def("SetDimension",&RigidCircleWallBoundingBox::SetDimension)
      ;

    //plane-wall
    class_<RigidPlaneWallBoundingBox, bases<BoundingBoxBaseType>, boost::noncopyable > 
      ( "RigidPlaneWallBoundingBox", 
	init<int, int, Vector, Vector, Vector, Vector, Vector>() )
      .def("SetAxisymmetric",&RigidPlaneWallBoundingBox::SetAxisymmetric)
      .def("SetDimension",&RigidPlaneWallBoundingBox::SetDimension)
      ;


     
  }

}  // namespace Python.

} // Namespace Kratos

