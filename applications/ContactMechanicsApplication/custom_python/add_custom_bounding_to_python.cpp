//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2016 $
//   Revision:            $Revision:                    0.0 $
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
#include "includes/kratos_parameters.h"

// Bounding Boxes
#include "custom_bounding/spatial_bounding_box.hpp"
#include "custom_bounding/plane_bounding_box.hpp"
#include "custom_bounding/sphere_bounding_box.hpp"
#include "custom_bounding/circle_bounding_box.hpp"
#include "custom_bounding/cylinder_bounding_box.hpp"
#include "custom_bounding/compound_noses_bounding_box.hpp"

namespace Kratos
{

namespace Python
{

  typedef SpatialBoundingBox                        BoundingBoxBaseType;
  typedef SpatialBoundingBox::Pointer                BoundingBoxPointer;
  typedef std::vector<SpatialBoundingBox::Pointer> BoundingBoxContainer;

  typedef SphereBoundingBox                   SphereBoundingBoxBaseType;


  void  AddCustomBoundingToPython()
  {

    using namespace boost::python;

    //plane-wall
    class_<PlaneBoundingBox, bases<BoundingBoxBaseType>, boost::noncopyable > 
      ( "PlaneBoundingBox", 
        init< Vector, Vector, Vector, int >() )
        .def(init< Parameters >())
        .def(init< Parameters& >())
      ;

    //sphere-wall
    class_<SphereBoundingBox, bases<BoundingBoxBaseType>, boost::noncopyable > 
      ( "SphereBoundingBox", 
	init< Vector, double, Vector, int >() )
        .def(init< Parameters >())
        .def(init< Parameters& >())
      ;

    //circle-wall
    class_<CircleBoundingBox, bases<SphereBoundingBoxBaseType>, boost::noncopyable > 
      ( "CircleBoundingBox", 
	init< Vector, double, Vector, int >() )
        .def(init< Parameters >())
        .def(init< Parameters& >())
        .def("CreateBoundingBoxBoundaryMesh",&CircleBoundingBox::CreateBoundingBoxBoundaryMesh)
      ;

    //cylinder-wall
    class_<CylinderBoundingBox, bases<BoundingBoxBaseType>, boost::noncopyable > 
      ( "CylinderBoundingBox", 
	init< Vector, Vector, double, Vector, int >() )
        .def(init< Parameters >())
        .def(init< Parameters& >())
        .def("CreateBoundingBoxBoundaryMesh",&CylinderBoundingBox::CreateBoundingBoxBoundaryMesh)
      ;
    
    //compound_noses-wall
    class_<CompoundNosesBoundingBox, bases<BoundingBoxBaseType>, boost::noncopyable > 
      ( "CompoundNosesBoundingBox", 
	init< Vector, Vector, Vector, Matrix, Vector, Vector, Vector, Vector, Vector, Matrix >() )
        .def(init< Parameters >())
        .def(init< Parameters& >())
      ;
    

  }

}  // namespace Python.

} // Namespace Kratos

