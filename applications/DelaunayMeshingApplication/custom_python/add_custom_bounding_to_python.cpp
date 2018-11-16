//
//   Project Name:        KratosDelaunayMeshingApplication $
//   Created by:          $Author:             JMCarbonell $
//   Last modified by:    $Co-Author:                      $
//   Date:                $Date:                April 2018 $
//   Revision:            $Revision:                   0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_python/add_custom_bounding_to_python.h"

// Bounding Boxes
#include "custom_bounding/spatial_bounding_box.hpp"

namespace Kratos
{

namespace Python
{
typedef SpatialBoundingBox::Pointer                BoundingBoxPointer;
typedef std::vector<SpatialBoundingBox::Pointer> BoundingBoxContainer;

void Push_Back_Bounding_Box( BoundingBoxContainer& ThisBoundingBoxContainer,
                             BoundingBoxPointer ThisBoundingBox )
{
  ThisBoundingBoxContainer.push_back( ThisBoundingBox );
}


void  AddCustomBoundingToPython(pybind11::module& m)
{

  namespace py = pybind11;

  //bounding box container
  py::class_<BoundingBoxContainer>(m, "BoundingBoxContainer")
      .def(py::init<>())
      .def("PushBack", Push_Back_Bounding_Box)
      ;

  //spatial bounding box
  py::class_<SpatialBoundingBox, typename SpatialBoundingBox::Pointer>
      (m, "SpatialBoundingBox")
      .def(py::init<Vector, Vector>())
      .def(py::init< Parameters >())
      .def(py::init< Parameters& >())
      .def("SetAxisymmetric",&SpatialBoundingBox::SetAxisymmetric)
      .def("SetDimension",&SpatialBoundingBox::SetDimension)
      .def("SetUpperPoint",&SpatialBoundingBox::SetUpperPoint)
      .def("SetLowerPoint",&SpatialBoundingBox::SetLowerPoint)
      .def("SetRigidBodyCenter",&SpatialBoundingBox::SetRigidBodyCenter)
      .def("CreateBoundingBoxBoundaryMesh",&SpatialBoundingBox::CreateBoundingBoxBoundaryMesh)
      ;

}

}  // namespace Python.

} // Namespace Kratos
