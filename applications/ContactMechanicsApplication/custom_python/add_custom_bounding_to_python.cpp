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

// Project includes
#include "custom_python/add_custom_bounding_to_python.h"

// Bounding Boxes
#include "custom_bounding/spatial_bounding_box.hpp"
#include "custom_bounding/plane_bounding_box.hpp"
#include "custom_bounding/sphere_bounding_box.hpp"
#include "custom_bounding/circle_bounding_box.hpp"
#include "custom_bounding/cylinder_bounding_box.hpp"
#include "custom_bounding/tube_bounding_box.hpp"
#include "custom_bounding/compound_noses_bounding_box.hpp"

namespace Kratos
{

namespace Python
{

void  AddCustomBoundingToPython(pybind11::module& m)
{

  namespace py = pybind11;

  //plane-wall
  py::class_<PlaneBoundingBox, typename PlaneBoundingBox::Pointer, SpatialBoundingBox>(m,"PlaneBoundingBox")
      .def(py::init< Vector, Vector, Vector, int >() )
      .def(py::init< Parameters >())
      .def(py::init< Parameters& >())
      .def("CreateBoundingBoxBoundaryMesh",&PlaneBoundingBox::CreateBoundingBoxBoundaryMesh)
      ;

  //sphere-wall
  py::class_<SphereBoundingBox, typename SphereBoundingBox::Pointer, SpatialBoundingBox>(m,"SphereBoundingBox")
      .def(py::init< Vector, double, Vector, int >() )
      .def(py::init< Parameters >())
      .def(py::init< Parameters& >())
      .def("CreateBoundingBoxBoundaryMesh",&SphereBoundingBox::CreateBoundingBoxBoundaryMesh)
      ;

  //circle-wall
  py::class_<CircleBoundingBox, typename CircleBoundingBox::Pointer, SphereBoundingBox>(m,"CircleBoundingBox")
      .def(py::init< Vector, double, Vector, int >() )
      .def(py::init< Parameters >())
      .def(py::init< Parameters& >())
      .def("CreateBoundingBoxBoundaryMesh",&CircleBoundingBox::CreateBoundingBoxBoundaryMesh)
      ;

  //cylinder-wall
  py::class_<CylinderBoundingBox, typename CylinderBoundingBox::Pointer, SpatialBoundingBox>(m,"CylinderBoundingBox")
      .def(py::init< Vector, Vector, double, Vector, int >())
      .def(py::init< Parameters >())
      .def(py::init< Parameters& >())
      .def("CreateBoundingBoxBoundaryMesh",&CylinderBoundingBox::CreateBoundingBoxBoundaryMesh)
      ;

  //tube-wall
  py::class_<TubeBoundingBox, typename TubeBoundingBox::Pointer, SpatialBoundingBox>(m,"TubeBoundingBox")
      .def(py::init< ModelPart&, double, int >())
      .def(py::init< ModelPart&, Parameters >())
      .def(py::init< ModelPart&, Parameters& >())
      .def("CreateBoundingBoxBoundaryMesh",&TubeBoundingBox::CreateBoundingBoxBoundaryMesh)
      ;

  //compound_noses-wall
  py::class_<CompoundNosesBoundingBox, typename CompoundNosesBoundingBox::Pointer, SpatialBoundingBox>(m,"CompoundNosesBoundingBox")
      .def(py::init< Vector, Vector, Vector, Matrix, Vector, Vector, Vector, Vector, Vector, Matrix >() )
      .def(py::init< Parameters >())
      .def(py::init< Parameters& >())
      .def("CreateBoundingBoxBoundaryMesh",&CompoundNosesBoundingBox::CreateBoundingBoxBoundaryMesh)
      ;


}

}  // namespace Python.

} // Namespace Kratos

