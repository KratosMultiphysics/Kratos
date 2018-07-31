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
#include "custom_python/add_custom_meshers_to_python.h"

// Meshers
#include "custom_meshers/triangular_mesh_2D_mesher.hpp"
#include "custom_meshers/tetrahedral_mesh_3D_mesher.hpp"

// Bounding Boxes


namespace Kratos
{

namespace Python
{

void  AddCustomMeshersToPython(pybind11::module& m)
{
  using namespace pybind11;

  //class that allos remeshing and adaptive refining (inserting and erasing nodes)
  class_<Mesher, typename Mesher::Pointer>(m,"Mesher")
      .def(init< >())
      .def("Initialize",&Mesher::Initialize)
      .def("InitializeMesher",&Mesher::InitializeMesher)
      .def("FinalizeMesher",&Mesher::FinalizeMesher)
      .def("SetMeshingParameters",&Mesher::SetMeshingParameters)
      .def("SetPreMeshingProcess",&Mesher::SetPreMeshingProcess)
      .def("SetPostMeshingProcess",&Mesher::SetPostMeshingProcess)
      .def("SetPreMeshingProcessVector",&Mesher::SetPreMeshingProcessVector)
      .def("SetPostMeshingProcessVector",&Mesher::SetPostMeshingProcessVector)
      .def("SetMesherUtilities",&Mesher::SetMesherUtilities)
      .def("SetDataTransferUtilities",&Mesher::SetDataTransferUtilities)
      .def("SetEchoLevel",&Mesher::SetEchoLevel)
      .def("ExecuteMeshing",&Mesher::ExecuteMeshing)
      ;

  //class that allows 3D adaptive remeshing (inserting and erasing nodes)
  class_<TetrahedralMesh3DMesher, typename TetrahedralMesh3DMesher::Pointer, Mesher>
      (m,"TetrahedralMesh3DMesher")
      .def(init< >())
      ;

  //class that allows 2D adaptive remeshing (inserting and erasing nodes)
  class_<TriangularMesh2DMesher, typename TriangularMesh2DMesher::Pointer, Mesher>
      (m,"TriangularMesh2DMesher")
      .def(init< >())
      ;

}

}  // namespace Python.

} // Namespace Kratos

