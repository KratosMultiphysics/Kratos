//
//   Project Name:        KratosPfemApplication     $
//   Created by:          $Author:      JMCarbonell $
//   Last modified by:    $Co-Author:               $
//   Date:                $Date:      February 2016 $
//   Revision:            $Revision:            0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_python/add_custom_modelers_to_python.h"

// Meshers
#include "custom_modelers/triangular_mesh_2D_modeler.hpp"
#include "custom_modelers/tetrahedral_mesh_3D_modeler.hpp"

// Bounding Boxes


namespace Kratos
{

namespace Python
{

void  AddCustomModelersToPython(pybind11::module& m)
{
  using namespace pybind11;
  
  //class that allows remeshing and adaptive refining (inserting and erasing nodes)
  class_<MeshModeler, typename MeshModeler::Pointer>(m,"MeshModeler")
      .def(init< >())
      .def("Initialize",&MeshModeler::Initialize)
      .def("InitializeMeshModeler",&MeshModeler::InitializeMeshModeler)
      .def("FinalizeMeshModeler",&MeshModeler::FinalizeMeshModeler)
      .def("SetMeshingParameters",&MeshModeler::SetMeshingParameters)
      .def("SetPreMeshingProcess",&MeshModeler::SetPreMeshingProcess)
      .def("SetPostMeshingProcess",&MeshModeler::SetPostMeshingProcess)
      .def("SetPreMeshingProcessVector",&MeshModeler::SetPreMeshingProcessVector)
      .def("SetPostMeshingProcessVector",&MeshModeler::SetPostMeshingProcessVector)
      .def("SetModelerUtilities",&MeshModeler::SetModelerUtilities)
      .def("SetDataTransferUtilities",&MeshModeler::SetDataTransferUtilities)
      .def("SetEchoLevel",&MeshModeler::SetEchoLevel)
      .def("ExecuteMeshing",&MeshModeler::ExecuteMeshing)
      ;

  //class that allows 3D adaptive remeshing (inserting and erasing nodes)
  class_<TetrahedralMesh3DModeler, typename TetrahedralMesh3DModeler::Pointer, MeshModeler>
      (m,"TetrahedralMesh3DModeler")
      .def(init< >())
      ;
    
  //class that allows 2D adaptive remeshing (inserting and erasing nodes)
  class_<TriangularMesh2DModeler, typename TriangularMesh2DModeler::Pointer, MeshModeler>
      (m,"TriangularMesh2DModeler")
      .def(init< >())
      ;

}

}  // namespace Python.

} // Namespace Kratos

