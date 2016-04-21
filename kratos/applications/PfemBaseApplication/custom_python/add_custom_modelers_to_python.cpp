//
//   Project Name:        KratosPfemBaseApplication $
//   Created by:          $Author:      JMCarbonell $
//   Last modified by:    $Co-Author:               $
//   Date:                $Date:      February 2016 $
//   Revision:            $Revision:            0.0 $
//
//


// System includes

// External includes
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_modelers_to_python.h"

// Meshers
#include "custom_modelers/triangular_mesh_2D_modeler.hpp"
#include "custom_modelers/tetrahedral_mesh_3D_modeler.hpp"

// Bounding Boxes


namespace Kratos
{

namespace Python
{

  typedef MeshModeler                        MeshModelerBaseType;
  typedef MeshModeler::Pointer                MeshModelerPointer;


  void  AddCustomModelersToPython()
  {
    using namespace boost::python;
  
    //class that allos remeshing and adaptive refining (inserting and erasing nodes)
    class_<MeshModeler, MeshModeler::Pointer, boost::noncopyable >
      ("MeshModeler", init< >())
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
      .def("GenerateMesh",&MeshModeler::GenerateMesh)
      .def("SetEchoLevel",&MeshModeler::SetEchoLevel)
      ;

    //class that allows 3D adaptive remeshing (inserting and erasing nodes)
    class_<TetrahedralMesh3DModeler, bases<MeshModelerBaseType>, boost::noncopyable >
      ("TetrahedralMesh3DModeler",init< >())
      ;
    
    //class that allows 2D adaptive remeshing (inserting and erasing nodes)
    class_<TriangularMesh2DModeler, bases<MeshModelerBaseType>, boost::noncopyable >
      ("TriangularMesh2DModeler",init< >())
      ;



  }

}  // namespace Python.

} // Namespace Kratos

