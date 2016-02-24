//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2013 $
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
#include "custom_python/add_custom_modelers_to_python.h"

// Meshers
#include "custom_modelers/contact_domain_2D_modeler.hpp"
//#include "custom_modelers/contact_domain_3D_modeler.hpp"

// Bounding Boxes


namespace Kratos
{

namespace Python
{

  typedef MeshModeler                        MeshModelerBaseType;
  typedef MeshModeler::Pointer                MeshModelerPointer;


  void GenerateTriangleContactMesh(ContactDomain2DModeler& Mesher, ModelPart& model_part, double penalty_parameter, double stability_parameter, bool friction_active, double nu_static, double nu_dynamic)
  {
    Mesher.GenerateContactMesh(model_part,
			       penalty_parameter, stability_parameter,
			       friction_active, nu_static, nu_dynamic);
  }
    
  void SetMeshingParameters(ContactDomain2DModeler& Mesher, ModelerUtilities::MeshingParameters::Pointer MeshingParameters, int MeshId)
  {
    Mesher.SetMeshingParameters(MeshingParameters, MeshId);
  }

  void  AddCustomModelersToPython()
  {

    using namespace boost::python;
    //class that allows 3D adaptive remeshing (inserting and erasing nodes)

    
    //class that allows 2D adaptive remeshing (inserting and erasing nodes)


    //class that allows 2D contact domain spatial search
    class_<ContactDomain2DModeler, bases<MeshModelerBaseType>, boost::noncopyable >
      ("ContactDomain2DModeler", init< >())
      .def("GenerateContactMesh",GenerateTriangleContactMesh)
      .def("SetMechingParamters",SetMeshingParameters)
      .def("TransferContactData",&ContactDomain2DModeler::TransferContactBoundaryData)
      .def("SetEchoLevel",&ContactDomain2DModeler::SetEchoLevel)
      ;
     
  }

}  // namespace Python.

} // Namespace Kratos

