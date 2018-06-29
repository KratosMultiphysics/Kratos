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

// Project includes
#include "custom_python/add_custom_processes_to_python.h"

//Processes
#include "custom_processes/contact_refine_mesh_boundary_process.hpp"
#include "custom_processes/set_mechanical_initial_state_process.hpp"

namespace Kratos
{

namespace Python
{

void  AddCustomProcessesToPython(pybind11::module& m)
{

  using namespace pybind11;

  typedef std::vector<SpatialBoundingBox::Pointer>   BoundingBoxContainer;

  // Mesh modeler process
  class_<ContactRefineMeshBoundaryProcess, ContactRefineMeshBoundaryProcess::Pointer, RefineMeshBoundaryProcess>(m,"ContactRefineMeshBoundary")
      .def(init<ModelPart&, BoundingBoxContainer&, ModelerUtilities::MeshingParameters&, int>())
      ;

  // Set initial mechanical state
  class_<SetMechanicalInitialStateProcess, SetMechanicalInitialStateProcess::Pointer, Process>(m,"SetMechanicalInitialStateProcess")
      .def(init<ModelPart&, Parameters>())
      .def(init< ModelPart&, Parameters >())
      .def("Execute", &SetMechanicalInitialStateProcess::Execute)
      ;
}

}  // namespace Python.

} // Namespace Kratos
