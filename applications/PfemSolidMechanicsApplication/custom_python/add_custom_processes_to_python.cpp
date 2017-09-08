//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2013 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes 
#include <boost/python.hpp>

// External includes 

// Project includes
#include "includes/node.h"
#include "includes/define.h"
#include "processes/process.h"

//Application includes
#include "custom_python/add_custom_processes_to_python.h"

//Processes
#include "custom_processes/contact_refine_mesh_boundary_process.hpp"
#include "custom_processes/set_mechanical_initial_state_process.hpp"

namespace Kratos
{
	
  namespace Python
  {
 	
    void  AddCustomProcessesToPython()
    {

      using namespace boost::python;
      typedef Process                                         ProcessBaseType;
      typedef RefineMeshBoundaryProcess             RefineMeshProcessBaseType;
      typedef std::vector<SpatialBoundingBox::Pointer>   BoundingBoxContainer;

      
      //**********MESH MODELLER PROCESS*********//

      class_<ContactRefineMeshBoundaryProcess, bases<RefineMeshProcessBaseType>, boost::noncopyable >
	(
	 "ContactRefineMeshBoundary", init<ModelPart&, BoundingBoxContainer&, ModelerUtilities::MeshingParameters&, int>()
	 )
	;

      
         // **** SET INITIAL MECHANICAL STATE **** //
         class_<SetMechanicalInitialStateProcess, bases<ProcessBaseType>, boost::noncopyable >
            (
             "SetMechanicalInitialStateProcess", init<ModelPart&, Parameters>()
            )
            .def(init< ModelPart&, Parameters >())
            .def("Execute", &SetMechanicalInitialStateProcess::Execute)
         ;


 

    }
 
  }  // namespace Python.

} // Namespace Kratos

