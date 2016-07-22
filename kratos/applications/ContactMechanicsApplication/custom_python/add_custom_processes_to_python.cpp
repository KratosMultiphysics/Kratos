//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2016 $
//   Revision:            $Revision:                    0.0 $
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
#include "custom_processes/contact_model_start_end_meshing_process.hpp"
#include "custom_processes/parametric_wall_contact_search_process.hpp"

namespace Kratos
{
	
  namespace Python
  {
 	
    void  AddCustomProcessesToPython()
    {

      using namespace boost::python;
      typedef Process                                              ProcessBaseType;
      typedef ModelStartEndMeshingProcess      ModelStartEndMeshingProcessBaseType;

      
      //**********MESH MODELLER PROCESS*********//

      class_<ContactModelStartEndMeshingProcess, bases< ModelStartEndMeshingProcessBaseType >, boost::noncopyable >
	("ContactModelMeshing", init<ModelPart&, Flags, int>())
	;

      class_<ParametricWallContactSearchProcess, bases< ProcessBaseType >, boost::noncopyable >
	("ParametricWallContactSearch", init<ModelPart&, SpatialBoundingBox::Pointer, Parameters>())
	;

    }
 
  }  // namespace Python.

} // Namespace Kratos

