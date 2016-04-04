//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 1.0 $
//

// External includes
#include <boost/python.hpp>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "spaces/ublas_space.h"

#include "custom_processes/displacement_table_interpolation_process.hpp"
#include "custom_processes/pressure_table_interpolation_process.hpp"
#include "custom_processes/force_table_interpolation_process.hpp"
#include "custom_processes/face_load_table_interpolation_process.hpp"
#include "custom_processes/normal_load_table_interpolation_process.hpp"
#include "custom_processes/tangential_load_table_interpolation_process.hpp"
#include "custom_processes/normal_flux_table_interpolation_process.hpp"


namespace Kratos
{
	
namespace Python
{

using namespace boost::python;

void  AddCustomProcessesToPython() 
{    
   
    class_< DisplacementTableInterpolationProcess, bases< Process >, boost::noncopyable > ( "DisplacementTableInterpolationProcess",
        init < ModelPart& >());
        
    class_< PressureTableInterpolationProcess, bases< Process >, boost::noncopyable > ( "PressureTableInterpolationProcess",
        init < ModelPart& >());

    class_< ForceTableInterpolationProcess, bases< Process >, boost::noncopyable > ( "ForceTableInterpolationProcess",
        init < ModelPart& >());

    class_< FaceLoadTableInterpolationProcess, bases< Process >, boost::noncopyable > ( "FaceLoadTableInterpolationProcess",
        init < ModelPart& >());

    class_< NormalLoadTableInterpolationProcess, bases< Process >, boost::noncopyable > ( "NormalLoadTableInterpolationProcess",
        init < ModelPart& >());

    class_< TangentialLoadTableInterpolationProcess, bases< Process >, boost::noncopyable > ( "TangentialLoadTableInterpolationProcess",
        init < ModelPart& >());

    class_< NormalFluxTableInterpolationProcess, bases< Process >, boost::noncopyable > ( "NormalFluxTableInterpolationProcess",
        init < ModelPart& >());
}

}  // namespace Python.
} // Namespace Kratos
