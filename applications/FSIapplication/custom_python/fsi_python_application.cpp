//
//   Project Name:        Kratos
//   Last modified by:    $Author: rrossi $
//   Date:                $Date: 2008-03-25 07:48:22 $
//   Revision:            $Revision: 1.5 $
//
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "fsi_application.h"
#include "custom_python/add_convergence_accelerators_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_mappers_to_python.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;



BOOST_PYTHON_MODULE(KratosFSIApplication)
{

    class_<KratosFSIApplication,
           KratosFSIApplication::Pointer,
           bases<KratosApplication>, boost::noncopyable >("KratosFSIApplication")
           ;

    AddCustomUtilitiesToPython();
    AddMappersToPython();
    AddConvergenceAcceleratorsToPython();

    //registering variables in python
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(CONVERGENCE_ACCELERATOR_ITERATION);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(MAPPER_SCALAR_PROJECTION_RHS);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(FICTITIOUS_FLUID_DENSITY);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(SCALAR_PROJECTED);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(FSI_INTERFACE_RESIDUAL_NORM);
//	KRATOS_REGISTER_IN_PYTHON_VARIABLE(PRESSURE_OLD_IT);
//	KRATOS_REGISTER_IN_PYTHON_VARIABLE(IS_INTERFACE);

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(MAPPER_VECTOR_PROJECTION_RHS);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(VAUX_EQ_TRACTION);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(VECTOR_PROJECTED);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(RELAXED_DISP);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(FSI_INTERFACE_RESIDUAL);

}


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
