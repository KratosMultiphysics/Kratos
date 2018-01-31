//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Andreas Winterstein (a.winterstein@tum.de)
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <boost/python.hpp>

// Project includes
#include "includes/define.h"
#include "ale_application.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;



BOOST_PYTHON_MODULE(KratosALEApplication)
{

    class_<KratosALEApplication,
           KratosALEApplication::Pointer,
           bases<KratosApplication>, boost::noncopyable >("KratosALEApplication")
           ;

    AddCustomStrategiesToPython();
    AddCustomUtilitiesToPython();

    // variables moved to core
    // MESH_VELOCITY currently put to the core since used in other applications
    // KRATOS_REGISTER_IN_PYTHON_VARIABLE(MESH_VELOCITY)
    // KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(MESH_DISPLACEMENT);
    // KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(MESH_REACTION);
    // KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(MESH_RHS);

}

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
