/*
 * Author: Miguel Angel Celigueta
 *
 *  maceli@cimne.upc.edu
 */

#if defined(KRATOS_PYTHON)

// System includes
#include <pybind11/pybind11.h>

// External includes

// Project includes
#include "includes/define_python.h"
#include "../dem_structures_coupling_application.h"
#include "../dem_structures_coupling_application_variables.h"
#include "add_custom_utilities_to_python.h"
#include "add_custom_processes_to_python.h"


namespace Kratos
{

namespace Python
{

using namespace pybind11;

PYBIND11_MODULE(KratosDemStructuresCouplingApplication, m)
{
    class_<KratosDemStructuresCouplingApplication, KratosDemStructuresCouplingApplication::Pointer, KratosApplication>(m, "KratosDemStructuresCouplingApplication")
        .def(init<>())
        ;

    //AddCustomStrategiesToPython(m);
    AddCustomUtilitiesToPython(m);
    AddCustomProcessesToPython(m);

    //registering variables in python

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, DEM_SURFACE_LOAD)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, BACKUP_LAST_STRUCTURAL_VELOCITY)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, BACKUP_LAST_STRUCTURAL_DISPLACEMENT)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, SMOOTHED_STRUCTURAL_VELOCITY)

}

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
