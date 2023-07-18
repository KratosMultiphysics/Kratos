

#if defined(KRATOS_PYTHON)

// System includes
#include <pybind11/pybind11.h>

// External includes

// Project includes
#include "includes/define_python.h"
#include "../DEMFEM_volume_coupling_application.h"
#include "../DEMFEM_volume_coupling_application_variables.h"
// #include "add_custom_utilities_to_python.h"
// #include "add_custom_processes_to_python.h"


namespace Kratos
{

namespace Python
{

using namespace pybind11;

PYBIND11_MODULE(FEMDEMVolumeCouplingApplication, m)
{
    class_<FEMDEMVolumeCouplingApplication, FEMDEMVolumeCouplingApplication::Pointer, KratosApplication>(m, "FEMDEMVolumeCouplingApplication")
        .def(init<>())
        ;

    // //AddCustomStrategiesToPython(m);
    // AddCustomUtilitiesToPython(m);
    // AddCustomProcessesToPython(m);

  
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NODAL_COUPLING_WEIGHT )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PARTICLE_COUPLING_WEIGHT )
}
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, DISPLACEMENT_MULTIPLIED_MASS)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, DEMFEM_VOLUME_COUPLING_FORCE)
  

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
