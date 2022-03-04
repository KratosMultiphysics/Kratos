// ==============================================================================
//  KratosOptimizationApplication
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
//
// =================================================================================

#if defined(KRATOS_PYTHON)

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <pybind11/pybind11.h>

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define_python.h"
#include "optimization_application.h"
#include "optimization_application_variables.h"
#include "custom_python/add_custom_controls_to_python.h"

// ==============================================================================

namespace Kratos {
namespace Python {

PYBIND11_MODULE(KratosOptimizationApplication, m)
{
    namespace py = pybind11;

    py::class_<KratosOptimizationApplication,
        KratosOptimizationApplication::Pointer,
        KratosApplication >(m, "KratosOptimizationApplication")
        .def(py::init<>())
        ;

    AddCustomControlsToPython(m);

    //registering variables in python

    // Optimization variables
    //strain energy
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, D_STRAIN_ENERGY_D_X);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, D_STRAIN_ENERGY_D_CX);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, D_STRAIN_ENERGY_D_T);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, D_STRAIN_ENERGY_D_CT);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, D_STRAIN_ENERGY_D_P);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, D_STRAIN_ENERGY_D_CP);   

    //mass
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, D_MASS_D_X);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, D_MASS_D_CX);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, D_MASS_D_T);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, D_MASS_D_CT);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, D_MASS_D_P);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, D_MASS_D_CP);   

    //eigenfrequency
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, D_EIGEN_FREQ_D_X);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, D_EIGEN_FREQ_D_CX);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, D_EIGEN_FREQ_D_T);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, D_EIGEN_FREQ_D_CT);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, D_EIGEN_FREQ_D_P);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, D_EIGEN_FREQ_D_CP);       

    //local_stress
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, D_LOCAL_STRESS_D_X);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, D_LOCAL_STRESS_D_CX);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, D_LOCAL_STRESS_D_T);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, D_LOCAL_STRESS_D_CT);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, D_LOCAL_STRESS_D_P);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, D_LOCAL_STRESS_D_CP);     
    

    //max_stress
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, D_MAX_STRESS_D_X);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, D_MAX_STRESS_D_CX);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, D_MAX_STRESS_D_T);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, D_MAX_STRESS_D_CT);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, D_MAX_STRESS_D_P);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, D_MAX_STRESS_D_CP);    


    // shape control
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, SHAPE_CONTROL);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, SHAPE_CONTROL_UPDATE);  
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, SHAPE_UPDATE);   
       

  }

}  // namespace Python.
}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
