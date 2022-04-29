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
#include "custom_python/add_custom_responses_to_python.h"
#include "custom_python/add_custom_optimization_algorithm_to_python.h"
#include "custom_python/add_custom_strategies_to_python.h"

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

    AddCustomResponsesToPython(m);
    AddCustomControlsToPython(m);
    AddCustomOptimizationAlgorithmToPython(m);
    AddCustomStrategiesToPython(m);

    //registering variables in python

    // Optimization variables

    //Auxilary field
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, AUXILIARY_FIELD);

    //linear function
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, D_LINEAR_D_X);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, D_LINEAR_D_CX); 

    //symmetry plane
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, D_PLANE_SYMMETRY_D_X);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, D_PLANE_SYMMETRY_D_CX);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, NEAREST_NEIGHBOUR_POINT); 
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NEAREST_NEIGHBOUR_DIST);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NEAREST_NEIGHBOUR_COND_ID);

    //strain energy
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, D_STRAIN_ENERGY_1_D_X);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, D_STRAIN_ENERGY_1_D_CX);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, D_STRAIN_ENERGY_2_D_X);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, D_STRAIN_ENERGY_2_D_CX);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, D_STRAIN_ENERGY_3_D_X);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, D_STRAIN_ENERGY_3_D_CX);        
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, D_STRAIN_ENERGY_D_PT);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, D_STRAIN_ENERGY_D_CT);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, D_STRAIN_ENERGY_D_PD);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, D_STRAIN_ENERGY_D_CD);   

    //mass
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, D_MASS_D_X);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, D_MASS_D_CX);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, D_MASS_D_PT);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, D_MASS_D_CT);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, D_MASS_D_PD);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, D_MASS_D_CD);   

    //eigenfrequency
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, D_EIGEN_FREQ_D_X);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, D_EIGEN_FREQ_D_CX);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, D_EIGEN_FREQ_D_PT);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, D_EIGEN_FREQ_D_CT);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, D_EIGEN_FREQ_D_PD);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, D_EIGEN_FREQ_D_CD);       

    //local_stress
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, D_LOCAL_STRESS_D_X);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, D_LOCAL_STRESS_D_CX);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, D_LOCAL_STRESS_D_PT);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, D_LOCAL_STRESS_D_CT);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, D_LOCAL_STRESS_D_PD);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, D_LOCAL_STRESS_D_CD);     
    

    //max_stress
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, D_MAX_STRESS_D_X);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, D_MAX_STRESS_D_CX);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, D_MAX_STRESS_D_PT);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, D_MAX_STRESS_D_CT);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, D_MAX_STRESS_D_PD);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, D_MAX_STRESS_D_CD);  
        

    // shape control
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, CX);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, D_CX);  
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, D_X); 

    // thickness control
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PT);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CT);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, FT);    
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, D_CT);  
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, D_PT);  

    // topology control
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CD);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PD);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, FD);    
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, D_CD);  
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, D_PD);          

    // For implicit vertex-morphing with Helmholtz PDE
	  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, HELMHOLTZ_MASS_MATRIX );
	  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, HELMHOLTZ_SURF_RADIUS_SHAPE );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, HELMHOLTZ_BULK_RADIUS_SHAPE );
	  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, COMPUTE_CONTROL_POINTS_SHAPE );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, HELMHOLTZ_SURF_POISSON_RATIO_SHAPE );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, HELMHOLTZ_BULK_POISSON_RATIO_SHAPE );
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, HELMHOLTZ_VARS_SHAPE);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, HELMHOLTZ_SOURCE_SHAPE);  

    // for thickness optimization
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, HELMHOLTZ_VAR_THICKNESS );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, HELMHOLTZ_SOURCE_THICKNESS ); 
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, HELMHOLTZ_RADIUS_THICKNESS );   

    // for topology optimization
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, HELMHOLTZ_VAR_DENSITY );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, HELMHOLTZ_SOURCE_DENSITY ); 
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, HELMHOLTZ_RADIUS_DENSITY );         
       

  }

}  // namespace Python.
}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
