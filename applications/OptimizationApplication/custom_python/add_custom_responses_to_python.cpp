// ==============================================================================
//  KratosOptimizationApplication
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
//
// =================================================================================

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------

// // ------------------------------------------------------------------------------
// // External includes
// // ------------------------------------------------------------------------------
#include <pybind11/stl.h>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"
#include "custom_python/add_custom_controls_to_python.h"
#include "custom_responses/shape_responses/plane_symmetry.h"
#include "custom_responses/shape_responses/linear.h"
#include "custom_responses/mass_opt_response.h"
#include "custom_responses/linear_strain_energy_opt_response.h"

// ==============================================================================

namespace Kratos {
namespace Python {



// ==============================================================================
void  AddCustomResponsesToPython(pybind11::module& m)
{
    namespace py = pybind11;
    // ================================================================
    // 
    // ================================================================
    py::class_<PlaneSymmetry >(m, "PlaneSymmetry")
        .def(py::init<std::string, Model&, Parameters&>())
        .def("Initialize", &PlaneSymmetry::Initialize)
        .def("CalculateValue", &PlaneSymmetry::CalculateValue)
        .def("CalculateGradient", &PlaneSymmetry::CalculateGradient)        
        ;     

    py::class_<Linear >(m, "Linear")
        .def(py::init<std::string, Model&, Parameters&>())
        .def("Initialize", &Linear::Initialize)
        .def("CalculateValue", &Linear::CalculateValue)
        .def("CalculateGradient", &Linear::CalculateGradient)        
        ; 

    py::class_<MassOptResponse >(m, "MassOptResponse")
        .def(py::init<std::string, Model&, Parameters&>())
        .def("Initialize", &MassOptResponse::Initialize)
        .def("CalculateValue", &MassOptResponse::CalculateValue)
        .def("CalculateGradient", &MassOptResponse::CalculateGradient)        
        ;  

    py::class_<LinearStrainEnergyOptResponse >(m, "LinearStrainEnergyOptResponse")
        .def(py::init<std::string, Model&, Parameters&>())
        .def("Initialize", &LinearStrainEnergyOptResponse::Initialize)
        .def("CalculateValue", &LinearStrainEnergyOptResponse::CalculateValue)
        .def("CalculateGradient", &LinearStrainEnergyOptResponse::CalculateGradient)        
        ;                                  
 
}

}  // namespace Python.
} // Namespace Kratos

