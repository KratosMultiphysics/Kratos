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
#include "custom_responses/stress_opt_response.h"
#include "custom_responses/interface_opt_response.h"
#include "custom_responses/partition_mass_opt_response.h"
#include "custom_responses/partition_interface_stress_opt_response.h"

// ==============================================================================

namespace Kratos {
namespace Python {



// ==============================================================================
void  AddCustomResponsesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;    
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;

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

    py::class_<InterfaceOptResponse >(m, "InterfaceOptResponse")
        .def(py::init<std::string, Model&, Parameters&>())
        .def("Initialize", &InterfaceOptResponse::Initialize)
        .def("CalculateValue", &InterfaceOptResponse::CalculateValue)
        .def("CalculateGradient", &InterfaceOptResponse::CalculateGradient)        
        ;     

    py::class_<PartitionMassOptResponse >(m, "PartitionMassOptResponse")
        .def(py::init<std::string, Model&, Parameters&>())
        .def("Initialize", &PartitionMassOptResponse::Initialize)
        .def("CalculateValue", &PartitionMassOptResponse::CalculateValue)
        .def("CalculateGradient", &PartitionMassOptResponse::CalculateGradient)        
        ;          

    py::class_<LinearStrainEnergyOptResponse >(m, "LinearStrainEnergyOptResponse")
        .def(py::init<std::string, Model&, Parameters&>())
        .def("Initialize", &LinearStrainEnergyOptResponse::Initialize)
        .def("CalculateValue", &LinearStrainEnergyOptResponse::CalculateValue)
        .def("CalculateGradient", &LinearStrainEnergyOptResponse::CalculateGradient)        
        ;   

    py::class_<StressOptResponse >(m, "StressOptResponse")
        .def(py::init<std::string, Model&, Parameters&, std::vector<LinearSolverType::Pointer>&>())
        .def("Initialize", &StressOptResponse::Initialize)
        .def("CalculateValue", &StressOptResponse::CalculateValue)
        .def("CalculateGradient", &StressOptResponse::CalculateGradient)        
        ;              

    py::class_<PartitionInterfaceStressOptResponse >(m, "PartitionInterfaceStressOptResponse")
        .def(py::init<std::string, Model&, Parameters&>())
        .def("Initialize", &PartitionInterfaceStressOptResponse::Initialize)
        .def("CalculateValue", &PartitionInterfaceStressOptResponse::CalculateValue)
        .def("CalculateGradient", &PartitionInterfaceStressOptResponse::CalculateGradient)        
        ;                                      
 
}

}  // namespace Python.
} // Namespace Kratos

