// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "sph_application_variables.h"

//Processes
#include "custom_processes/neighbours_search_process.h"
#include "custom_processes/compute_kernel_correction_process.h"
#include "custom_processes/compute_volume_process.h"


namespace Kratos::Python {

void  AddCustomProcessesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    /// Processes
    py::class_<NeighboursSearchProcess, NeighboursSearchProcess::Pointer, Process>(m,"NeighboursSearchProcess")
        .def(py::init<ModelPart&, Parameters>())
        ;
    
    py::class_<ComputeKernelCorrectionProcess, ComputeKernelCorrectionProcess::Pointer, Process>(m,"ComputeKernelCorrectionProcess")
        .def(py::init<ModelPart&, Parameters>())
        ;
    
    py::class_<ComputeVolumeProcess, ComputeVolumeProcess::Pointer, Process>(m,"ComputeVolumeProcess")
        .def(py::init<ModelPart&, Parameters>())
        ;

}

}  // namespace Kratos::Python


