//
//  Main authors:    Miguel Angel Celigueta   maceli@cimne.upc.edu
//
//


// System includes

// External includes

// Project includes
#include "custom_python/add_custom_processes_to_python.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "custom_processes/advance_and_rotate_part_of_structure_process.h"

namespace Kratos
{

namespace Python
{

void AddCustomProcessesToPython(pybind11::module& m)
{
    using namespace pybind11;

    pybind11::class_<AdvanceAndRotatePartOfStructureProcess, AdvanceAndRotatePartOfStructureProcess::Pointer, Process>
    (m,"AdvanceAndRotatePartOfStructureProcess")
    .def(pybind11::init < ModelPart&, const double, const double, const double, const double, const double, const double, const double >())
    .def(pybind11::init< ModelPart&, Parameters& >())
    ;

}

} // namespace Python.

} // Namespace Kratos
