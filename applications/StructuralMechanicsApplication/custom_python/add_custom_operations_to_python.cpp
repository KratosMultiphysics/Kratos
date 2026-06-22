// --- Structural Includes ---
#include "custom_python/add_custom_operations_to_python.hpp" // AddCustomOperationsToPython
#include "custom_operations/insert_pre_tension_operation.hpp" // InsertPreTensionOperation


namespace Kratos::Python {


void AddCustomOperationsToPython(pybind11::module& rModule) {
    pybind11::class_<
        InsertPreTensionOperation,
        InsertPreTensionOperation::Pointer,
        Operation>(
            rModule,
            "InsertPreTensionOperation")
                .def(pybind11::init<Model&, Parameters>())
                ;
}


} // namespace Kratos::Python
