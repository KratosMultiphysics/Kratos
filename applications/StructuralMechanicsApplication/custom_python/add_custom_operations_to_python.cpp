// --- Structural Includes ---
#include "custom_python/add_custom_operations_to_python.hpp" // AddCustomOperationsToPython
#include "custom_operations/insert_pretension_operation.hpp" // InsertPretensionOperation


namespace Kratos::Python {


void AddCustomOperationsToPython(pybind11::module& rModule) {
    pybind11::class_<InsertPretensionOperation, InsertPretensionOperation::Pointer, Operation>(rModule, "InsertPretensionOperation")
        .def(pybind11::init<Model&, Parameters>())
        ;
}


} // namespace Kratos::Python
