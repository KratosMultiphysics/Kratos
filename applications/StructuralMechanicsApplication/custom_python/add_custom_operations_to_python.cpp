// --- Structural Includes ---
#include "custom_python/add_custom_operations_to_python.hpp" // AddCustomOperationsToPython
#include "custom_operations/insert_pre_tension_operation.hpp" // InsertNeumannPreTensionOperation, InsertDirichletPreTensionOperation


namespace Kratos::Python {


void AddCustomOperationsToPython(pybind11::module& rModule) {
    pybind11::class_<
        InsertDirichletPreTensionOperation,
        InsertDirichletPreTensionOperation::Pointer,
        Operation>(
            rModule,
            "InsertDirichletPreTensionOperation")
                .def(pybind11::init<Model&, Parameters>())
                ;
    pybind11::class_<
        InsertNeumannPreTensionOperation,
        InsertNeumannPreTensionOperation::Pointer,
        Operation>(
            rModule,
            "InsertNeumannPreTensionOperation")
                .def(pybind11::init<Model&, Parameters>())
                ;
}


} // namespace Kratos::Python
