---
description: "Use when writing or reviewing Python scripts or pybind11 bindings in Kratos Multiphysics core or applications. Covers PEP 8 style, Kratos imports, KratosUnittest, and pybind11 binding conventions."
applyTo: ["**/*.py", "**/custom_python/**/*.cpp"]
---

# Python and pybind11 Conventions — Kratos Multiphysics

## Python Style

- Follow **PEP 8** for all new Python code.
- Use `snake_case` for all function names, variable names, and module names.
- Import from `KratosMultiphysics` and the specific application module:

```python
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanics
# Or the relevant application, e.g.:
import KratosMultiphysics.FluidDynamicsApplication as FluidDynamics
```

- Python analysis stages and processes must inherit from the appropriate Kratos base classes.
- Do **not** write Python tests without inheriting from `KratosUnittest.TestCase` (see `testing.instructions.md`).

## pybind11 Bindings

Bindings live under each application's `custom_python/` directory (or `kratos/python/` for core).

- Binding files follow the naming pattern: `add_custom_<type>_to_python.cpp`
  (e.g., `add_custom_processes_to_python.cpp`, `add_custom_elements_to_python.cpp`).
- Keep binding files **thin** — business logic belongs in C++ classes, not in binding code.
- Expose classes with `py::class_<Derived, Derived::Pointer, Base>`, including the correct base class and pointer type.
- Maintain naming and signature consistency with neighboring binding files.
- Provide docstrings for all **newly introduced** public bindings; avoid rewriting unrelated legacy ones.
- Specify `py::return_value_policy` explicitly when returning references or pointers.
- For Processes, expose `Execute()` and any relevant stage methods (`ExecuteInitialize`, `ExecuteInitializeSolutionStep`, etc.).

```cpp
namespace py = pybind11;

void AddCustomProcessesToPython(py::module& m) {
    py::class_<MyCustomProcess, MyCustomProcess::Pointer, Process>(m, "MyCustomProcess")
        .def(py::init<ModelPart&, Parameters>(),
             py::arg("model_part"), py::arg("parameters"),
             "Initializes MyCustomProcess with a ModelPart and JSON parameters.")
        .def("Execute", &MyCustomProcess::Execute,
             "Executes the process.")
        .def("ExecuteInitialize", &MyCustomProcess::ExecuteInitialize,
             "Runs the initialization stage.")
        .def("ExecuteInitializeSolutionStep", &MyCustomProcess::ExecuteInitializeSolutionStep,
             "Runs the per-step initialization stage.")
        ;
}
```

The main module file (`<app>_python_application.cpp`) calls all `AddCustom*ToPython(m)` functions and registers application variables.
