---
description: "Use when writing or reviewing Python scripts or pybind11 bindings in applications/GeoMechanicsApplication. Covers Python process wrappers, analysis stage registration, import aliases, and the binding style of custom_python/."
applyTo: ["applications/GeoMechanicsApplication/**/*.py", "applications/GeoMechanicsApplication/custom_python/**/*.cpp"]
---

# Python and binding conventions: GeoMechanicsApplication

The general rules in `/agents/instructions/python-conventions.instructions.md` apply (PEP 8,
`snake_case`, thin bindings). This file adds the GeoMechanics specifics.

## Python modules

- Python modules live in `python_scripts/` and are installed as
  `KratosMultiphysics.GeoMechanicsApplication.<module>` by `CMakeLists.txt`. No CMake change is
  needed for a new module.
- Imports: most modules use `import KratosMultiphysics` together with
  `import KratosMultiphysics.GeoMechanicsApplication as KratosGeo`. Some newer modules use
  `as Core` and `as Geo`. Use the aliases the file you edit already uses; in a new file use
  `KratosMultiphysics` and `KratosGeo`.
- The analysis stage is `geomechanics_analysis.GeoMechanicsAnalysis`, registered in
  `python_registry_lists.py` (`python_stages_to_be_registered`). Register new Python modelers,
  operations, processes, stages or orchestrators in the matching list of that file.
- Solvers: `geomechanics_U_Pw_solver.py`, `geomechanics_Pw_solver.py`,
  `geomechanics_T_solver.py`, created through `geomechanics_solvers_wrapper.py`. A new solver
  type needs an entry in the wrapper.

## Python process wrappers

A C++ process is used from `ProjectParameters.json` through a small module in `python_scripts/`
whose name is the `snake_case` form of the process name. It has one `Factory` function, as in
`apply_k0_procedure_process.py`:

```python
import KratosMultiphysics as Core
import KratosMultiphysics.GeoMechanicsApplication as Geo


def Factory(settings, model):
    if not isinstance(settings, Core.Parameters):
        raise TypeError("expected input shall be a Parameters object, encapsulating a json string")

    return Geo.ApplyK0ProcedureProcess(model, settings["Parameters"])
```

The JSON entry that uses it:

```json
{
    "python_module": "apply_k0_procedure_process",
    "kratos_module": "KratosMultiphysics.GeoMechanicsApplication",
    "process_name": "ApplyK0ProcedureProcess",
    "Parameters": { ... }
}
```

Keep the wrapper free of logic. Parameter validation and defaults belong in the C++ constructor.

## pybind11 bindings (`custom_python/`)

- Processes go in `add_custom_processes_to_python.cpp`, utilities in
  `add_custom_utilities_to_python.cpp`, strategies, schemes and builders in
  `add_custom_strategies_to_python.cpp`, retention laws in `add_retention_laws_to_python.cpp`.
  Variables are exposed in `geo_mechanics_python_application.cpp` with
  `KRATOS_REGISTER_IN_PYTHON_VARIABLE`.
- Processes bind with the base class and the constructor that the Python wrapper calls:

  ```cpp
  py::class_<ApplyK0ProcedureProcess, ApplyK0ProcedureProcess::Pointer, Process>(
      m, "ApplyK0ProcedureProcess")
      .def(py::init<Model&, Parameters>());
  ```

  Most GeoMechanics processes take `Model&` (or `ModelPart&`) and `Parameters`; match the
  C++ constructor.
- Docstrings: the root convention asks for docstrings on newly introduced bindings. None of the
  current GeoMechanics bindings has one. Add a docstring to new bindings; do not flag the
  existing ones in reviews.

