---
description: "Use when writing or reviewing C++ code in Kratos Multiphysics core or applications. Covers Kratos naming conventions, error handling macros, memory management, and patterns for Elements, Conditions, and Processes."
applyTo: ["**/*.cpp", "**/*.h", "**/*.hpp"]
---

# C++ Coding Conventions — Kratos Multiphysics

Follow Kratos conventions as used throughout `kratos/` and `applications/`.

## Naming

| Symbol | Convention | Example |
|--------|-----------|---------|
| Classes / types | `PascalCase` | `TotalLagrangianElement` |
| Kratos interface methods | `PascalCase` | `CalculateLocalSystem`, `GetDofList` |
| Free functions / local vars | `snake_case` | `compute_volume` |
| Member variables | `m` prefix | `mThickness`, `mConstitutiveLaw` |
| Reference parameters | `r` prefix | `rModelPart`, `rCurrentProcessInfo` |
| Pointer parameters | `p` prefix | `pElement`, `pNode` |
| Constants / Kratos variables | follow existing project style | `DISPLACEMENT`, `TEMPERATURE` |

## File Layout

- Always use `#pragma once` in header files (no macro guards).
- Always initialize variables at declaration.
- Use `const` aggressively.
- Prefer `const auto&` for range-based loops over containers.

## Error Handling and Logging

Use Kratos macros — **never** use `std::cout` in production code.

```cpp
void MyProcess::Execute() {
    KRATOS_TRY

    KRATOS_ERROR_IF(mrModelPart.NumberOfNodes() == 0)
        << "ModelPart has no nodes." << std::endl;

    KRATOS_INFO("MyProcess") << "Starting execution." << std::endl;
    KRATOS_WARNING("MyProcess") << "Unusual condition encountered." << std::endl;

    // ... implementation ...

    KRATOS_CATCH("")
}
```

| Macro | Purpose |
|-------|---------|
| `KRATOS_TRY` / `KRATOS_CATCH("")` | Wrap method bodies for exception context |
| `KRATOS_ERROR_IF(cond)` | Throw on condition |
| `KRATOS_ERROR_IF_NOT(cond)` | Throw if condition is false |
| `KRATOS_ERROR` | Unconditional throw |
| `KRATOS_INFO("tag")` | Info-level logging |
| `KRATOS_WARNING("tag")` | Warning-level logging |

## Memory Management

- **Never** use raw `new` / `delete`.
- Use `Kratos::shared_ptr`, `Kratos::unique_ptr`, and Kratos pointer macros/type aliases.
- Use Kratos type aliases instead of raw STL/primitive types where available: `IndexType`, `SizeType`, `MatrixType`, `VectorType`, etc.
- Do **not** use `using namespace std;` globally.

## Elements and Conditions

Always override the required virtual methods:

```cpp
class MyCustomElement : public Element {
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(MyCustomElement);

    // Required overrides
    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                               VectorType& rRightHandSideVector,
                               const ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdVector(EquationIdVectorType& rResult,
                          const ProcessInfo& rCurrentProcessInfo) const override;

    void GetDofList(DofsVectorType& rElementalDofList,
                    const ProcessInfo& rCurrentProcessInfo) const override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;
};
```

## Processes

Inherit from `Process` and implement the appropriate stage methods:

```cpp
class MyCustomProcess : public Process {
public:
    KRATOS_CLASS_POINTER_DEFINITION(MyCustomProcess);

    MyCustomProcess(ModelPart& rModelPart, Parameters rParameters)
        : mrModelPart(rModelPart), mParameters(rParameters) {}

    void ExecuteInitialize() override;
    void ExecuteBeforeSolutionLoop() override;
    void ExecuteInitializeSolutionStep() override;
    void Execute() override;  // for single-shot processes
    void ExecuteFinalizeSolutionStep() override;
    void ExecuteFinalize() override;

private:
    ModelPart& mrModelPart;
    Parameters mParameters;
};
```

- Prefer **data-driven design** using Kratos `Parameters` (JSON-based) for process configuration.
- Register new components in the application's `RegisterComponents()` method inside `<app>_application.cpp`.
