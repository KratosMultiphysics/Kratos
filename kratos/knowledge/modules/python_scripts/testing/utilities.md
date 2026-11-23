---
type: Python Module
title: utilities
resource: python_scripts/testing/utilities.py#L1-L335
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/kratosmultiphysics-as-kratos
    resolved_by: tree-sitter
    confidence: exact
  - target: external/kratosmultiphysics-kratos-utilities-as-kratos-utils
    resolved_by: tree-sitter
    confidence: exact
  - target: external/kratosmultiphysics
    resolved_by: tree-sitter
    confidence: exact
  - target: external/os
    resolved_by: tree-sitter
    confidence: exact
  - target: external/sys
    resolved_by: tree-sitter
    confidence: exact
  - target: external/subprocess
    resolved_by: tree-sitter
    confidence: exact
  - target: external/pathlib
    resolved_by: tree-sitter
    confidence: exact
  - target: external/kratosmultiphysics-mpi
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [GetPython3Command](../../../functions/python_scripts/testing/utilities/GetPython3Command.md)
- [ReadModelPart](../../../functions/python_scripts/testing/utilities/ReadModelPart.md)
- [ReadSerialModelPart](../../../functions/python_scripts/testing/utilities/ReadSerialModelPart.md)
- [ReadDistributedModelPart](../../../functions/python_scripts/testing/utilities/ReadDistributedModelPart.md)
- [PrintTestHeader](../../../functions/python_scripts/testing/utilities/PrintTestHeader.md)
- [PrintTestFooter](../../../functions/python_scripts/testing/utilities/PrintTestFooter.md)
- [PrintTestSummary](../../../functions/python_scripts/testing/utilities/PrintTestSummary.md)
- [Commander](../../../classes/python_scripts/testing/utilities/Commander.md)
- [__init__](../../../functions/python_scripts/testing/utilities/Commander/__init__.md)
- [PrintOutput](../../../functions/python_scripts/testing/utilities/Commander/PrintOutput.md)
- [TestToAppName](../../../functions/python_scripts/testing/utilities/Commander/TestToAppName.md)
- [MPITestToAppName](../../../functions/python_scripts/testing/utilities/Commander/MPITestToAppName.md)
- [_RunTest](../../../functions/python_scripts/testing/utilities/Commander/_RunTest.md)
- [RunPythonTests](../../../functions/python_scripts/testing/utilities/Commander/RunPythonTests.md)
- [RunCppTests](../../../functions/python_scripts/testing/utilities/Commander/RunCppTests.md)
- [RunMPIPythonTests](../../../functions/python_scripts/testing/utilities/Commander/RunMPIPythonTests.md)
- [RunMPICppTests](../../../functions/python_scripts/testing/utilities/Commander/RunMPICppTests.md)

# Imports

- `KratosMultiphysics as Kratos`
- `KratosMultiphysics.kratos_utilities as kratos_utils`
- `KratosMultiphysics`
- `os`
- `sys`
- `subprocess`
- `pathlib`
- `KratosMultiphysics.mpi`

# Member of

- [KratosMultiphysics](../../../packages/KratosMultiphysics.md)