---
type: Python Module
title: test_vtu_output
resource: tests/test_vtu_output.py#L1-L641
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/math
    resolved_by: tree-sitter
    confidence: exact
  - target: external/typing
    resolved_by: tree-sitter
    confidence: exact
  - target: external/pathlib
    resolved_by: tree-sitter
    confidence: exact
  - target: external/json
    resolved_by: tree-sitter
    confidence: exact
  - target: external/filecmp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/xml-etree-elementtree-as-et
    resolved_by: tree-sitter
    confidence: exact
  - target: external/kratosmultiphysics-as-kratos
    resolved_by: tree-sitter
    confidence: exact
  - target: external/kratosmultiphysics-kratosunittest-as-kratos-unittest
    resolved_by: tree-sitter
    confidence: exact
  - target: external/kratosmultiphysics-kratos-utilities-as-kratos-utils
    resolved_by: tree-sitter
    confidence: exact
  - target: external/kratosmultiphysics-compare-two-files-check-process
    resolved_by: tree-sitter
    confidence: exact
  - target: external/test-vtk-output-process
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [TestVtuOutputBase](../../classes/tests/test_vtu_output/TestVtuOutputBase.md)
- [SetSolution](../../functions/tests/test_vtu_output/TestVtuOutputBase/SetSolution.md)
- [setUpClass](../../functions/tests/test_vtu_output/TestVtuOutputBase/setUpClass.md)
- [WriteVtu](../../functions/tests/test_vtu_output/TestVtuOutputBase/WriteVtu.md)
- [test_WriteMeshAscii](../../functions/tests/test_vtu_output/TestVtuOutputBase/test_WriteMeshAscii.md)
- [test_WriteMeshBinary](../../functions/tests/test_vtu_output/TestVtuOutputBase/test_WriteMeshBinary.md)
- [test_WriteMeshRaw](../../functions/tests/test_vtu_output/TestVtuOutputBase/test_WriteMeshRaw.md)
- [test_WriteMeshCompressedRaw](../../functions/tests/test_vtu_output/TestVtuOutputBase/test_WriteMeshCompressedRaw.md)
- [test_WriteMeshAsciiWithError](../../functions/tests/test_vtu_output/TestVtuOutputBase/test_WriteMeshAsciiWithError.md)
- [test_WriteMeshBinaryWithError](../../functions/tests/test_vtu_output/TestVtuOutputBase/test_WriteMeshBinaryWithError.md)
- [test_WriteMeshRawWithError](../../functions/tests/test_vtu_output/TestVtuOutputBase/test_WriteMeshRawWithError.md)
- [test_WriteMeshCompressedRawWithError](../../functions/tests/test_vtu_output/TestVtuOutputBase/test_WriteMeshCompressedRawWithError.md)
- [Check](../../functions/tests/test_vtu_output/TestVtuOutputBase/Check.md)
- [check_file](../../functions/tests/test_vtu_output/check_file.md)
- [TestVtuOutput2D](../../classes/tests/test_vtu_output/TestVtuOutput2D.md)
- [setUpClass](../../functions/tests/test_vtu_output/TestVtuOutput2D/setUpClass.md)
- [tearDown](../../functions/tests/test_vtu_output/TestVtuOutput2D/tearDown.md)
- [TestVtuOutput3D](../../classes/tests/test_vtu_output/TestVtuOutput3D.md)
- [setUpClass](../../functions/tests/test_vtu_output/TestVtuOutput3D/setUpClass.md)
- [tearDown](../../functions/tests/test_vtu_output/TestVtuOutput3D/tearDown.md)
- [TestVtuOutput](../../classes/tests/test_vtu_output/TestVtuOutput.md)
- [setUpClass](../../functions/tests/test_vtu_output/TestVtuOutput/setUpClass.md)
- [fill_nodes](../../functions/tests/test_vtu_output/fill_nodes.md)
- [add_variables](../../functions/tests/test_vtu_output/add_variables.md)
- [setUp](../../functions/tests/test_vtu_output/TestVtuOutput/setUp.md)
- [test_GetOutputContainerList](../../functions/tests/test_vtu_output/TestVtuOutput/test_GetOutputContainerList.md)
- [test_PointVariableAddition](../../functions/tests/test_vtu_output/TestVtuOutput/test_PointVariableAddition.md)
- [test_CellVariableAddition](../../functions/tests/test_vtu_output/TestVtuOutput/test_CellVariableAddition.md)
- [GetUnstructuredGridList](../../functions/tests/test_vtu_output/TestVtuOutput/GetUnstructuredGridList.md)
- [GetUnstructuredGridName](../../functions/tests/test_vtu_output/TestVtuOutput/GetUnstructuredGridName.md)
- [GetNodes](../../functions/tests/test_vtu_output/TestVtuOutput/GetNodes.md)
- [GetNumberOfNodes](../../functions/tests/test_vtu_output/TestVtuOutput/GetNumberOfNodes.md)
- [CheckDataArray](../../functions/tests/test_vtu_output/TestVtuOutput/CheckDataArray.md)
- [CheckVtuFile](../../functions/tests/test_vtu_output/TestVtuOutput/CheckVtuFile.md)
- [CheckVtmFile](../../functions/tests/test_vtu_output/TestVtuOutput/CheckVtmFile.md)
- [CheckSeriesFile](../../functions/tests/test_vtu_output/TestVtuOutput/CheckSeriesFile.md)
- [CheckGaussVtuFile](../../functions/tests/test_vtu_output/TestVtuOutput/CheckGaussVtuFile.md)
- [__GetUnstructuredGridList](../../functions/tests/test_vtu_output/TestVtuOutput/__GetUnstructuredGridList.md)
- [tearDownClass](../../functions/tests/test_vtu_output/TestVtuOutput/tearDownClass.md)

# Imports

- `math`
- `typing`
- `pathlib`
- `json`
- `filecmp`
- `xml.etree.ElementTree as ET`
- `KratosMultiphysics as Kratos`
- `KratosMultiphysics.KratosUnittest as kratos_unittest`
- `KratosMultiphysics.kratos_utilities as kratos_utils`
- `KratosMultiphysics.compare_two_files_check_process`
- `test_vtk_output_process`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)