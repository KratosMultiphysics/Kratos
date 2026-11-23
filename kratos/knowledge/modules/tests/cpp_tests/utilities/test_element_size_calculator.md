---
type: C++ Module
title: test_element_size_calculator
resource: tests/cpp_tests/utilities/test_element_size_calculator.cpp#L1-L534
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/limits
    resolved_by: tree-sitter
    confidence: exact
  - target: external/testing-testing-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-triangle-2d-3-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-triangle-2d-6-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-quadrilateral-2d-4-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-quadrilateral-2d-9-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-tetrahedra-3d-4-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-tetrahedra-3d-10-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-prism-3d-6-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-hexahedra-3d-8-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-hexahedra-3d-27-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/tests-cpp-tests-geometries-test-geometry-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-element-size-calculator-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [RunElementSizeCalculatorDerivativesTest](../../../../functions/tests/cpp_tests/utilities/test_element_size_calculator/RunElementSizeCalculatorDerivativesTest.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_element_size_calculator/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_element_size_calculator/KRATOS_TEST_CASE_IN_SUITE-2.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_element_size_calculator/KRATOS_TEST_CASE_IN_SUITE-3.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_element_size_calculator/KRATOS_TEST_CASE_IN_SUITE-4.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_element_size_calculator/KRATOS_TEST_CASE_IN_SUITE-5.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_element_size_calculator/KRATOS_TEST_CASE_IN_SUITE-6.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_element_size_calculator/KRATOS_TEST_CASE_IN_SUITE-7.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_element_size_calculator/KRATOS_TEST_CASE_IN_SUITE-8.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_element_size_calculator/KRATOS_TEST_CASE_IN_SUITE-9.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_element_size_calculator/KRATOS_TEST_CASE_IN_SUITE-10.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_element_size_calculator/KRATOS_TEST_CASE_IN_SUITE-11.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_element_size_calculator/KRATOS_TEST_CASE_IN_SUITE-12.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_element_size_calculator/KRATOS_TEST_CASE_IN_SUITE-13.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_element_size_calculator/KRATOS_TEST_CASE_IN_SUITE-14.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_element_size_calculator/KRATOS_TEST_CASE_IN_SUITE-15.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_element_size_calculator/KRATOS_TEST_CASE_IN_SUITE-16.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_element_size_calculator/KRATOS_TEST_CASE_IN_SUITE-17.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_element_size_calculator/KRATOS_TEST_CASE_IN_SUITE-18.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_element_size_calculator/KRATOS_TEST_CASE_IN_SUITE-19.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_element_size_calculator/KRATOS_TEST_CASE_IN_SUITE-20.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_element_size_calculator/KRATOS_TEST_CASE_IN_SUITE-21.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_element_size_calculator/KRATOS_TEST_CASE_IN_SUITE-22.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_element_size_calculator/KRATOS_TEST_CASE_IN_SUITE-23.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_element_size_calculator/KRATOS_TEST_CASE_IN_SUITE-24.md)

# Imports

- `<limits>`
- `"testing/testing.h"`
- `"geometries/triangle_2d_3.h"`
- `"geometries/triangle_2d_6.h"`
- `"geometries/quadrilateral_2d_4.h"`
- `"geometries/quadrilateral_2d_9.h"`
- `"geometries/tetrahedra_3d_4.h"`
- `"geometries/tetrahedra_3d_10.h"`
- `"geometries/prism_3d_6.h"`
- `"geometries/hexahedra_3d_8.h"`
- `"geometries/hexahedra_3d_27.h"`
- `"tests/cpp_tests/geometries/test_geometry.h"`
- `"utilities/element_size_calculator.h"`

# Member of

- [KratosMultiphysics](../../../../packages/KratosMultiphysics.md)