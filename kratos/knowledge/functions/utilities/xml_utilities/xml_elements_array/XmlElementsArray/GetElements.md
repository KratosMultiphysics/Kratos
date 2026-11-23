---
type: C++ Method
title: GetElements
resource: utilities/xml_utilities/xml_elements_array.cpp#L32-L35
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/InitializeDatabases
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/CheckGPValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_xml_utilities/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_cad_tessellation_modeler/TestCadTessellationModeler/test_cad_tessellation_modeler
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::vector<XmlElement::Pointer> XmlElementsArray::GetElements() const`

# Called by

- [InitializeDatabases](../../../../../functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/InitializeDatabases.md)
- [CheckGPValues](../../../../../functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/CheckGPValues.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../../functions/tests/cpp_tests/utilities/test_xml_utilities/KRATOS_TEST_CASE_IN_SUITE.md)
- [test_cad_tessellation_modeler](../../../../../functions/tests/test_cad_tessellation_modeler/TestCadTessellationModeler/test_cad_tessellation_modeler.md)