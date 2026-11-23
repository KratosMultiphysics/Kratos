---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: tests/cpp_tests/sources/test_kratos_parameters.cpp#L304-L321
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/cpp_tests/sources/test_kratos_parameters/GetJSONString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/WriteJsonString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetDouble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/PrettyPrintJsonString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_kratos_parameters/GetJSONStringPrettyOut
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/IsBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_kratos_parameters/GetJSONStringPrettyOutAfterChange
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_kratos_parameters/GetJSONStringWrongType
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_kratos_parameters/GetJSONStringDefaults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_kratos_parameters/GetJSONStringWrongSpelling
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_kratos_parameters/GetJSONStringWrongLevel2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/RecursivelyValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_kratos_parameters/GetJSONStringFourLevels
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_kratos_parameters/GetJSONStringForLevelsVariation
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_kratos_parameters/GetJSONStringForLevelsWrongVariation
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_kratos_parameters/GetJSONStringForLevelsDefaults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/IsEquivalentTo
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/HasSameKeysAndTypeOfValuesAs
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetDouble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_kratos_parameters/GetJSONStringExpectedValidationOutput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_kratos_parameters/GetJSONStringIncompleteWithExtraParameter
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddMissingParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/RecursivelyAddMissingParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/ValidateDefaults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_kratos_parameters/GetJSONStringIncomplete
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/RecursivelyValidateDefaults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddEmptyValue
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddEmptyArray
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/RemoveValue
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/IsInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/IsDouble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/IsString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/IsStringArray
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/IsVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/IsMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddDouble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetVectorArray
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetStringArray
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetStringArray
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_kratos_parameters/GetIncludedJSONString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_kratos_parameters/GetIncludedJSONLevel2String
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_kratos_parameters/GetJSONStringWithIncludes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_kratos_parameters/GetCircularIncludeJSONString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/exception/Exception/what
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(KratosParameters, KratosCoreFastSuite)`

# Calls

- [GetJSONString](../../../../../functions/tests/cpp_tests/sources/test_kratos_parameters/GetJSONString.md)
- [WriteJsonString](../../../../../functions/sources/kratos_parameters/Parameters/WriteJsonString.md)
- [GetInt](../../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)
- [GetDouble](../../../../../functions/sources/kratos_parameters/Parameters/GetDouble.md)
- [GetBool](../../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)
- [GetString](../../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [PrettyPrintJsonString](../../../../../functions/sources/kratos_parameters/Parameters/PrettyPrintJsonString.md)
- [GetJSONStringPrettyOut](../../../../../functions/tests/cpp_tests/sources/test_kratos_parameters/GetJSONStringPrettyOut.md)
- [IsBool](../../../../../functions/sources/kratos_parameters/Parameters/IsBool.md)
- [SetString](../../../../../functions/sources/kratos_parameters/Parameters/SetString.md)
- [GetJSONStringPrettyOutAfterChange](../../../../../functions/tests/cpp_tests/sources/test_kratos_parameters/GetJSONStringPrettyOutAfterChange.md)
- [SetInt](../../../../../functions/sources/kratos_parameters/Parameters/SetInt.md)
- [GetJSONStringWrongType](../../../../../functions/tests/cpp_tests/sources/test_kratos_parameters/GetJSONStringWrongType.md)
- [GetJSONStringDefaults](../../../../../functions/tests/cpp_tests/sources/test_kratos_parameters/GetJSONStringDefaults.md)
- [ValidateAndAssignDefaults](../../../../../functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults.md)
- [GetJSONStringWrongSpelling](../../../../../functions/tests/cpp_tests/sources/test_kratos_parameters/GetJSONStringWrongSpelling.md)
- [GetJSONStringWrongLevel2](../../../../../functions/tests/cpp_tests/sources/test_kratos_parameters/GetJSONStringWrongLevel2.md)
- [RecursivelyValidateAndAssignDefaults](../../../../../functions/sources/kratos_parameters/Parameters/RecursivelyValidateAndAssignDefaults.md)
- [GetJSONStringFourLevels](../../../../../functions/tests/cpp_tests/sources/test_kratos_parameters/GetJSONStringFourLevels.md)
- [GetJSONStringForLevelsVariation](../../../../../functions/tests/cpp_tests/sources/test_kratos_parameters/GetJSONStringForLevelsVariation.md)
- [GetJSONStringForLevelsWrongVariation](../../../../../functions/tests/cpp_tests/sources/test_kratos_parameters/GetJSONStringForLevelsWrongVariation.md)
- [GetJSONStringForLevelsDefaults](../../../../../functions/tests/cpp_tests/sources/test_kratos_parameters/GetJSONStringForLevelsDefaults.md)
- [IsEquivalentTo](../../../../../functions/sources/kratos_parameters/Parameters/IsEquivalentTo.md)
- [HasSameKeysAndTypeOfValuesAs](../../../../../functions/sources/kratos_parameters/Parameters/HasSameKeysAndTypeOfValuesAs.md)
- [SetDouble](../../../../../functions/sources/kratos_parameters/Parameters/SetDouble.md)
- [GetJSONStringExpectedValidationOutput](../../../../../functions/tests/cpp_tests/sources/test_kratos_parameters/GetJSONStringExpectedValidationOutput.md)
- [GetJSONStringIncompleteWithExtraParameter](../../../../../functions/tests/cpp_tests/sources/test_kratos_parameters/GetJSONStringIncompleteWithExtraParameter.md)
- [AddMissingParameters](../../../../../functions/sources/kratos_parameters/Parameters/AddMissingParameters.md)
- [RecursivelyAddMissingParameters](../../../../../functions/sources/kratos_parameters/Parameters/RecursivelyAddMissingParameters.md)
- [ValidateDefaults](../../../../../functions/sources/kratos_parameters/Parameters/ValidateDefaults.md)
- [GetJSONStringIncomplete](../../../../../functions/tests/cpp_tests/sources/test_kratos_parameters/GetJSONStringIncomplete.md)
- [RecursivelyValidateDefaults](../../../../../functions/sources/kratos_parameters/Parameters/RecursivelyValidateDefaults.md)
- [AddEmptyValue](../../../../../functions/sources/kratos_parameters/Parameters/AddEmptyValue.md)
- [AddEmptyArray](../../../../../functions/sources/kratos_parameters/Parameters/AddEmptyArray.md)
- [RemoveValue](../../../../../functions/sources/kratos_parameters/Parameters/RemoveValue.md)
- [IsInt](../../../../../functions/sources/kratos_parameters/Parameters/IsInt.md)
- [IsDouble](../../../../../functions/sources/kratos_parameters/Parameters/IsDouble.md)
- [IsString](../../../../../functions/sources/kratos_parameters/Parameters/IsString.md)
- [IsStringArray](../../../../../functions/sources/kratos_parameters/Parameters/IsStringArray.md)
- [IsVector](../../../../../functions/sources/kratos_parameters/Parameters/IsVector.md)
- [IsMatrix](../../../../../functions/sources/kratos_parameters/Parameters/IsMatrix.md)
- [SetBool](../../../../../functions/sources/kratos_parameters/Parameters/SetBool.md)
- [SetVector](../../../../../functions/sources/kratos_parameters/Parameters/SetVector.md)
- [SetMatrix](../../../../../functions/sources/kratos_parameters/Parameters/SetMatrix.md)
- [AddInt](../../../../../functions/sources/kratos_parameters/Parameters/AddInt.md)
- [AddDouble](../../../../../functions/sources/kratos_parameters/Parameters/AddDouble.md)
- [AddBool](../../../../../functions/sources/kratos_parameters/Parameters/AddBool.md)
- [AddString](../../../../../functions/sources/kratos_parameters/Parameters/AddString.md)
- [GetVectorArray](../../../../../functions/sources/kratos_parameters/Parameters/GetVectorArray.md)
- [GetStringArray](../../../../../functions/sources/kratos_parameters/Parameters/GetStringArray.md)
- [SetStringArray](../../../../../functions/sources/kratos_parameters/Parameters/SetStringArray.md)
- [GetIncludedJSONString](../../../../../functions/tests/cpp_tests/sources/test_kratos_parameters/GetIncludedJSONString.md)
- [GetIncludedJSONLevel2String](../../../../../functions/tests/cpp_tests/sources/test_kratos_parameters/GetIncludedJSONLevel2String.md)
- [GetJSONStringWithIncludes](../../../../../functions/tests/cpp_tests/sources/test_kratos_parameters/GetJSONStringWithIncludes.md)
- [GetCircularIncludeJSONString](../../../../../functions/tests/cpp_tests/sources/test_kratos_parameters/GetCircularIncludeJSONString.md)
- [what](../../../../../functions/sources/exception/Exception/what.md)