---
type: C++ Method
title: GenericFunctionUtility
resource: utilities/function_parser_utility.cpp#L255-L273
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetDouble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/function_parser_utility/BasicGenericFunctionUtility/InitializeParser
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`GenericFunctionUtility::GenericFunctionUtility( const std::string& rFunctionBody, Parameters LocalSystem ) : BasicGenericFunctionUtility(rFunctionBody)`

# Calls

- [GetDouble](../../../../functions/sources/kratos_parameters/Parameters/GetDouble.md)
- [InitializeParser](../../../../functions/utilities/function_parser_utility/BasicGenericFunctionUtility/InitializeParser.md)