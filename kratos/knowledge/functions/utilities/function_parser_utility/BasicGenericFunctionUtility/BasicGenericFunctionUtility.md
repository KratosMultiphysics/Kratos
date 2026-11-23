---
type: C++ Method
title: BasicGenericFunctionUtility
resource: utilities/function_parser_utility.cpp#L51-L66
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/string_utilities/ReplaceAllSubstrings
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/function_parser_utility/BasicGenericFunctionUtility/InitializeParser
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`BasicGenericFunctionUtility::BasicGenericFunctionUtility(const std::string& rFunctionBody) : mFunctionBody(rFunctionBody)`

# Calls

- [ReplaceAllSubstrings](../../../../functions/utilities/string_utilities/ReplaceAllSubstrings.md)
- [InitializeParser](../../../../functions/utilities/function_parser_utility/BasicGenericFunctionUtility/InitializeParser.md)