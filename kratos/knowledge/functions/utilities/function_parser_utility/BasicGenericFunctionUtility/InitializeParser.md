---
type: C++ Method
title: InitializeParser
resource: utilities/function_parser_utility.cpp#L156-L250
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/string_utilities/ContainsPartialString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/string_utilities/ReplaceAllSubstrings
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/string_utilities/SplitStringByDelimiter
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/function_parser_utility/BasicGenericFunctionUtility/BasicGenericFunctionUtility
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/function_parser_utility/GenericFunctionUtility/GenericFunctionUtility
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void BasicGenericFunctionUtility::InitializeParser()`

# Calls

- [ContainsPartialString](../../../../functions/utilities/string_utilities/ContainsPartialString.md)
- [ReplaceAllSubstrings](../../../../functions/utilities/string_utilities/ReplaceAllSubstrings.md)
- [SplitStringByDelimiter](../../../../functions/utilities/string_utilities/SplitStringByDelimiter.md)

# Called by

- [BasicGenericFunctionUtility](../../../../functions/utilities/function_parser_utility/BasicGenericFunctionUtility/BasicGenericFunctionUtility.md)
- [GenericFunctionUtility](../../../../functions/utilities/function_parser_utility/GenericFunctionUtility/GenericFunctionUtility.md)