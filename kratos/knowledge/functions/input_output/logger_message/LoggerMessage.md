---
type: C++ Function
title: LoggerMessage
resource: input_output/logger_message.h#L142-L143
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/logger_message/DistributedFilter/FromRoot
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`explicit LoggerMessage(std::string const& TheLabel) : mLabel(TheLabel), mLevel(1), mSeverity(Severity::INFO), mCategory(Category::STATUS), mMessageSource(), mDistributedFilter(DistributedFilter::FromRoot())`

# Calls

- [FromRoot](../../../functions/input_output/logger_message/DistributedFilter/FromRoot.md)