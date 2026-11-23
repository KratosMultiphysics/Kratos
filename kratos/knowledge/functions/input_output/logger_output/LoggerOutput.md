---
type: C++ Function
title: LoggerOutput
resource: input_output/logger_output.h#L70-L81
generated:
  by: okf-rs/0.7.0
---

# Signature

`explicit LoggerOutput(std::ostream& rOutputStream) : mpStream(&rOutputStream), mMaxLevel(1), mSeverity(LoggerMessage::Severity::INFO), mCategory(LoggerMessage::Category::STATUS)`