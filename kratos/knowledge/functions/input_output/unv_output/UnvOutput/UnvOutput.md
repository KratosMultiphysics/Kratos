---
type: C++ Method
title: UnvOutput
resource: input_output/unv_output.cpp#L28-L31
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/python_scripts/unv_output_process/UnvOutputProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`UnvOutput::UnvOutput(Kratos::ModelPart &modelPart, const std::string &outFileWithoutExtension) : mrOutputModelPart(modelPart), mOutputFileName(outFileWithoutExtension + ".unv")`

# Called by

- [__init__](../../../../functions/python_scripts/unv_output_process/UnvOutputProcess/__init__.md)