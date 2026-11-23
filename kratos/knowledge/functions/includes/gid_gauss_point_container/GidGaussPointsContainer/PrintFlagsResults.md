---
type: C++ Method
title: PrintFlagsResults
resource: includes/gid_gauss_point_container.h#L93-L122
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/includes/gid_io/PrintFlagsOnGaussPoints
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`virtual void PrintFlagsResults( GiD_FILE ResultFile, const Kratos::Flags& rFlag, const std::string& rFlagName, const ModelPart& rModelPart, const double SolutionTag )`

# Called by

- [PrintFlagsOnGaussPoints](../../../../functions/includes/gid_io/PrintFlagsOnGaussPoints.md)