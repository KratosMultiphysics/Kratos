---
type: C++ Function
title: WriteNodalFlags
resource: includes/gid_io.h#L847-L865
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/timer/Timer/Start
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/timer/Timer/Stop
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python/add_io_to_python/AddIOToPython
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/gid_output_process/GiDOutputProcess/__write_nodal_flags
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void WriteNodalFlags( const Kratos::Flags& rFlag, const std::string& rFlagName, const NodesContainerType& rNodes, const double SolutionTag )`

# Calls

- [Start](../../../functions/utilities/timer/Timer/Start.md)
- [Stop](../../../functions/utilities/timer/Timer/Stop.md)

# Called by

- [AddIOToPython](../../../functions/python/add_io_to_python/AddIOToPython.md)
- [__write_nodal_flags](../../../functions/python_scripts/gid_output_process/GiDOutputProcess/__write_nodal_flags.md)