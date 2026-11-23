---
type: C++ Method
title: CreateDistanceGetterFunctor
resource: processes/apply_ray_casting_process.cpp#L513-L525
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/apply_ray_casting_interface_recognition_process/ApplyRayCastingInterfaceRecognitionProcess<TDim>/ApplyRayCastingInterfaceRecognitionProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/apply_ray_casting_process/ApplyRayCastingProcess<TDim>/ApplyRayCastingProcess
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::function<double&(Node& rNode, const Variable<double>& rDistanceVariable)> ApplyRayCastingProcess<TDim>::CreateDistanceGetterFunctor() const`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)

# Called by

- [ApplyRayCastingInterfaceRecognitionProcess](../../../../functions/processes/apply_ray_casting_interface_recognition_process/ApplyRayCastingInterfaceRecognitionProcess<TDim>/ApplyRayCastingInterfaceRecognitionProcess.md)
- [ApplyRayCastingProcess](../../../../functions/processes/apply_ray_casting_process/ApplyRayCastingProcess<TDim>/ApplyRayCastingProcess.md)