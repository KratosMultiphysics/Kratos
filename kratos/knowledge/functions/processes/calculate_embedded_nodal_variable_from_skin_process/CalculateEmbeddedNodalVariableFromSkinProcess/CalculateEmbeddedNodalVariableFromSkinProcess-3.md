---
type: C++ Method
title: CalculateEmbeddedNodalVariableFromSkinProcess
resource: processes/calculate_embedded_nodal_variable_from_skin_process.h#L661-L676
visibility: private
generated:
  by: okf-rs/0.7.0
---

# Signature

`CalculateEmbeddedNodalVariableFromSkinProcess( ModelPart &rBaseModelPart, ModelPart &rSkinModelPart, Parameters rSettings) : CalculateEmbeddedNodalVariableFromSkinProcess( rBaseModelPart, rSkinModelPart, rSettings["linear_solver_settings"], KratosComponents<Variable<TVarType>>::Get(rSettings["skin_variable_name"].GetString()), KratosComponents<Variable<TVarType>>::Get(rSettings["embedded_nodal_variable_name"].GetString()), rSettings["gradient_penalty_coefficient"].GetDouble(), rSettings["buffer_position"].GetInt(), rSettings["aux_model_part_name"].GetString(), rSettings["echo_level"].GetInt())`