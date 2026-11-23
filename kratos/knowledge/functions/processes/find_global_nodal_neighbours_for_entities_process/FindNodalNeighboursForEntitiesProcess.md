---
type: C++ Function
title: FindNodalNeighboursForEntitiesProcess
resource: processes/find_global_nodal_neighbours_for_entities_process.h#L57-L64
generated:
  by: okf-rs/0.7.0
---

# Signature

`FindNodalNeighboursForEntitiesProcess( ModelPart& rModelPart, const Variable<GlobalPointersVector<NodeType>>& rOutputVariable) : mrModelPart(rModelPart), mrDataCommunicator(rModelPart.GetCommunicator().GetDataCommunicator()), mrOutputVariable(rOutputVariable)`