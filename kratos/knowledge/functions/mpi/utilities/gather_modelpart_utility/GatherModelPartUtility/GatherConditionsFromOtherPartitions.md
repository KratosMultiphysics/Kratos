---
type: C++ Method
title: GatherConditionsFromOtherPartitions
resource: mpi/utilities/gather_modelpart_utility.cpp#L274-L298
visibility: private
generated:
  by: okf-rs/0.7.0
---

# Signature

`void GatherModelPartUtility::GatherConditionsFromOtherPartitions( ModelPart& rModelPart, const std::map<int, std::vector<std::size_t>>& rConditionsToBring, const bool CallExecuteAfterBringingEntities, const int EchoLevel )`