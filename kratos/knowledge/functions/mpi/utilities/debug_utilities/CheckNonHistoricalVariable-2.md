---
type: C++ Function
title: CheckNonHistoricalVariable
resource: mpi/utilities/debug_utilities.h#L117-L155
generated:
  by: okf-rs/0.7.0
---

# Signature

`static void CheckNonHistoricalVariable( ModelPart & rModelPart, const TContainerType & rContainer, const TVarType & rVariable, GlobalPointerCommunicator<typename TContainerType::data_type> & rPointerCommunicator, std::unordered_map<int, GlobalPointer<typename TContainerType::data_type>> & gp_map)`