---
type: C++ Method
title: WriteData
resource: input_output/vtu_output.cpp#L1368-L1387
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/vtu_output/VtuOutput/WriteIntegrationPointData
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void VtuOutput::WriteData( std::vector<std::pair<std::string, std::string>>& rBlockFileNameInfo, TXmlElementDataWrapperCreateFunctor&& rElementDataWrapperCreateFunctor, TXmlElementDataWrapperAppendFunctor&& rElementDataWrapperAppendFunctor, UnstructuredGridData& rUnstructuredGridData, const std::string& rOutputPrefix, const IndexType Step) const`

# Calls

- [WriteIntegrationPointData](../../../../functions/input_output/vtu_output/VtuOutput/WriteIntegrationPointData.md)