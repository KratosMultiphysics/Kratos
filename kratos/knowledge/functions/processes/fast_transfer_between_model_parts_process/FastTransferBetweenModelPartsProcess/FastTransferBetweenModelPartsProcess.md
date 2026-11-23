---
type: C++ Method
title: FastTransferBetweenModelPartsProcess
resource: processes/fast_transfer_between_model_parts_process.cpp#L25-L44
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/input_output/vtk_output/VtkOutput/WriteModelPartWithoutNodesToFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/skin_detection_process/SkinDetectionProcess/ExecuteInitialize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/skin_detection_process/SkinDetectionProcess/ExecuteInitializeSolutionStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/processes/test_fast_transfer_between_model_parts_process/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/model_part_combination_utilities/ModelPartCombinationUtilities/RecursiveAddEntities
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`FastTransferBetweenModelPartsProcess::FastTransferBetweenModelPartsProcess( ModelPart& rDestinationModelPart, ModelPart& rOriginModelPart, const EntityTransfered Entity, const Flags Flag, const bool ReplicateEntities ) : Process(), mrDestinationModelPart(rDestinationModelPart), mrOriginModelPart(rOriginModelPart), mEntity(Entity), mFlag(Flag)`

# Called by

- [WriteModelPartWithoutNodesToFile](../../../../functions/input_output/vtk_output/VtkOutput/WriteModelPartWithoutNodesToFile.md)
- [ExecuteInitialize](../../../../functions/python_scripts/skin_detection_process/SkinDetectionProcess/ExecuteInitialize.md)
- [ExecuteInitializeSolutionStep](../../../../functions/python_scripts/skin_detection_process/SkinDetectionProcess/ExecuteInitializeSolutionStep.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/processes/test_fast_transfer_between_model_parts_process/KRATOS_TEST_CASE_IN_SUITE.md)
- [RecursiveAddEntities](../../../../functions/utilities/model_part_combination_utilities/ModelPartCombinationUtilities/RecursiveAddEntities.md)