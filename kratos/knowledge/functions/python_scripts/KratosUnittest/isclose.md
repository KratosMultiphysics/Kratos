---
type: Python Function
title: isclose
resource: python_scripts/KratosUnittest.py#L274-L278
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/mpi/tests/test_gather_modelpart_utility/TestGatherModelPartUtility/testGatherOnMaster
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_gather_modelpart_utility/TestGatherModelPartUtility/testScatterFromMaster
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/KratosUnittest/TestCase/assertIsClose
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__CompareEnsightSolutionFile
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):`

# Called by

- [testGatherOnMaster](../../../functions/mpi/tests/test_gather_modelpart_utility/TestGatherModelPartUtility/testGatherOnMaster.md)
- [testScatterFromMaster](../../../functions/mpi/tests/test_gather_modelpart_utility/TestGatherModelPartUtility/testScatterFromMaster.md)
- [assertIsClose](../../../functions/python_scripts/KratosUnittest/TestCase/assertIsClose.md)
- [__CompareEnsightSolutionFile](../../../functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__CompareEnsightSolutionFile.md)