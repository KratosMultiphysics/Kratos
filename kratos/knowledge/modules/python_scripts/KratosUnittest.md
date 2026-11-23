---
type: Python Module
title: KratosUnittest
resource: python_scripts/KratosUnittest.py#L1-L310
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/kratosmultiphysics-as-km
    resolved_by: tree-sitter
    confidence: exact
  - target: external/kratosmultiphysics
    resolved_by: tree-sitter
    confidence: exact
  - target: external/kratosmultiphysics-kratos-utilities
    resolved_by: tree-sitter
    confidence: exact
  - target: external/unittest
    resolved_by: tree-sitter
    confidence: exact
  - target: external/unittest-util
    resolved_by: tree-sitter
    confidence: exact
  - target: external/contextlib
    resolved_by: tree-sitter
    confidence: exact
  - target: external/argparse
    resolved_by: tree-sitter
    confidence: exact
  - target: external/sys
    resolved_by: tree-sitter
    confidence: exact
  - target: external/os
    resolved_by: tree-sitter
    confidence: exact
  - target: external/time
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [TestLoader](../../classes/python_scripts/KratosUnittest/TestLoader.md)
- [loadTestsFromTestCases](../../functions/python_scripts/KratosUnittest/TestLoader/loadTestsFromTestCases.md)
- [TestCase](../../classes/python_scripts/KratosUnittest/TestCase.md)
- [run](../../functions/python_scripts/KratosUnittest/TestCase/run.md)
- [skipTestIfApplicationsNotAvailable](../../functions/python_scripts/KratosUnittest/TestCase/skipTestIfApplicationsNotAvailable.md)
- [assertEqualTolerance](../../functions/python_scripts/KratosUnittest/TestCase/assertEqualTolerance.md)
- [assertIsClose](../../functions/python_scripts/KratosUnittest/TestCase/assertIsClose.md)
- [assertVectorAlmostEqual](../../functions/python_scripts/KratosUnittest/TestCase/assertVectorAlmostEqual.md)
- [LazyErrMsg](../../classes/python_scripts/KratosUnittest/LazyErrMsg.md)
- [__init__](../../functions/python_scripts/KratosUnittest/LazyErrMsg/__init__.md)
- [__str__](../../functions/python_scripts/KratosUnittest/LazyErrMsg/__str__.md)
- [assertMatrixAlmostEqual](../../functions/python_scripts/KratosUnittest/TestCase/assertMatrixAlmostEqual.md)
- [LazyDimErrMsg](../../classes/python_scripts/KratosUnittest/LazyDimErrMsg.md)
- [__init__](../../functions/python_scripts/KratosUnittest/LazyDimErrMsg/__init__.md)
- [__str__](../../functions/python_scripts/KratosUnittest/LazyDimErrMsg/__str__.md)
- [LazyValErrMsg](../../classes/python_scripts/KratosUnittest/LazyValErrMsg.md)
- [__init__](../../functions/python_scripts/KratosUnittest/LazyValErrMsg/__init__.md)
- [__str__](../../functions/python_scripts/KratosUnittest/LazyValErrMsg/__str__.md)
- [KratosTextTestResult](../../classes/python_scripts/KratosUnittest/KratosTextTestResult.md)
- [__init__](../../functions/python_scripts/KratosUnittest/KratosTextTestResult/__init__.md)
- [startTest](../../functions/python_scripts/KratosUnittest/KratosTextTestResult/startTest.md)
- [addError](../../functions/python_scripts/KratosUnittest/KratosTextTestResult/addError.md)
- [addExpectedFailure](../../functions/python_scripts/KratosUnittest/KratosTextTestResult/addExpectedFailure.md)
- [addFailure](../../functions/python_scripts/KratosUnittest/KratosTextTestResult/addFailure.md)
- [addSkip](../../functions/python_scripts/KratosUnittest/KratosTextTestResult/addSkip.md)
- [addSuccess](../../functions/python_scripts/KratosUnittest/KratosTextTestResult/addSuccess.md)
- [addUnexpectedSuccess](../../functions/python_scripts/KratosUnittest/KratosTextTestResult/addUnexpectedSuccess.md)
- [getDescription](../../functions/python_scripts/KratosUnittest/KratosTextTestResult/getDescription.md)
- [skipIfApplicationsNotAvailable](../../functions/python_scripts/KratosUnittest/skipIfApplicationsNotAvailable.md)
- [SupressConsoleOutput](../../functions/python_scripts/KratosUnittest/SupressConsoleOutput.md)
- [SupressConsoleError](../../functions/python_scripts/KratosUnittest/SupressConsoleError.md)
- [SupressAllConsole](../../functions/python_scripts/KratosUnittest/SupressAllConsole.md)
- [main](../../functions/python_scripts/KratosUnittest/main.md)
- [runTests](../../functions/python_scripts/KratosUnittest/runTests.md)
- [isclose](../../functions/python_scripts/KratosUnittest/isclose.md)
- [WorkFolderScope](../../classes/python_scripts/KratosUnittest/WorkFolderScope.md)
- [__init__](../../functions/python_scripts/KratosUnittest/WorkFolderScope/__init__.md)
- [__enter__](../../functions/python_scripts/KratosUnittest/WorkFolderScope/__enter__.md)
- [__exit__](../../functions/python_scripts/KratosUnittest/WorkFolderScope/__exit__.md)

# Imports

- `KratosMultiphysics as KM`
- `KratosMultiphysics`
- `KratosMultiphysics.kratos_utilities`
- `unittest`
- `unittest.util`
- `contextlib`
- `argparse`
- `sys`
- `os`
- `time`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)