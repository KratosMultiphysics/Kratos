from KratosMultiphysics import *
Tester.SetVerbosity(Tester.Verbosity.TESTS_OUTPUTS)
Tester.RunTestCases("TestGraphConstruction")
Tester.RunTestCases("TestGraphContiguousRowConstruction")
#Tester.RunTestCases("TestDistributedGraphConstruction")
Tester.RunTestCases("TestCsrMatrixAssemble")
Tester.RunTestCases("TestSystemVectorAssembly")

Tester.RunTestCases("TestCSRConstruction")

Tester.RunTestCases("TestSpMV")


#Tester.RunTestCases("TestPerformanceBenchmarkSparseGraph")
#Tester.RunTestCases("TestPerformanceBenchmarkSparseContiguousRowGraph")
#Tester.RunTestSuite("KratosCoreFastSuite")
