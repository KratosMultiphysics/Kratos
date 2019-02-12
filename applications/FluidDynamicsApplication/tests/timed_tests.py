# using an idea from https://hakibenita.com/timing-tests-in-python-for-fun-and-profit

import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication
import test_FluidDynamicsApplication as fluid_tests
import KratosMultiphysics.KratosUnittest as KratosUnittest

import run_cpp_unit_tests
import subprocess

import time
from operator import itemgetter

class TimeLoggingTestResult(KratosUnittest.TextTestResult):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._test_timings = []
        self._total_elapsed = 0.0

    def startTest(self, test):
        self._start_time = time.time()
        self.stream.write("\nRunning {}\n".format(self.getDescription(test)))
        super().startTest(test)

    def addSuccess(self, test):
        elapsed = time.time() - self._start_time
        name = self.getDescription(test)
        self._test_timings.append((name,elapsed))
        self._total_elapsed += elapsed
        self.stream.write("\n{} ({:.03}s)\n".format(name, elapsed))
        super().addSuccess(test)

    def getTestTimings(self):
        return self._test_timings

    def getTestTotalElapsed(self):
        return self._total_elapsed

class TimeLoggingTextTestRunner(KratosUnittest.TextTestRunner):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def run(self, tests):
        result = super().run(tests)

        timings = result.getTestTimings()
        total_elapsed = result.getTestTotalElapsed()
        for name,elapsed in sorted(timings, key=itemgetter(1)):
            self.stream.writeln("{} ({:.03f}s, {:.02f}%)".format(name, elapsed, 100.*(elapsed/total_elapsed)))
        self.stream.writeln("Total ({:.03f}s, {:.02f}%)".format(total_elapsed, 100.))

        return result

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosMultiphysics.Logger.PrintInfo("Unittests", "\nRunning cpp unit tests ...")
    run_cpp_unit_tests.run()
    KratosMultiphysics.Logger.PrintInfo("Unittests", "Finished running cpp unit tests!")

    KratosMultiphysics.Logger.PrintInfo("Unittests", "\nRunning mpi python tests ...")
    try:
        import KratosMultiphysics.mpi as KratosMPI
        import KratosMultiphysics.MetisApplication as MetisApplication
        import KratosMultiphysics.TrilinosApplication as TrilinosApplication
        p = subprocess.Popen(["mpiexec", "-np", "2", "python3", "test_FluidDynamicsApplication_mpi.py"], stdout=subprocess.PIPE)
        p.wait()
        KratosMultiphysics.Logger.PrintInfo("Unittests", "Finished mpi python tests!")
    except ImportError:
        KratosMultiphysics.Logger.PrintInfo("Unittests", "mpi is not available!")

    KratosMultiphysics.Logger.PrintInfo("Unittests", "\nRunning python tests ...")
    #KratosUnittest.runTests(AssembleTestSuites())
    #KratosUnittest.TextTestRunner(verbosity=0, buffer=True).run(AssembleTestSuites()["all"])
    suites = fluid_tests.AssembleTestSuites()
    TimeLoggingTextTestRunner(verbosity=1, buffer=True, resultclass=TimeLoggingTestResult).run(suites["all"])
    #TimeLoggingTextTestRunner(verbosity=1, buffer=True, resultclass=TimeLoggingTestResult).run(suites["validation"])
    KratosMultiphysics.Logger.PrintInfo("Unittests", "Finished python tests!")