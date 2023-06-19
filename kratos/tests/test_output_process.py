import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestOutputProcess(KratosUnittest.TestCase):

    def testOutputProcessSubclass(self):
        output_process = KM.OutputProcess()
        self.assertTrue(issubclass(type(output_process), KM.OutputProcess))
        self.assertTrue(issubclass(type(output_process), KM.Process))

    def testPythonOutputProcessSubclass(self):
        class PythonOutputProcess(KM.OutputProcess):
            pass
        output_process = PythonOutputProcess()
        self.assertTrue(issubclass(type(output_process), KM.OutputProcess))
        self.assertTrue(issubclass(type(output_process), KM.Process))
    
    def testDerivedPythonOutputProcessSubclass(self):
        class PythonOutputProcess(KM.OutputProcess):
            pass
        class DerivedPythonOutputProcess(PythonOutputProcess):
            pass
        output_process = DerivedPythonOutputProcess()
        self.assertTrue(issubclass(type(output_process), KM.OutputProcess))
        self.assertTrue(issubclass(type(output_process), KM.Process))

if __name__ == "__main__":
    KratosUnittest.main()
