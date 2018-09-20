from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.KratosUnittest as KratosUnittest

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestObjectPrinting(KratosUnittest.TestCase):
    def test_Properties_printing(self):
        pass
    def test_ModelPart_printing(self):
        model_part = KratosMultiphysics.ModelPart("Main")
        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("test_model_part_io_read"))
        model_part_io.ReadModelPart(model_part)

        print(model_part)




if __name__ == '__main__':
    KratosUnittest.main()

