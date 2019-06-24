from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.CoSimulationApplication.coupling_interface_data import CouplingInterfaceData

# The expected definitions are here to make the handling of the
# multiline-stings easier (no need to deal with indentation)
coupling_interface_data_str = '''CouplingInterfaceData:
	ModelPart: "mp_4_test"
	IsDistributed: False
	Variable: "DISPLACEMENT" (Vector with dimension: 2)
	Location: "node_historical"
'''

class TestCouplingInterfaceData(KratosUnittest.TestCase):

    def test_str(self):
        settings = KM.Parameters("""{
            "model_part_name" : "mp_4_test",
            "variable_name"   : "DISPLACEMENT",
            "dimension"       : 2
        }""")

        model = KM.Model()
        model.CreateModelPart("mp_4_test")

        coupling_data = CouplingInterfaceData(settings, model)

        self.assertMultiLineEqual(str(coupling_data), coupling_interface_data_str)




if __name__ == '__main__':
    KratosUnittest.main()
