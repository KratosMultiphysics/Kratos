from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.CoSimulationApplication.coupling_interface_data import CouplingInterfaceData
from KratosMultiphysics.CoSimulationApplication.coupling_interface_residual import CouplingInterfaceResidual

import numpy as np

class TestCouplingInterfaceResidual(KratosUnittest.TestCase):
    def setUp(self):
        self.model = KM.Model()
        self.model_part = self.model.CreateModelPart("default")
        self.model_part.AddNodalSolutionStepVariable(KM.PRESSURE)
        self.model_part.AddNodalSolutionStepVariable(KM.VELOCITY)
        self.dimension = 3
        self.model_part.ProcessInfo[KM.DOMAIN_SIZE] = self.dimension

        self.num_nodes = 11

        for i in range(self.num_nodes):
            self.model_part.CreateNewNode(i+1, i*0.1, 0.0, 0.0)

        data_settings_scalar = KM.Parameters("""{
            "model_part_name" : "default",
            "variable_name"   : "PRESSURE"
        }""")
        self.interface_data_scalar = CouplingInterfaceData(data_settings_scalar, self.model)
        self.interface_data_scalar.Initialize()

        data_settings_vector = KM.Parameters("""{
            "model_part_name" : "default",
            "variable_name"   : "VELOCITY",
            "dimension"       : %d
        }""" % self.dimension )
        self.interface_data_vector = CouplingInterfaceData(data_settings_vector, self.model)
        self.interface_data_vector.Initialize()

    def test_interface_residual(self):
        # setting initial values
        init_data_scalar, init_data_vector = self.__SetData(1)

        interface_res_scalar = CouplingInterfaceResidual(self.interface_data_scalar)
        interface_res_vector = CouplingInterfaceResidual(self.interface_data_vector)

        # checking just to make sure the IntefaceData gives the correct results
        self.__CheckData(init_data_scalar, self.interface_data_scalar.GetData())
        self.__CheckData(init_data_vector, self.interface_data_vector.GetData())

        # updating values to compute residual
        updated_data_scalar, updated_data_vector = self.__SetData(3)

        # checking just to make sure the IntefaceData gives the correct results
        self.__CheckData(updated_data_scalar, self.interface_data_scalar.GetData())
        self.__CheckData(updated_data_vector, self.interface_data_vector.GetData())

        scalar_res = interface_res_scalar.GetResidual()
        vector_res = interface_res_vector.GetResidual()

        self.assertEqual(len(scalar_res), self.num_nodes)
        self.assertEqual(len(vector_res), self.num_nodes*self.dimension)

        exp_res_scalar = updated_data_scalar - init_data_scalar
        exp_res_vector = updated_data_vector - init_data_vector

        self.__CheckData(exp_res_scalar, scalar_res)
        self.__CheckData(exp_res_vector, vector_res)

        self.assertAlmostEqual(7.9280616416775764, interface_res_scalar.GetNorm(), 12)
        self.assertAlmostEqual(352.02261291002316, interface_res_vector.GetNorm(), 12)

        # Updating data
        scalar_res = interface_res_scalar.UpdateReferenceData()
        vector_res = interface_res_vector.UpdateReferenceData()

        # updating values to compute residual
        updated_data_scalar_2, updated_data_vector_2 = self.__SetData(15)

        scalar_res_2 = interface_res_scalar.GetResidual()
        vector_res_2 = interface_res_vector.GetResidual()

        exp_res_scalar_2 = updated_data_scalar_2 - updated_data_scalar
        exp_res_vector_2 = updated_data_vector_2 - updated_data_vector

        self.__CheckData(exp_res_scalar_2, scalar_res_2)
        self.__CheckData(exp_res_vector_2, vector_res_2)

        self.assertAlmostEqual(66.778072196306383, interface_res_scalar.GetNorm(), 12)
        self.assertAlmostEqual(2112.1356774601381, interface_res_vector.GetNorm(), 12)


    def __SetData(self, offset_factor):
        for node in self.model_part.Nodes:
            idx = node.Id + offset_factor - 1 # -1 bcs the node-ids already have an offset of 1
            node.SetSolutionStepValue(KM.PRESSURE, 0, NodeScalarValue(idx))
            node.SetSolutionStepValue(KM.VELOCITY, 0, NodeVectorValue(idx))

        data_scalar = np.array([NodeScalarValue(i+offset_factor) for i in range(self.num_nodes)])
        data_vector = np.array([v for i in range(self.num_nodes) for v in NodeVectorValue(i+offset_factor)]) # double list comprehension to get the values in a consecutive array

        return data_scalar, data_vector

    def __CheckData(self, exp_data, data):
        self.assertEqual(len(exp_data), len(data))

        for exp_val, val in zip(exp_data, data):
            self.assertAlmostEqual(exp_val, val)

def NodeScalarValue(the_id):
    return the_id**1.5

def NodeVectorValue(the_id):
    return [the_id*14.7, the_id*19.2-10.5, the_id*303.9]

if __name__ == '__main__':
    KratosUnittest.main()
