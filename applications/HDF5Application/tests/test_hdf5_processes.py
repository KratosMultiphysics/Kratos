import os
from KratosMultiphysics import *
from KratosMultiphysics.HDF5Application import *
import KratosMultiphysics.KratosUnittest as KratosUnittest

from single_mesh_temporal_output_process import Factory as OutputFactory
from single_mesh_temporal_input_process import Factory as TimeInputFactory
from initialization_from_hdf5_process import Factory as InitializationFactory

class TestHDF5Processes(KratosUnittest.TestCase):

    def setUp(self):
        self.model = Model()
        self.buffer_size = 2

    def tearDown(self):
        self._remove_h5_files("WriteModelPart")

    def testSingleMeshTemporalOutputInput(self):
        """
        Output ModelPart using SingleMeshTemporalOutputProcess.
        Use output file to set values in a new ModelPart with SingleMeshTemporalInputProcess.
        """
        write_model_part = self._CreateNewModelPart("WriteModelPartForTemporalInput")
        self._InitializeModelPart(write_model_part)

        read_model_part = self._CreateNewModelPart("ReadModelPart")

        output_settings = Parameters(r"""{
            "Parameters" : {
                "model_part_name" : "WriteModelPartForTemporalInput",
                "nodal_solution_step_data_settings" : {
                    "list_of_variables": ["VELOCITY", "DENSITY"]
                },
                "nodal_data_value_settings" : {
                    "list_of_variables": ["PRESSURE"]
                },
                "element_data_value_settings" : {
                    "list_of_variables": ["TEMPERATURE"]
                },
                "output_time_settings" : {
                    "output_step_frequency": 1
                }
            }
        }""")
        output_process = OutputFactory(output_settings, self.model)

        input_settings = Parameters(r"""{
            "Parameters" : {
                "model_part_name" : "ReadModelPart",
                "nodal_solution_step_data_settings" : {
                    "list_of_variables": ["VELOCITY", "DENSITY"]
                },
                "nodal_data_value_settings" : {
                    "list_of_variables": ["PRESSURE"]
                },
                "element_data_value_settings" : {
                    "list_of_variables": ["TEMPERATURE"]
                },
                "file_name": "WriteModelPartForTemporalInput"
            }
        }""")
        input_process = TimeInputFactory(input_settings, self.model)

        output_process.ExecuteBeforeSolutionLoop()

        for step in range(1,3):
            # Run step on write model part and dump results
            self._SimulateTimeStep(write_model_part, 0.1*step)
            output_process.ExecuteFinalizeSolutionStep()

            # Read results on read model part
            read_model_part.CloneTimeStep(0.1*step)
            input_process.ExecuteInitializeSolutionStep()

            for read_node,write_node in zip(read_model_part.Nodes, write_model_part.Nodes):
                self.assertEqual(read_node.GetSolutionStepValue(VELOCITY_X,0), write_node.GetSolutionStepValue(VELOCITY_X,0))
                self.assertEqual(read_node.GetSolutionStepValue(VELOCITY_Y,0), write_node.GetSolutionStepValue(VELOCITY_Y,0))
                self.assertEqual(read_node.GetSolutionStepValue(VELOCITY_Z,0), write_node.GetSolutionStepValue(VELOCITY_Z,0))
                self.assertEqual(read_node.GetSolutionStepValue(DENSITY,0),    write_node.GetSolutionStepValue(DENSITY,0))
                # reaction is not written, so it should not be updated
                self.assertEqual(read_node.GetSolutionStepValue(REACTION_X,0), 0.0)
                self.assertEqual(read_node.GetSolutionStepValue(REACTION_Y,0), 0.0)
                self.assertEqual(read_node.GetSolutionStepValue(REACTION_Z,0), 0.0)

                self.assertEqual(read_node.GetValue(PRESSURE), write_node.GetValue(PRESSURE))

            for read_element,write_element in zip(read_model_part.Elements, write_model_part.Elements):
                self.assertEqual(read_element.GetValue(TEMPERATURE), write_element.GetValue(TEMPERATURE))

    def testSingleMeshTemporalOutputInitialization(self):
        """
        Output ModelPart using SingleMeshTemporalOutputProcess.
        Use output file to initialize (restart) a ModelPart with initialization_from_hdf5_process.
        """
        write_model_part = self._CreateNewModelPart("WriteModelPartForInitialization")
        self._InitializeModelPart(write_model_part)

        read_model_part = self._CreateNewModelPart("ReadModelPart")

        output_settings = Parameters(r"""{
            "Parameters" : {
                "model_part_name" : "WriteModelPartForInitialization",
                "nodal_solution_step_data_settings" : {
                    "list_of_variables": ["VELOCITY", "DENSITY"]
                },
                "nodal_data_value_settings" : {
                    "list_of_variables": ["PRESSURE"]
                },
                "element_data_value_settings" : {
                    "list_of_variables": ["TEMPERATURE"]
                },
                "output_time_settings" : {
                    "output_step_frequency": 1
                }
            }
        }""")
        output_process = OutputFactory(output_settings, self.model)

        input_settings = Parameters(r"""{
            "Parameters" : {
                "model_part_name" : "ReadModelPart",
                "nodal_solution_step_data_settings" : {
                    "list_of_variables": ["VELOCITY", "DENSITY"]
                },
                "nodal_data_value_settings" : {
                    "list_of_variables": ["PRESSURE"]
                },
                "element_data_value_settings" : {
                    "list_of_variables": ["TEMPERATURE"]
                },
                "file_name": "WriteModelPartForInitialization-0.2000"
            }
        }""")
        input_process = InitializationFactory(input_settings, self.model)

        output_process.ExecuteBeforeSolutionLoop()

        for step in range(1,3):
            # Run step on write model part and dump results
            self._SimulateTimeStep(write_model_part, 0.1*step)
            output_process.ExecuteFinalizeSolutionStep()

        # Now initialize the read model part and compare the two model parts
        input_process.ExecuteInitialize()

        for read_node,write_node in zip(read_model_part.Nodes, write_model_part.Nodes):
            self.assertEqual(read_node.GetSolutionStepValue(VELOCITY_X,0), write_node.GetSolutionStepValue(VELOCITY_X,0))
            self.assertEqual(read_node.GetSolutionStepValue(VELOCITY_Y,0), write_node.GetSolutionStepValue(VELOCITY_Y,0))
            self.assertEqual(read_node.GetSolutionStepValue(VELOCITY_Z,0), write_node.GetSolutionStepValue(VELOCITY_Z,0))
            self.assertEqual(read_node.GetSolutionStepValue(DENSITY,0),    write_node.GetSolutionStepValue(DENSITY,0))
            # reaction is not written, so it should not be updated
            self.assertEqual(read_node.GetSolutionStepValue(REACTION_X,0), 0.0)
            self.assertEqual(read_node.GetSolutionStepValue(REACTION_Y,0), 0.0)
            self.assertEqual(read_node.GetSolutionStepValue(REACTION_Z,0), 0.0)

            self.assertEqual(read_node.GetValue(PRESSURE), write_node.GetValue(PRESSURE))

        for read_element,write_element in zip(read_model_part.Elements, write_model_part.Elements):
            self.assertEqual(read_element.GetValue(TEMPERATURE), write_element.GetValue(TEMPERATURE))

    def _CreateNewModelPart(self,label):
        model_part = self.model.CreateModelPart(label,self.buffer_size)
        model_part.AddNodalSolutionStepVariable(VELOCITY)
        model_part.AddNodalSolutionStepVariable(REACTION)
        model_part.AddNodalSolutionStepVariable(DENSITY)

        model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        model_part.CreateNewNode(3, 1.0, 1.0, 0.0)
        model_part.CreateNewNode(4, 0.0, 0.0, 1.0)

        properties = model_part.GetProperties(1,0)

        model_part.CreateNewElement("Element3D4N", 1, [1, 2, 3, 4], properties)
        model_part.CreateNewCondition("SurfaceCondition3D3N", 1, [2, 3, 4], properties)

        return model_part

    def _InitializeModelPart(self, model_part):

        for node in model_part.Nodes:
            node.SetSolutionStepValue(VELOCITY_X, 0, 10. + node.X)
            node.SetSolutionStepValue(VELOCITY_Y, 0, 10. + node.Y)
            node.SetSolutionStepValue(VELOCITY_Z, 0, 10. + node.Z)
            node.SetSolutionStepValue(REACTION_X, 0, 20. + node.X)
            node.SetSolutionStepValue(REACTION_Y, 0, 20. + node.Y)
            node.SetSolutionStepValue(REACTION_Z, 0, 20. + node.Z)
            node.SetSolutionStepValue(DENSITY, 0, 5.*node.Id)

            node.SetValue(PRESSURE, 50. - node.Id)

        for element in model_part.Elements:
            element.SetValue(TEMPERATURE,3. + element.Id)



    def _SimulateTimeStep(self, model_part, new_time):
        model_part.CloneTimeStep(new_time)
        update_factor = 1+new_time
        for node in model_part.Nodes:
            u = node.GetSolutionStepValue(VELOCITY,0)
            f = node.GetSolutionStepValue(REACTION,0)
            rho = node.GetSolutionStepValue(DENSITY,0)

            node.SetSolutionStepValue(VELOCITY, 0, u*update_factor)
            node.SetSolutionStepValue(REACTION, 0, f*update_factor)
            node.SetSolutionStepValue(DENSITY, 0, rho*update_factor)

            p = node.GetValue(PRESSURE)
            node.SetValue(PRESSURE, p*update_factor)

        for element in model_part.Elements:
            t = element.GetValue(TEMPERATURE)
            element.SetValue(TEMPERATURE,t*update_factor)


    def _remove_h5_files(self, model_part_name):
        for name in os.listdir():
            if name.find(model_part_name) == 0:
                kratos_utilities.DeleteFileIfExisting(name)

if __name__ == "__main__":
    KratosUnittest.main()