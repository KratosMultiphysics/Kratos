import os
import KratosMultiphysics
from KratosMultiphysics import *
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utilities

from KratosMultiphysics.HDF5Application.single_mesh_temporal_output_process import Factory as OutputFactory
from KratosMultiphysics.HDF5Application.single_mesh_xdmf_output_process import Factory as XdmfOutputFactory
from KratosMultiphysics.HDF5Application.single_mesh_temporal_input_process import Factory as TimeInputFactory
from KratosMultiphysics.HDF5Application.initialization_from_hdf5_process import Factory as InitializationFactory


def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestHDF5Processes(KratosUnittest.TestCase):

    def setUp(self):
        self.model = Model()
        self.buffer_size = 2

    def tearDown(self):
        self._remove_h5_files("WriteModelPart")

    def test_SingleMeshTemporalOutputInput(self):
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

    def test_SingleMeshXdmfOutputProcess(self):
        """
        Output ModelPart using SingleMeshXdmfOutputProcess.
        """
        import warnings
        try:
            with warnings.catch_warnings():
                # suppressing an import-related warning from h5py
                # problem appears when using it in a test with python >=3.6
                warnings.simplefilter('ignore', category=ImportWarning)
                import h5py
        except:
            self.skipTest("h5py not available")

        write_model_part = self._CreateNewModelPart("ModelPartXDMFOutput")
        self._InitializeModelPart(write_model_part)

        # simulate the existance of "old" files, i.e. from previous simulations
        kratos_utilities.DeleteDirectoryIfExisting("ModelPartXDMFOutput__h5_files")
        os.makedirs("ModelPartXDMFOutput__h5_files")
        with open(os.path.join("ModelPartXDMFOutput__h5_files", "ModelPartXDMFOutput.h5"), "w") as h5_file:
            pass
        with open("ModelPartXDMFOutput.xdmf", "w") as xdmf_file:
            pass

        output_settings = Parameters(r"""{
            "Parameters" : {
                "model_part_name" : "ModelPartXDMFOutput",
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
        output_process = XdmfOutputFactory(output_settings, self.model)

        self.assertTrue(os.path.isdir("ModelPartXDMFOutput__h5_files"))

        output_process.ExecuteBeforeSolutionLoop()

        for step in range(1,3):
            # Run step on write model part and dump results
            self._SimulateTimeStep(write_model_part, 0.1*step)
            output_process.ExecuteFinalizeSolutionStep()

            # the data in the h5-files is being checked in other tests
            # here only testing the xdmf-file
            ref_file_name = "ModelPartXDMFOutput_" + str(step) + ".xdmf"
            Check("ModelPartXDMFOutput.xdmf", os.path.join("xdmf_reference_files", ref_file_name))

        kratos_utilities.DeleteFileIfExisting("ModelPartXDMFOutput.xdmf")
        kratos_utilities.DeleteDirectoryIfExisting("ModelPartXDMFOutput__h5_files")

    def test_SingleMeshTemporalOutputInitialization(self):
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

def Check(output_file, reference_file):
    import KratosMultiphysics.compare_two_files_check_process as compare_process

    ## Settings string in json format
    params = KratosMultiphysics.Parameters("""{
        "reference_file_name" : \"""" + GetFilePath(reference_file) + """\",
        "output_file_name"    : \"""" + output_file + """\"
    }""")

    cmp_process = compare_process.CompareTwoFilesCheckProcess(params)

    cmp_process.ExecuteInitialize()
    cmp_process.ExecuteBeforeSolutionLoop()
    cmp_process.ExecuteInitializeSolutionStep()
    cmp_process.ExecuteFinalizeSolutionStep()
    cmp_process.ExecuteBeforeOutputStep()
    cmp_process.ExecuteAfterOutputStep()
    cmp_process.ExecuteFinalize()

if __name__ == "__main__":
    KratosUnittest.main()