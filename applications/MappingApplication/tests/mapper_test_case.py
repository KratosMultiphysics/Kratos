import KratosMultiphysics as KM
from KratosMultiphysics import KratosUnittest
data_comm = KM.Testing.GetDefaultDataCommunicator()
import os
from KratosMultiphysics import from_json_check_result_process
from KratosMultiphysics import json_output_process
from KratosMultiphysics import vtk_output_process

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), "mdpa_files", fileName)

class MapperTestCase(KratosUnittest.TestCase):
    """This test case is designed for performing multiple test with the same modelparts
    this way the partitioning has to be done only once
    The values in the ModelParts are re-initialized after every
    """
    @classmethod
    def setUpModelParts(cls, mdpa_file_name_origin, mdpa_file_name_destination):
        cls.input_file_origin      = GetFilePath(mdpa_file_name_origin)
        cls.input_file_destination = GetFilePath(mdpa_file_name_destination)

        cls.current_model = KM.Model()
        cls.model_part_origin = cls.current_model.CreateModelPart("origin")
        cls.model_part_destination = cls.current_model.CreateModelPart("destination")

        # list of variables involved in the Mapper-Tests
        cls.model_part_origin.AddNodalSolutionStepVariable(KM.PRESSURE)
        cls.model_part_origin.AddNodalSolutionStepVariable(KM.FORCE)
        cls.model_part_origin.AddNodalSolutionStepVariable(KM.DISPLACEMENT)

        cls.model_part_destination.AddNodalSolutionStepVariable(KM.TEMPERATURE)
        cls.model_part_destination.AddNodalSolutionStepVariable(KM.VELOCITY)
        cls.model_part_destination.AddNodalSolutionStepVariable(KM.REACTION)
        cls.model_part_destination.AddNodalSolutionStepVariable(KM.MESH_DISPLACEMENT)

        cls.model_part_origin.ProcessInfo[KM.DOMAIN_SIZE] = 3 # needed for the partitioner!
        cls.model_part_destination.ProcessInfo[KM.DOMAIN_SIZE] = 3 # needed for the partitioner!
        cls.model_part_origin.ProcessInfo[KM.TIME] = 0.0 # needed for the check-processes
        cls.model_part_destination.ProcessInfo[KM.TIME] = 0.0 # needed for the check-processes
        cls.model_part_origin.ProcessInfo[KM.DELTA_TIME] = 1.0 # needed for the check-processes
        cls.model_part_destination.ProcessInfo[KM.DELTA_TIME] = 1.0 # needed for the check-processes

        if data_comm.IsDistributed():
            ReadDistributedModelPart(cls.model_part_origin, cls.input_file_origin)
            ReadDistributedModelPart(cls.model_part_destination, cls.input_file_destination)
        else:
            ReadModelPart(cls.model_part_origin, cls.input_file_origin)
            ReadModelPart(cls.model_part_destination, cls.input_file_destination)

    def setUp(self):
        # reset the ModelPart
        # initialize it with random values, such that I am not accidentially
        # checking against 0.0
        default_scalar = -12345.6789
        default_vector = KM.Vector([111.222, -222.999, 333.444])

        for node in self.model_part_origin.Nodes:
            node.SetSolutionStepValue(KM.PRESSURE, default_scalar)
            node.SetSolutionStepValue(KM.FORCE, default_vector)
            node.SetValue(KM.PRESSURE, default_scalar)
            node.SetValue(KM.FORCE, default_vector)

        for element in self.model_part_origin.Elements:
            element.SetValue(KM.PRESSURE, default_scalar)
            element.SetValue(KM.FORCE, default_vector)

        for condition in self.model_part_origin.Conditions:
            condition.SetValue(KM.PRESSURE, default_scalar)
            condition.SetValue(KM.FORCE, default_vector)


        for node in self.model_part_destination.Nodes:
            node.SetSolutionStepValue(KM.TEMPERATURE, default_scalar)
            node.SetSolutionStepValue(KM.VELOCITY, default_vector)
            node.SetValue(KM.TEMPERATURE, default_scalar)
            node.SetValue(KM.VELOCITY, default_vector)

        for element in self.model_part_destination.Elements:
            element.SetValue(KM.TEMPERATURE, default_scalar)
            element.SetValue(KM.VELOCITY, default_vector)

        for condition in self.model_part_destination.Conditions:
            condition.SetValue(KM.TEMPERATURE, default_scalar)
            condition.SetValue(KM.VELOCITY, default_vector)

def ReadModelPart(model_part, mdpa_file_name):
    import_flags = KM.ModelPartIO.READ | KM.ModelPartIO.SKIP_TIMER

    KM.ModelPartIO(mdpa_file_name, import_flags).ReadModelPart(model_part)

def ReadDistributedModelPart(model_part, mdpa_file_name):
    from KratosMultiphysics.mpi import distributed_import_model_part_utility
    model_part.AddNodalSolutionStepVariable(KM.PARTITION_INDEX)

    importer_settings = KM.Parameters("""{
        "model_import_settings": {
            "input_type": "mdpa",
            "input_filename": \"""" + mdpa_file_name + """\",
            "partition_in_memory" : true
        },
        "echo_level" : 0
    }""")

    model_part_import_util = distributed_import_model_part_utility.DistributedImportModelPartUtility(model_part, importer_settings)
    model_part_import_util.ImportModelPart()
    model_part_import_util.CreateCommunicators()

def GetFullModelPartName(model_part):
    full_name = model_part.Name
    if model_part.IsSubModelPart():
        full_name = GetFullModelPartName(model_part.GetParentModelPart()) + "." + full_name

    return full_name

def OutputReferenceSolution(model_part, variable, file_name):
    full_model_part_name = GetFullModelPartName(model_part)

    if data_comm.IsDistributed():
        raise Exception("Writing of reference results in not possible in MPI!")
    KM.Logger.PrintWarning('BladeMappingTests', 'Writing reference solution for ModelPart "{}"; Variable "{}"; FileName "{}"'.format(full_model_part_name, variable.Name(), file_name))

    output_parameters = KM.Parameters("""{
        "output_variables"     : [],
        "output_file_name"     : "",
        "model_part_name"      : "",
        "time_frequency"       : 0.00
    }""")

    output_parameters["output_variables"].Append(variable.Name())
    output_parameters["output_file_name"].SetString(file_name + ".json")
    output_parameters["model_part_name"].SetString(full_model_part_name)

    output_proc = json_output_process.JsonOutputProcess(model_part.GetModel(), output_parameters)
    output_proc.ExecuteInitialize()
    output_proc.ExecuteBeforeSolutionLoop()
    output_proc.ExecuteFinalizeSolutionStep()

def CheckHistoricalNonUniformValues(model_part, variable, file_name, output_reference_solution=False):
    if output_reference_solution:
        OutputReferenceSolution(model_part, variable, file_name)
    else:
        full_model_part_name = GetFullModelPartName(model_part)

        check_parameters = KM.Parameters("""{
            "check_variables"           : [],
            "input_file_name"           : "",
            "model_part_name"           : "",
            "tolerance"                 : 1e-6,
            "relative_tolerance"        : 1e-9,
            "time_frequency"            : 0.00,
            "check_only_local_entities" : false
        }""")

        check_parameters["check_variables"].Append(variable.Name())
        check_parameters["input_file_name"].SetString(file_name + ".json")
        check_parameters["model_part_name"].SetString(full_model_part_name)

        check_proc = from_json_check_result_process.FromJsonCheckResultProcess(model_part.GetModel(), check_parameters)
        check_proc.ExecuteInitialize()
        check_proc.ExecuteFinalizeSolutionStep()

def VtkOutputNodesHistorical(model_part, variable, prefix=""):
    vtk_parameters = KM.Parameters("""{
        "model_part_name"                    : \"""" + model_part.Name + """\",
        "file_format"                        : "binary",
        "output_control_type"                : "step",
        "output_sub_model_parts"             : false,
        "folder_name"                        : \"""" + "VTK_Output_" + prefix + """\",
        "save_output_files_in_folder"        : true,
            "custom_name_prefix"             : \"""" + prefix + "_" + """\",
        "nodal_solution_step_data_variables" : [\"""" + variable.Name() + """\"]
    }""")

    vtk_output_process.VtkOutputProcess(model_part.GetModel(),vtk_parameters).PrintOutput()
