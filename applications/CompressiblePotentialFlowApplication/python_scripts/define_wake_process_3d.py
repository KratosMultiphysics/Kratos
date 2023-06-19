import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp
import math
from KratosMultiphysics.gid_output_process import GiDOutputProcess
import time as time

def DotProduct(A,B):
    result = 0
    for i,j in zip(A,B):
        result += i*j
    return result

def CrossProduct(A, B):
    C = KratosMultiphysics.Vector(3)
    C[0] = A[1]*B[2]-A[2]*B[1]
    C[1] = A[2]*B[0]-A[0]*B[2]
    C[2] = A[0]*B[1]-A[1]*B[0]
    return C

def Factory(settings, Model):
    if(not isinstance(settings, KratosMultiphysics.Parameters)):
        raise Exception(
            "expected input shall be a Parameters object, encapsulating a json string")

    return DefineWakeProcess3D(Model, settings["Parameters"])

class DefineWakeProcess3D(KratosMultiphysics.Process):
    def __init__(self, Model, settings):
        KratosMultiphysics.Process.__init__(self)

        # Check default settings
        default_settings = KratosMultiphysics.Parameters(r'''{
            "model_part_name": "",
            "body_model_part_name": "",
            "wake_stl_file_name" : "",
            "output_wake": false,
            "wake_process_cpp_parameters":    {
                "tolerance"                     : 1e-9,
                "wake_normal"                   : [0.0,0.0,1.0],
                "wake_direction"                : [1.0,0.0,0.0],
                "switch_wake_normal"            : false,
                "count_elements_number"         : false,
                "write_elements_ids_to_file"    : false,
                "shed_wake_from_trailing_edge"  : false,
                "shedded_wake_distance"         : 12.5,
                "shedded_wake_element_size"     : 0.2,
                "echo_level": 1
            }
        }''')
        settings.ValidateAndAssignDefaults(default_settings)

        self.model = Model

        self.output_wake = settings["output_wake"].GetBool()
        self.wake_process_cpp_parameters = settings["wake_process_cpp_parameters"]

        self.epsilon = self.wake_process_cpp_parameters["tolerance"].GetDouble()
        self.shedded_wake_distance = self.wake_process_cpp_parameters["shedded_wake_distance"].GetDouble()
        self.shedded_wake_element_size = self.wake_process_cpp_parameters["shedded_wake_element_size"].GetDouble()
        self.echo_level = self.wake_process_cpp_parameters["echo_level"].GetInt()

        # This is a reference value for the wake normal
        self.wake_normal = self.wake_process_cpp_parameters["wake_normal"].GetVector()
        if( abs(DotProduct(self.wake_normal,self.wake_normal) - 1) > self.epsilon ):
            raise Exception('The wake normal should be a unitary vector')

        trailing_edge_model_part_name = settings["model_part_name"].GetString()
        if trailing_edge_model_part_name == "":
            err_msg = "Empty model_part_name in DefineWakeProcess3D\n"
            err_msg += "Please specify the model part that contains the trailing edge nodes"
            raise Exception(err_msg)
        self.trailing_edge_model_part = Model[trailing_edge_model_part_name]

        self.fluid_model_part = self.trailing_edge_model_part.GetRootModelPart()
        self.fluid_model_part.ProcessInfo.SetValue(CPFApp.WAKE_NORMAL,self.wake_normal)

        if self.wake_process_cpp_parameters.Has("wake_direction"):
            self.wake_direction = self.wake_process_cpp_parameters["wake_direction"].GetVector()
        else:
            self.wake_direction = self.fluid_model_part.ProcessInfo.GetValue(CPFApp.FREE_STREAM_VELOCITY_DIRECTION)
            self.wake_process_cpp_parameters.AddEmptyValue("wake_direction")
            self.wake_process_cpp_parameters["wake_direction"].SetVector(self.wake_direction)

        body_model_part_name = settings["body_model_part_name"].GetString()
        if body_model_part_name == "":
            err_msg = "Empty model_part_name in DefineWakeProcess3D\n"
            err_msg += "Please specify the model part that contains the body nodes"
            raise Exception(err_msg)
        self.body_model_part = Model[body_model_part_name]

        self.shed_wake_from_trailing_edge = self.wake_process_cpp_parameters["shed_wake_from_trailing_edge"].GetBool()
        if self.shed_wake_from_trailing_edge:
            warn_msg = 'Generating the wake automatically from the trailing edge.'
            KratosMultiphysics.Logger.PrintWarning('::[DefineWakeProcess3D]::', warn_msg)

        self.wake_stl_file_name = settings["wake_stl_file_name"].GetString()
        if self.wake_stl_file_name == "":
            self.shed_wake_from_trailing_edge = True
            warn_msg = 'Empty wake_stl_file_name in DefineWakeProcess3D,'
            warn_msg += ' generating the wake automatically from the trailing edge.'
            KratosMultiphysics.Logger.PrintWarning('::[DefineWakeProcess3D]::', warn_msg)

    def ExecuteInitialize(self):
        # If stl available, read wake from stl and create the wake model part
        start_time = time.time()
        self.__CreateWakeModelPart()
        exe_time = time.time() - start_time
        KratosMultiphysics.Logger.PrintInfo(
            'DefineWakeProcess3D', 'Executing __CreateWakeModelPart took ', round(exe_time, 2), ' sec')
        KratosMultiphysics.Logger.PrintInfo(
            'DefineWakeProcess3D', 'Executing __CreateWakeModelPart took ', round(exe_time/60, 2), ' min')

        start_time = time.time()
        CPFApp.Define3DWakeProcess(self.trailing_edge_model_part, self.body_model_part,
                                   self.wake_model_part, self.wake_process_cpp_parameters).ExecuteInitialize()
        exe_time = time.time() - start_time
        KratosMultiphysics.Logger.PrintInfo(
            'DefineWakeProcess3D', 'Executing Define3DWakeProcess took ', round(exe_time, 2), ' sec')
        KratosMultiphysics.Logger.PrintInfo(
            'DefineWakeProcess3D', 'Executing Define3DWakeProcess took ', round(exe_time/60, 2), ' min')

        # # Output the wake in GiD for visualization
        if(self.output_wake):
            self.__VisualizeWake()

    # This function imports the stl file containing the wake and creates the wake model part out of it.
    # TODO: implement an automatic generation of the wake
    def __CreateWakeModelPart(self):
        self.wake_model_part = self.model.CreateModelPart("wake_model_part")
        self.dummy_property = self.wake_model_part.Properties[0]
        self.node_id = 1
        self.elem_id = 1
        if not self.shed_wake_from_trailing_edge:
            self.__ReadWakeStlModelFromFile()

    def __ReadWakeStlModelFromFile(self):
        KratosMultiphysics.Logger.PrintInfo('DefineWakeProcess3D', 'Reading wake from stl file')
        from stl import mesh #this requires numpy-stl
        wake_stl_mesh = mesh.Mesh.from_multi_file(self.wake_stl_file_name)

        z = 0.0#-1e-4

        # Looping over stl meshes
        for stl_mesh in wake_stl_mesh:
            for vertex in stl_mesh.points:
                node1 = self.__AddNodeToWakeModelPart(float(vertex[0]), float(vertex[1]), float(vertex[2]) + z )
                node2 = self.__AddNodeToWakeModelPart(float(vertex[3]), float(vertex[4]), float(vertex[5]) + z )
                node3 = self.__AddNodeToWakeModelPart(float(vertex[6]), float(vertex[7]), float(vertex[8]) + z )

                side1 = node2 - node1
                side2 = node3 - node1
                face_normal = CrossProduct(side1,side2)

                normal_projection = DotProduct(face_normal, self.wake_normal)

                if normal_projection >  0.0:
                    self.__AddElementToWakeModelPart(node1.Id, node2.Id, node3.Id)
                else:
                    self.__AddElementToWakeModelPart(node1.Id, node3.Id, node2.Id)

    def __AddNodeToWakeModelPart(self, x, y, z):
        node = self.wake_model_part.CreateNewNode(self.node_id, x, y, z)
        self.node_id +=1
        return node

    def __AddElementToWakeModelPart(self, id1, id2, id3):
        self.wake_model_part.CreateNewElement("Element3D3N", self.elem_id,  [id1, id2, id3], self.dummy_property)
        self.elem_id += 1

    def __VisualizeWake(self):
        # To visualize the wake
        number_of_nodes = self.fluid_model_part.NumberOfNodes()
        number_of_elements = self.fluid_model_part.NumberOfElements()

        node_id = number_of_nodes + 1
        for node in self.wake_model_part.Nodes:
            node.Id = node_id
            node.SetValue(KratosMultiphysics.REACTION_WATER_PRESSURE, 1.0)
            node_id += 1

        counter = number_of_elements + 1
        for elem in self.wake_model_part.Elements:
            elem.Id = counter
            counter +=1

        output_file = "representation_of_wake"
        gid_output =  GiDOutputProcess(self.wake_model_part,
                                output_file,
                                KratosMultiphysics.Parameters("""
                                    {
                                        "result_file_configuration": {
                                            "gidpost_flags": {
                                                "GiDPostMode": "GiD_PostAscii",
                                                "WriteDeformedMeshFlag": "WriteUndeformed",
                                                "WriteConditionsFlag": "WriteConditions",
                                                "MultiFileFlag": "SingleFile"
                                            },
                                            "file_label": "time",
                                            "output_control_type": "step",
                                            "output_frequency": 1.0,
                                            "body_output": true,
                                            "node_output": false,
                                            "skin_output": false,
                                            "plane_output": [],
                                            "nodal_results": [],
                                            "nodal_nonhistorical_results": ["REACTION_WATER_PRESSURE"],
                                            "nodal_flags_results": [],
                                            "gauss_point_results": [],
                                            "additional_list_files": []
                                        }
                                    }
                                    """)
                                )

        gid_output.ExecuteInitialize()
        gid_output.ExecuteBeforeSolutionLoop()
        gid_output.ExecuteInitializeSolutionStep()
        gid_output.PrintOutput()
        gid_output.ExecuteFinalizeSolutionStep()
        gid_output.ExecuteFinalize()

    def ExecuteFinalizeSolutionStep(self):
        if not self.fluid_model_part.HasSubModelPart("wake_elements_model_part"):
            raise Exception("Fluid model part does not have a wake_elements_model_part")
        else:
            self.wake_sub_model_part = self.fluid_model_part.GetSubModelPart("wake_elements_model_part")

        CPFApp.PotentialFlowUtilities.CheckIfWakeConditionsAreFulfilled3D(self.wake_sub_model_part, 1e-1, self.echo_level-1)
        CPFApp.PotentialFlowUtilities.ComputePotentialJump3D(self.wake_sub_model_part)