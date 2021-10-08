import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp
import KratosMultiphysics.MeshingApplication as MeshingApplication
import math
from KratosMultiphysics.gid_output_process import GiDOutputProcess
from KratosMultiphysics.vtk_output_process import VtkOutputProcess
import time as time

def DotProduct(A,B):
    result = 0
    for i,j in zip(A,B):
        result += i*j
    return result

sections=[-1]
def GetSectionName(section):
    if section ==-1:
        return 'Wake3D_Wake_Auto1'
    if section == 0:
        return 'Middle_Airfoil'
    else:
        return "Section_"+str(section)

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
            "kutta_condition_rotation_angle" : 0.0,
            "wake_process_cpp_parameters":    {
                "tolerance"                     : 1e-9,
                "wake_normal"                   : [0.0,0.0,1.0],
                "wake_direction"                : [1.0,0.0,0.0],
                "is_sharp_trailing_edge"        : true,
                "switch_wake_normal"            : false,
                "count_elements_number"         : false,
                "write_elements_ids_to_file"    : false,
                "shed_wake_from_trailing_edge"  : false,
                "shedded_wake_distance"         : 12.5,
                "shedded_wake_element_size"     : 0.2,
                "echo_level": 1
            },
            "refinement_iterations": 0,
            "target_wake_h" : 1.0
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
        print('wake_normal = ', self.wake_normal)
        print('dot = ', DotProduct(self.wake_normal,self.wake_normal))
        # if( abs(DotProduct(self.wake_normal,self.wake_normal) - 1) > self.epsilon ):
        #     raise Exception('The wake normal should be a unitary vector')

        trailing_edge_model_part_name = settings["model_part_name"].GetString()
        if trailing_edge_model_part_name == "":
            err_msg = "Empty model_part_name in DefineWakeProcess3D\n"
            err_msg += "Please specify the model part that contains the trailing edge nodes"
            raise Exception(err_msg)
        self.trailing_edge_model_part = Model[trailing_edge_model_part_name]

        self.fluid_model_part = self.trailing_edge_model_part.GetRootModelPart()
        self.fluid_model_part.ProcessInfo.SetValue(CPFApp.WAKE_NORMAL,self.wake_normal)

        self.rotation_angle = settings["kutta_condition_rotation_angle"].GetDouble()
        self.fluid_model_part.ProcessInfo.SetValue(CPFApp.ROTATION_ANGLE, self.rotation_angle)

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

        self.refinement_iterations = settings["refinement_iterations"].GetInt()
        self.target_h_wake = settings["target_wake_h"].GetDouble()

    def ExecuteInitialize(self):
        # If stl available, read wake from stl and create the wake model part
        start_time = time.time()
        self.__CreateWakeModelPart()
        exe_time = time.time() - start_time
        KratosMultiphysics.Logger.PrintInfo(
            'DefineWakeProcess3D', 'Executing __CreateWakeModelPart took ', round(exe_time, 2), ' sec')
        KratosMultiphysics.Logger.PrintInfo(
            'DefineWakeProcess3D', 'Executing __CreateWakeModelPart took ', round(exe_time/60, 2), ' min')

        # for section in sections:
        #     section_model_part = self.body_model_part.GetRootModelPart().GetSubModelPart(GetSectionName(section))
        #     for condition in section_model_part.Conditions:
        #         condition.Set(KratosMultiphysics.TO_ERASE)
        #     self.body_model_part.GetRootModelPart().RemoveSubModelPart(GetSectionName(section))
        #     self.body_model_part.GetRootModelPart().RemoveConditionsFromAllLevels(KratosMultiphysics.TO_ERASE)

        start_time = time.time()
        # self.trailing_edge_model_part = self.fluid_model_part.CreateSubModelPart("Wake3D_Wake_Auto1")
        CPFApp.Define3DWakeProcess(self.trailing_edge_model_part, self.body_model_part,
                                   self.wake_model_part, self.wake_process_cpp_parameters).ExecuteInitialize()
        exe_time = time.time() - start_time
        KratosMultiphysics.Logger.PrintInfo(
            'DefineWakeProcess3D', 'Executing Define3DWakeProcess took ', round(exe_time, 2), ' sec')
        KratosMultiphysics.Logger.PrintInfo(
            'DefineWakeProcess3D', 'Executing Define3DWakeProcess took ', round(exe_time/60, 2), ' min')

        # if self.refinement_iterations > 0:
        #     for section in sections:
        #         section_model_part = self.body_model_part.GetRootModelPart().GetSubModelPart(GetSectionName(section))
        #         for condition in section_model_part.Conditions:
        #             condition.Set(KratosMultiphysics.TO_ERASE)
        #         self.body_model_part.GetRootModelPart().RemoveSubModelPart(GetSectionName(section))
        #     self.body_model_part.GetRootModelPart().RemoveConditionsFromAllLevels(KratosMultiphysics.TO_ERASE)

        #self._BlockDomain()
        self.number_of_sweeps = 1
        self.remove_modelparts = True

        # Refinement iteration loop
        for _ in range(self.refinement_iterations):
            print('self.target_h_wake = ', self.target_h_wake)
            # Option 1 - Fix domain, remesh elements intersected by wake only wake
            start_time = time.time()
            self._BlockDomain()
            exe_time = time.time() - start_time
            KratosMultiphysics.Logger.PrintInfo(
                'DefineWakeProcess3D', 'Executing _BlockDomain took ', round(exe_time, 2), ' sec')

            # Option 2 - Fix wake and wing, remesh the rest of the elements (typically to reduce number of nodes)
            # self._BlockWake()

            # Option 3 - Remesh everything, setting a metric for the wake and domain
            # self._CalculateMetricWake()

            # Option 4 - Remesh everything, setting a metric for the wake and domain
            # start_time = time.time()
            # self._BlockBodyAndRefineWake()
            # exe_time = time.time() - start_time
            # KratosMultiphysics.Logger.PrintInfo(
            #     'DefineWakeProcess3D', 'Executing _BlockBodyAndRefineWake took ', round(exe_time, 2), ' sec')
            # KratosMultiphysics.Logger.PrintInfo(
            # 'DefineWakeProcess3D', 'Executing _BlockBodyAndRefineWake took ', round(exe_time/60, 2), ' min')

            # Block body nodes
            start_time = time.time()
            CPFApp.PotentialFlowUtilities.BlockBodyNodes(self.body_model_part)
            exe_time = time.time() - start_time
            KratosMultiphysics.Logger.PrintInfo(
                'DefineWakeProcess3D', 'Executing BlockBodyNodes took ', round(exe_time, 2), ' sec')

            # Refine
            self._CallMMG()
            # self.trailing_edge_model_part = self.fluid_model_part.CreateSubModelPart("Wake3D_Wake_Auto1")
            CPFApp.Define3DWakeProcess(self.trailing_edge_model_part, self.body_model_part,
                                   self.wake_model_part, self.wake_process_cpp_parameters).ExecuteInitialize()
            self.target_h_wake /= 2.0
            # Set target for next iteration
            # if self.target_h_wake < 0.3:
            #     self.number_of_sweeps = 1
            #     self.target_h_wake /= 2.0
            # else:
            #     self.target_h_wake -= 0.2
            if self.target_h_wake < 3:
                self.number_of_sweeps = 0

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

        z = 0.003#-1.44e-3#-2.87e-3#-1e-4 #-1.44e-3
        y = 0.0
        x = -0.001#-3e-4

        # Looping over stl meshes
        for stl_mesh in wake_stl_mesh:
            for vertex in stl_mesh.points:
                node1 = self.__AddNodeToWakeModelPart(float(vertex[0]) + x, float(vertex[1]) + y, float(vertex[2]) + z )
                node2 = self.__AddNodeToWakeModelPart(float(vertex[3]) + x, float(vertex[4]) + y, float(vertex[5]) + z )
                node3 = self.__AddNodeToWakeModelPart(float(vertex[6]) + x, float(vertex[7]) + y, float(vertex[8]) + z )

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
                                            "output_interval": 1.0,
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

        vtk_output = VtkOutputProcess(self.model, KratosMultiphysics.Parameters("""
                                    {
                "model_part_name"                             : "wake_model_part",
                "output_control_type"                         : "step",
                "output_interval"                             : 1,
                "file_format"                                 : "ascii",
                "output_precision"                            : 7,
                "output_sub_model_parts"                      : false,
                "output_path"                                 : "vtk_wake_output",
                "save_output_files_in_folder"                 : true,
                "nodal_solution_step_data_variables"          : [],
                "nodal_data_value_variables"                  : ["REACTION_WATER_PRESSURE"],
                "element_data_value_variables"                : [],
                "condition_data_value_variables"              : [],
                "gauss_point_variables_extrapolated_to_nodes" : []
                                    }
                                    """)
                                )

        vtk_output.ExecuteInitialize()
        vtk_output.ExecuteBeforeSolutionLoop()
        vtk_output.ExecuteInitializeSolutionStep()
        vtk_output.PrintOutput()
        vtk_output.ExecuteFinalizeSolutionStep()
        vtk_output.ExecuteFinalize()

    def ExecuteFinalizeSolutionStep(self):
        if not self.fluid_model_part.HasSubModelPart("wake_elements_model_part"):
            raise Exception("Fluid model part does not have a wake_elements_model_part")
        else:
            self.wake_sub_model_part = self.fluid_model_part.GetSubModelPart("wake_elements_model_part")

    #     CPFApp.PotentialFlowUtilities.CheckIfWakeConditionsAreFulfilled3D(self.wake_sub_model_part, 1e-1, self.echo_level-1)
        CPFApp.PotentialFlowUtilities.ComputePotentialJump3D(self.wake_sub_model_part)

    def _BlockDomain(self):

        find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(self.body_model_part.GetRootModelPart())
        find_nodal_h.Execute()

        # for elem in self.body_model_part.GetRootModelPart().Elements:
        #     elem.Set(KratosMultiphysics.BLOCKED)

        if not self.fluid_model_part.HasSubModelPart("wake_elements_model_part"):
            raise Exception("Fluid model part does not have a wake_elements_model_part")
        else:
            self.wake_sub_model_part = self.fluid_model_part.GetSubModelPart("wake_elements_model_part")

        CPFApp.PotentialFlowUtilities.SetRefinementLevel(self.wake_sub_model_part, self.target_h_wake, self.number_of_sweeps)

        '''
        for node in self.body_model_part.GetRootModelPart().Nodes:
            node.Set(KratosMultiphysics.BLOCKED)
            this_h = node.GetValue(KratosMultiphysics.NODAL_H)
            node.SetValue(MeshingApplication.METRIC_SCALAR, this_h*1e6)
            node.SetValue(CPFApp.DEACTIVATED_WAKE, 0)

        selected_node_counter = 0
        # For 0 sweeps the loop is simple
        if self.number_of_sweeps < 1:
            for node in self.wake_sub_model_part.Nodes:
                this_h = node.GetValue(KratosMultiphysics.NODAL_H)
                if this_h > self.target_h_wake:
                    selected_node_counter += 1
                    node.Set(KratosMultiphysics.BLOCKED, False)
                    node.SetValue(MeshingApplication.METRIC_SCALAR, self.target_h_wake)

            print('Number of refined nodes = ', selected_node_counter)

        else:
            node_marker = 5
            with open("nodes_to_be_refined.dat", 'w') as node_file:
                with open("elements_to_be_refined.dat", 'w') as elem_file:
                    selected_element_counter = 0
                    for elem in self.body_model_part.GetRootModelPart().Elements:
                        if elem.GetValue(CPFApp.WAKE):
                            selected_element = False
                            for node in elem.GetNodes():
                                this_h = node.GetValue(KratosMultiphysics.NODAL_H)
                                if this_h > self.target_h_wake:
                                    if not selected_element:
                                        elem_file.write('{0:15d}\n'.format(elem.Id))
                                        selected_element_counter += 1
                                        selected_element = True
                                        elem.SetValue(CPFApp.DEACTIVATED_WAKE, 10)

                                    if node.GetValue(CPFApp.DEACTIVATED_WAKE) != node_marker:
                                        node_file.write('{0:15d}\n'.format(node.Id))
                                        selected_node_counter += 1
                                        node.SetValue(CPFApp.DEACTIVATED_WAKE, node_marker)
                                    for elem_node in elem.GetNodes():
                                            elem_node.SetValue(CPFApp.DEACTIVATED_WAKE, node_marker)


                                    # if node.Is(KratosMultiphysics.BLOCKED):
                                    #     selected_node_counter += 1
                                    #     node_file.write('{0:15d}\n'.format(node.Id))
                                    node.Set(KratosMultiphysics.BLOCKED, False)
                                    node.SetValue(MeshingApplication.METRIC_SCALAR, self.target_h_wake)

                    print('Number of refined elements = ', selected_element_counter)
                    print('Number of refined nodes = ', selected_node_counter)

                    for _ in range(self.number_of_sweeps):
                        print('node_marker = ', node_marker)
                        for elem in self.body_model_part.GetRootModelPart().Elements:
                            if (elem.GetValue(CPFApp.DEACTIVATED_WAKE) != 10):
                                selected_element = False
                                for node in elem.GetNodes():
                                    if (abs(node.GetValue(CPFApp.DEACTIVATED_WAKE) - node_marker) < 1e-3):
                                        for elem_node in elem.GetNodes():
                                            this_node_h = elem_node.GetValue(KratosMultiphysics.NODAL_H)
                                            if this_node_h > self.target_h_wake:
                                                if elem_node.Is(KratosMultiphysics.BLOCKED):
                                                    node_file.write('{0:15d}\n'.format(elem_node.Id))
                                                    selected_node_counter += 1
                                                elem_node.SetValue(MeshingApplication.METRIC_SCALAR, self.target_h_wake)
                                                elem_node.Set(KratosMultiphysics.BLOCKED, False)


                                        if not selected_element:
                                            elem.Set(KratosMultiphysics.BLOCKED, False)
                                            elem_file.write('{0:15d}\n'.format(elem.Id))
                                            selected_element_counter += 1
                                            selected_element = True
                                            elem.SetValue(CPFApp.DEACTIVATED_WAKE, 10)

                        # Marking nodes for next iteration
                        for node in self.body_model_part.GetRootModelPart().Nodes:
                            if node.IsNot(KratosMultiphysics.BLOCKED):
                                node.SetValue(CPFApp.DEACTIVATED_WAKE, node_marker + 1)

                        node_marker += 1


            print('Number of refined elements = ', selected_element_counter)
            print('Number of refined nodes = ', selected_node_counter)
        '''

        if self.remove_modelparts:
            print('Removing modelparts')
            self.body_model_part.GetRootModelPart().RemoveSubModelPart("trailing_edge_elements_model_part")
            self.body_model_part.GetRootModelPart().RemoveSubModelPart("wake_elements_model_part")
            # self.body_model_part.GetRootModelPart().RemoveSubModelPart("wake_model_part")
            # self.body_model_part.GetRootModelPart().RemoveSubModelPart("Wake3D_Wake_Auto1")

    # Option 4 - Remesh everything, setting a metric for the wake and domain
    def _BlockBodyAndRefineWake(self):

        # Compute nodal H
        find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(self.body_model_part.GetRootModelPart())
        find_nodal_h.Execute()

        # Unblock all nodes and set their H
        for node in self.body_model_part.GetRootModelPart().Nodes:
            node.Set(KratosMultiphysics.BLOCKED, False)
            this_h = node.GetValue(KratosMultiphysics.NODAL_H)
            node.SetValue(MeshingApplication.METRIC_SCALAR, this_h)

        # Set a finer H at the wake
        for elem in self.body_model_part.GetRootModelPart().Elements:
            if elem.GetValue(CPFApp.WAKE):
                for node in elem.GetNodes():
                    this_h = node.GetValue(KratosMultiphysics.NODAL_H)
                    if this_h > self.target_h_wake:
                        node.SetValue(MeshingApplication.METRIC_SCALAR, self.target_h_wake)

        # Block all nodes from the body
        for node in self.body_model_part.Nodes:
            node.Set(KratosMultiphysics.BLOCKED)
            this_h = node.GetValue(KratosMultiphysics.NODAL_H)
            node.SetValue(MeshingApplication.METRIC_SCALAR, this_h*1e6)

        if self.remove_modelparts:
            print('Removing modelparts')
            self.body_model_part.GetRootModelPart().RemoveSubModelPart("trailing_edge_elements_model_part")
            self.body_model_part.GetRootModelPart().RemoveSubModelPart("wake_elements_model_part")
            # self.body_model_part.GetRootModelPart().RemoveSubModelPart("wake_model_part")
            # self.body_model_part.GetRootModelPart().RemoveSubModelPart("Wake3D_Wake_Auto1")
            self.main_model_part.RemoveSubModelPart('wake_sub_model_part')
            self.main_model_part.RemoveSubModelPart('trailing_edge_sub_model_part')
            self.main_model_part.RemoveSubModelPart('fluid_computational_model_part')

    def _CallMMG(self):

        ini_time=time.time()

        mmg_parameters = KratosMultiphysics.Parameters("""
        {
            "discretization_type"              : "STANDARD",
            "save_external_files"              : false,
            "initialize_entities"              : false,
            "preserve_flags"                   : false,
            "save_mdpa_file"                       : false,
            "echo_level"                       : 3,
            "interpolate_nodal_values"             : false
        }
        """)

        mmg_process =MeshingApplication.MmgProcess3D(self.body_model_part.GetRootModelPart(), mmg_parameters)
        mmg_process.Execute()

        KratosMultiphysics.Logger.PrintInfo('DefineWakeProcess','Remesh time: ',time.time()-ini_time)