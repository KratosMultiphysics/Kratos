import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp
import KratosMultiphysics.MeshingApplication as MeshingApplication
import math
from KratosMultiphysics.gid_output_process import GiDOutputProcess
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
            "wake_normal": [0.0,0.0,1.0],
            "refinement_iterations": 0,
            "target_wake_h" : 1.0,
            "output_wake": false,
            "epsilon": 1e-9
        }''')
        settings.ValidateAndAssignDefaults(default_settings)

        self.model = Model

        trailing_edge_model_part_name = settings["model_part_name"].GetString()
        if trailing_edge_model_part_name == "":
            err_msg = "Empty model_part_name in DefineWakeProcess3D\n"
            err_msg += "Please specify the model part that contains the trailing edge nodes"
            raise Exception(err_msg)
        self.trailing_edge_model_part = Model[trailing_edge_model_part_name]

        body_model_part_name = settings["body_model_part_name"].GetString()
        if body_model_part_name == "":
            err_msg = "Empty model_part_name in DefineWakeProcess3D\n"
            err_msg += "Please specify the model part that contains the body nodes"
            raise Exception(err_msg)
        self.body_model_part = Model[body_model_part_name]

        self.wake_stl_file_name = settings["wake_stl_file_name"].GetString()
        if self.wake_stl_file_name == "":
            err_msg = "Empty wake_stl_file_name in DefineWakeProcess3D\n"
            err_msg += "Please specify the stl file name that contains the wake surface nodes"
            raise Exception(err_msg)

        self.epsilon = settings["epsilon"].GetDouble()
        self.output_wake = settings["output_wake"].GetBool()

        # For now plane wake surfaces are considered.
        # TODO: Generalize this to curved wake surfaces
        self.wake_normal = settings["wake_normal"].GetVector()
        if( abs(DotProduct(self.wake_normal,self.wake_normal) - 1) > self.epsilon ):
            raise Exception('The wake normal should be a unitary vector')

        self.fluid_model_part = self.trailing_edge_model_part.GetRootModelPart()
        self.fluid_model_part.ProcessInfo.SetValue(CPFApp.WAKE_NORMAL,self.wake_normal)

        self.refinement_iterations = settings["refinement_iterations"].GetInt()
        self.target_h_wake = settings["target_wake_h"].GetDouble()

    def ExecuteInitialize(self):

        # Read wake from stl and create the wake model part
        self.__CreateWakeModelPart()

        # for section in sections:
        #     section_model_part = self.body_model_part.GetRootModelPart().GetSubModelPart(GetSectionName(section))
        #     for condition in section_model_part.Conditions:
        #         condition.Set(KratosMultiphysics.TO_ERASE)
        #     self.body_model_part.GetRootModelPart().RemoveSubModelPart(GetSectionName(section))
        #     self.body_model_part.GetRootModelPart().RemoveConditionsFromAllLevels(KratosMultiphysics.TO_ERASE)

        start_time = time.time()
        # self.trailing_edge_model_part = self.fluid_model_part.CreateSubModelPart("Wake3D_Wake_Auto1")
        number_of_nodes = self.trailing_edge_model_part.NumberOfNodes()
        print('number_of_nodes = ', number_of_nodes)
        CPFApp.Define3DWakeProcess(self.trailing_edge_model_part, self.body_model_part, self.wake_model_part, self.epsilon, self.wake_normal).ExecuteInitialize()
        exe_time = time.time() - start_time
        print('Executing Define3DWakeProcess took ' + str(round(exe_time, 2)) + ' sec')
        print('Executing Define3DWakeProcess took ' + str(round(exe_time/60, 2)) + ' min')

        if self.refinement_iterations > 0:
            for section in sections:
                section_model_part = self.body_model_part.GetRootModelPart().GetSubModelPart(GetSectionName(section))
                for condition in section_model_part.Conditions:
                    condition.Set(KratosMultiphysics.TO_ERASE)
                self.body_model_part.GetRootModelPart().RemoveSubModelPart(GetSectionName(section))
            self.body_model_part.GetRootModelPart().RemoveConditionsFromAllLevels(KratosMultiphysics.TO_ERASE)

        #self._BlockDomain()
        self.number_of_sweeps = 2
        self.remove_modelparts = True

        for _ in range(self.refinement_iterations):
            print('self.target_h_wake = ', self.target_h_wake)
            # Option 1 - Fix domain, remesh elements intersected by wake only wake
            self._BlockDomain()
            # Option 2 - Fix wake and wing, remesh the rest of the elements (typically to reduce number of nodes)
            # self._BlockWake()
            # Option 3 - Remesh everything, setting a metric for the wake and domain
            # self._CalculateMetricWake()
            self._CallMMG()
            self.trailing_edge_model_part = self.fluid_model_part.CreateSubModelPart("Wake3D_Wake_Auto1")
            CPFApp.Define3DWakeProcess(self.trailing_edge_model_part, self.body_model_part, self.wake_model_part, self.epsilon, self.wake_normal).ExecuteInitialize()
            #self.target_h_wake /= 2.0
            if self.target_h_wake < 0.3:
                self.number_of_sweeps = 1
                self.target_h_wake /= 2.0
            else:
                self.target_h_wake -= 0.2

        # self.remove_modelparts = False
        # self._BlockDomain()

        # self.__SetWakeAndSpanDirections()
        # # Save the trailing edge and wing tip nodes for further computations
        # self.__MarkTrailingEdgeNodes()
        # # Save the lower surface normals to help mark kutta elements later on
        # self.__ComputeLowerSurfaceNormals()
        # # Read wake from stl and create the wake model part
        # self.__CreateWakeModelPart()
        # # Check which elements are cut and mark them as wake
        # self.__MarkWakeElements()
        # # Mark the elements touching the trailing edge from below as kutta
        # self.__MarkKuttaElements()
        # # Output the wake in GiD for visualization
        # if(self.output_wake):
        #     self.__VisualizeWake()

    def __SetWakeAndSpanDirections(self):
        free_stream_velocity = self.fluid_model_part.ProcessInfo.GetValue(CPFApp.FREE_STREAM_VELOCITY)
        if(free_stream_velocity.Size() != 3):
            raise Exception('The free stream velocity should be a vector with 3 components!')
        self.wake_direction = KratosMultiphysics.Vector(3)
        vnorm = math.sqrt(
            free_stream_velocity[0]**2 + free_stream_velocity[1]**2 + free_stream_velocity[2]**2)
        self.wake_direction[0] = free_stream_velocity[0]/vnorm
        self.wake_direction[1] = free_stream_velocity[1]/vnorm
        self.wake_direction[2] = free_stream_velocity[2]/vnorm

        # span_direction = wake_normal x wake_direction
        self.span_direction = KratosMultiphysics.Vector(3)
        self.span_direction[0] = self.wake_normal[1] * self.wake_direction[2] - self.wake_normal[2] * self.wake_direction[1]
        self.span_direction[1] = self.wake_normal[2] * self.wake_direction[0] - self.wake_normal[0] * self.wake_direction[2]
        self.span_direction[0] = self.wake_normal[0] * self.wake_direction[1] - self.wake_normal[1] * self.wake_direction[0]

    def __MarkTrailingEdgeNodes(self):
        # Mark trailing edge nodes
        for node in self.trailing_edge_model_part.Nodes:
            node.SetValue(CPFApp.TRAILING_EDGE, True)

    def __ComputeLowerSurfaceNormals(self):
        for cond in self.body_model_part.Conditions:
            # The surface normal points outisde the domain
            surface_normal = cond.GetGeometry().Normal()
            norm = math.sqrt(surface_normal[0]**2 + surface_normal[1]**2 + surface_normal[2]**2)
            if abs(norm) < self.epsilon:
                raise Exception('The norm of the condition ', cond.Id , ' should be larger than 0.')
            # Normalizing normal vector
            surface_normal /= norm
            projection = DotProduct(surface_normal, self.wake_normal)
            # The surface normal in the same direction as the wake normal belongs to the lower_surface
            if(projection > 0.0):
                for node in cond.GetNodes():
                    node.SetValue(KratosMultiphysics.NORMAL,surface_normal)
                    node.SetValue(CPFApp.LOWER_SURFACE, True)
            else:
                for node in cond.GetNodes():
                    node.SetValue(CPFApp.UPPER_SURFACE, True)

    # This function imports the stl file containing the wake and creates the wake model part out of it.
    # TODO: implement an automatic generation of the wake
    def __CreateWakeModelPart(self):
        from stl import mesh #this requires numpy-stl
        wake_stl_mesh = mesh.Mesh.from_multi_file(self.wake_stl_file_name)
        self.wake_model_part = self.model.CreateModelPart("wake_model_part")

        dummy_property = self.wake_model_part.Properties[0]
        node_id = 1
        elem_id = 1
        #z= -1e-5
        z= 0.0#-1e-4

        # Looping over stl meshes
        for stl_mesh in wake_stl_mesh:
            for vertex in stl_mesh.points:
                node1 = self.wake_model_part.CreateNewNode(node_id, float(vertex[0]), float(vertex[1]), float(vertex[2]) + z )
                node_id+=1
                node2 = self.wake_model_part.CreateNewNode(node_id, float(vertex[3]), float(vertex[4]), float(vertex[5]) + z )
                node_id+=1
                node3 = self.wake_model_part.CreateNewNode(node_id, float(vertex[6]), float(vertex[7]), float(vertex[8]) + z )
                node_id+=1

                self.wake_model_part.CreateNewElement("Element3D3N", elem_id,  [
                                              node1.Id, node2.Id, node3.Id], dummy_property)
                elem_id += 1

        # compute_wake_normal = False
        # for elem in self.wake_model_part.Elements:
        #     if not compute_wake_normal:
        #         #print(elem.GetGeometry().Normal())
        #         wake_normal = elem.GetGeometry().Normal()
        #         magnitude = -math.sqrt(DotProduct(wake_normal,wake_normal))
        #         compute_wake_normal = True
        #         for i in range(len(self.wake_normal)):
        #             self.wake_normal[i] = wake_normal[i]  / magnitude
        #             print("{:.15f}".format(self.wake_normal[i]))
        #         print('self.wake_normal = ', self.wake_normal)



    # Check which elements are cut and mark them as wake
    def __MarkWakeElements(self):
        KratosMultiphysics.Logger.PrintInfo('...Selecting wake elements...')

        # Mark cut elements and compute elemental distances
        # Attention: Note that in this process a negative distance is assigned to nodes
        # laying on the wake. In 2D it is done viceversa.
        distance_calculator = KratosMultiphysics.CalculateDistanceToSkinProcess3D(
            self.fluid_model_part, self.wake_model_part)
        distance_calculator.Execute()

        # Lists to store trailing edge and wake elements id
        self.trailing_edge_element_id_list = []
        self.wake_element_id_list = []

        for elem in self.fluid_model_part.Elements:
            # Mark and save the elements touching the trailing edge
            self.__MarkTrailingEdgeElement(elem)

            # Cut elements are wake
            if(elem.Is(KratosMultiphysics.TO_SPLIT)):
                # Mark wake element
                elem.SetValue(CPFApp.WAKE, True)
                self.wake_element_id_list.append(elem.Id)
                # Save wake elemental distances
                wake_elemental_distances = elem.GetValue(KratosMultiphysics.ELEMENTAL_DISTANCES)
                # Check tolerances
                for i in range(len(wake_elemental_distances)):
                    if(abs(wake_elemental_distances[i]) < self.epsilon ):
                        if(wake_elemental_distances[i] < 0.0):
                            wake_elemental_distances[i] = -self.epsilon
                        else:
                            wake_elemental_distances[i] = self.epsilon
                # Save wake elemental distances
                elem.SetValue(CPFApp.WAKE_ELEMENTAL_DISTANCES,wake_elemental_distances)
                # Save wake nodal distances
                counter = 0
                for node in elem.GetNodes():
                    node.SetValue(CPFApp.WAKE_DISTANCE,wake_elemental_distances[counter])
                    counter += 1

        self.__SaveTrailingEdgeElements()
        self.__SaveWakeElements()
        KratosMultiphysics.Logger.PrintInfo('...Selecting wake elements finished...')

    def __MarkTrailingEdgeElement(self, elem):
        # This function marks the elements touching the trailing
        # edge and saves them in the trailing_edge_element_id_list
        # for further computations
        for elnode in elem.GetNodes():
            if(elnode.GetValue(CPFApp.TRAILING_EDGE)):
                elem.SetValue(CPFApp.TRAILING_EDGE_ELEMENT, True)
                self.trailing_edge_element_id_list.append(elem.Id)
                break

    def __SaveTrailingEdgeElements(self):
        # This function stores the trailing edge elements in its submodelpart.
        if(self.fluid_model_part.HasSubModelPart("trailing_edge_elements_model_part")):
            for elem in self.trailing_edge_elements_model_part.Elements:
                elem.Set(KratosMultiphysics.TO_ERASE)
            self.trailing_edge_elements_model_part.RemoveElements(KratosMultiphysics.TO_ERASE)
        else:
            self.trailing_edge_elements_model_part = self.fluid_model_part.CreateSubModelPart("trailing_edge_elements_model_part")
        self.trailing_edge_elements_model_part.AddElements(self.trailing_edge_element_id_list)

    def __SaveWakeElements(self):
        # This function stores the wake elements in its submodelpart.
        if(self.fluid_model_part.HasSubModelPart("wake_elements_model_part")):
            self.wake_sub_model_part = self.fluid_model_part.GetSubModelPart("wake_elements_model_part")
            for elem in self.wake_sub_model_part.Elements:
                elem.Set(KratosMultiphysics.TO_ERASE)
            self.wake_sub_model_part.RemoveElements(KratosMultiphysics.TO_ERASE)
        else:
            self.wake_sub_model_part = self.fluid_model_part.CreateSubModelPart("wake_elements_model_part")
        self.wake_sub_model_part.AddElements(self.wake_element_id_list)

    def __MarkKuttaElements(self):
        # This function selects the kutta elements. Kutta elements
        # are touching the trailing edge from below.

        # Loop over elements touching the trailing edge
        for elem in self.trailing_edge_elements_model_part.Elements:
            trailing_edge_node, number_of_non_te_nodes = self.__GetATrailingEdgeNodeAndNumberOfNonTENodes(elem)
            nodal_distances = self.__ComputeNodalDistancesToWakeAndLowerSurface(elem, trailing_edge_node, number_of_non_te_nodes)
            self.__CheckIfKuttaElement(elem, nodal_distances, number_of_non_te_nodes)

    def __GetATrailingEdgeNodeAndNumberOfNonTENodes(self,elem):
        # This function returns a trailing edge node (note that
        # an element may have more than one trailing edge node)
        # and the number of nodes that are not trailing edge.
        number_of_te_nodes = 0
        for elnode in elem.GetNodes():
            if(elnode.GetValue(CPFApp.TRAILING_EDGE)):
                trailing_edge_node = elnode
                number_of_te_nodes += 1

        number_of_non_te_nodes = 4 - number_of_te_nodes

        return trailing_edge_node, number_of_non_te_nodes

    def __ComputeNodalDistancesToWakeAndLowerSurface(self, elem, trailing_edge_node, number_of_non_te_nodes):
        # This function computes the distance of the element nodes
        # to the wake and the wing lower surface

        # Only computing the distances of the nodes that are not trailing edge
        nodal_distances_to_te = KratosMultiphysics.Vector(number_of_non_te_nodes)
        counter = 0
        for elnode in elem.GetNodes():
            # Looping only over non traling edge nodes
            if not (elnode.GetValue(CPFApp.TRAILING_EDGE)):
                # Compute the distance vector from the trailing edge to the node
                distance_vector = elnode - trailing_edge_node

                # Compute the distance in the free stream direction
                free_stream_direction_distance = DotProduct(distance_vector, self.wake_direction)

                if(elnode.GetValue(CPFApp.UPPER_SURFACE)):
                    distance = self.epsilon
                elif(elnode.GetValue(CPFApp.LOWER_SURFACE)):
                    distance = -self.epsilon
                # Node laying either above or below the lower surface
                elif(free_stream_direction_distance < 0.0):
                    # Compute the distance in the lower surface normal direction
                    distance = DotProduct(distance_vector, trailing_edge_node.GetValue(KratosMultiphysics.NORMAL))
                # Node laying either above or below the wake
                else:
                    # Compute the distance in the wake normal direction
                    distance = DotProduct(distance_vector, self.wake_normal)

                # Nodes laying on the wake or on the lower surface have a negative distance
                if(abs(distance) < self.epsilon):
                    distance = -self.epsilon

                nodal_distances_to_te[counter] = distance
                counter += 1

        return nodal_distances_to_te

    def __CheckIfKuttaElement(self, elem, nodal_distances, number_of_non_te_nodes):
        # This function checks whether the element is kutta

        # Count number of nodes above and below the wake and lower surface
        number_of_nodes_with_positive_distance, number_of_nodes_with_negative_distance = self.__CountNodes(nodal_distances)

        # Elements with all non trailing edge nodes below the wake and the lower surface are kutta
        if(number_of_nodes_with_negative_distance > number_of_non_te_nodes - 1):
            elem.SetValue(CPFApp.KUTTA, True)
            elem.SetValue(CPFApp.WAKE, False)
            self.wake_sub_model_part.RemoveElement(elem)
        # Elements with nodes above and below the wake are wake elements
        elif(number_of_nodes_with_positive_distance > 0 and number_of_nodes_with_negative_distance > 0):
            # Wake elements touching the trailing edge are marked as structure
            # TODO: change STRUCTURE to a more meaningful variable name
            elem.Set(KratosMultiphysics.STRUCTURE)
            # Updating distance values
            if(elem.GetValue(CPFApp.WAKE)):
                wake_elemental_distances = elem.GetValue(CPFApp.WAKE_ELEMENTAL_DISTANCES)
            else:
                KratosMultiphysics.Logger.PrintInfo('...Setting non cut element to structure...', elem.Id)
                elem.SetValue(CPFApp.WAKE, True)
                wake_elemental_distances = KratosMultiphysics.Vector(4)
            counter = 0
            counter2 = 0
            for elnode in elem.GetNodes():
                if elnode.GetValue(CPFApp.TRAILING_EDGE):
                    # Setting the distance at the trailing edge to positive
                    elnode.SetValue(CPFApp.WAKE_DISTANCE,self.epsilon)
                    wake_elemental_distances[counter] = self.epsilon
                else:
                    elnode.SetValue(CPFApp.WAKE_DISTANCE,nodal_distances[counter2])
                    wake_elemental_distances[counter] = nodal_distances[counter2]
                    counter2 +=1
                counter +=1
            elem.SetValue(CPFApp.WAKE_ELEMENTAL_DISTANCES,wake_elemental_distances)
            pass
        # Elements with all non trailing edge nodes above the wake and the lower surface are normal
        elif(number_of_nodes_with_positive_distance > number_of_non_te_nodes - 1):
            elem.SetValue(CPFApp.WAKE, False)
            self.wake_sub_model_part.RemoveElement(elem)

    def __CountNodes(self, distances_to_te):
        # This function counts the number of nodes that are above and below
        # the wake and the lower surface

        # Initialize counters
        number_of_nodes_with_positive_distance = 0
        number_of_nodes_with_negative_distance = 0

        # Count how many element nodes are above and below the wake
        for nodal_distance_to_wake in distances_to_te:
            if(nodal_distance_to_wake < 0.0):
                number_of_nodes_with_negative_distance += 1
            else:
                number_of_nodes_with_positive_distance += 1

        return number_of_nodes_with_positive_distance, number_of_nodes_with_negative_distance

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

    # def ExecuteFinalizeSolutionStep(self):
    #     if not self.fluid_model_part.HasSubModelPart("wake_elements_model_part"):
    #         raise Exception("Fluid model part does not have a wake_elements_model_part")
    #     else:
    #         self.wake_sub_model_part = self.fluid_model_part.GetSubModelPart("wake_elements_model_part")

    #     elements = self.wake_sub_model_part.NumberOfElements()
    #     nodes = self.wake_sub_model_part.NumberOfNodes()

    #     print('Number of wake elements = ', elements)
    #     print('Number of wake nodes = ', nodes)
    #     print('ratio = ', elements / nodes)

    #     CPFApp.PotentialFlowUtilities.CheckIfWakeConditionsAreFulfilled3D(self.wake_sub_model_part, 1e-1, 1)
    #     CPFApp.PotentialFlowUtilities.CheckIfPressureEqualityWakeConditionsAreFulfilled3D(self.wake_sub_model_part, 1e-1, 1)

    def _BlockDomain(self):

        find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(self.body_model_part.GetRootModelPart())
        find_nodal_h.Execute()

        # for elem in self.body_model_part.GetRootModelPart().Elements:
        #     elem.Set(KratosMultiphysics.BLOCKED)

        for node in self.body_model_part.GetRootModelPart().Nodes:
            node.Set(KratosMultiphysics.BLOCKED)
            this_h = node.GetValue(KratosMultiphysics.NODAL_H)
            node.SetValue(MeshingApplication.METRIC_SCALAR, this_h*1e6)

        node_marker = 5
        with open("nodes_to_be_refined.dat", 'w') as node_file:
            with open("elements_to_be_refined.dat", 'w') as elem_file:
                selected_element_counter = 0
                selected_node_counter = 0
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

                                # print(this_h)
                                #print(elem.Id)
                                elem.Set(KratosMultiphysics.BLOCKED, False)
                                node.Set(KratosMultiphysics.BLOCKED, False)
                                node.SetValue(MeshingApplication.METRIC_SCALAR, self.target_h_wake)

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

        if self.remove_modelparts:
            self.body_model_part.GetRootModelPart().RemoveSubModelPart("trailing_edge_elements_model_part")
            self.body_model_part.GetRootModelPart().RemoveSubModelPart("wake_elements_model_part")
            self.body_model_part.GetRootModelPart().RemoveSubModelPart("Wake3D_Wake_Auto1")


    def _BlockWake(self):

        find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(self.body_model_part.GetRootModelPart())
        find_nodal_h.Execute()

        for node in self.body_model_part.GetRootModelPart().Nodes:
            this_h = node.GetValue(KratosMultiphysics.NODAL_H)*5
            node.SetValue(MeshingApplication.METRIC_SCALAR, this_h)

        for node in self.body_model_part.Nodes:
            node.Set(KratosMultiphysics.BLOCKED)

        for elem in self.body_model_part.GetRootModelPart().Elements:
            if elem.GetValue(CPFApp.WAKE):
                for node in elem.GetNodes():
                    node.Set(KratosMultiphysics.BLOCKED)

        self.body_model_part.GetRootModelPart().RemoveSubModelPart("trailing_edge_elements_model_part")
        self.body_model_part.GetRootModelPart().RemoveSubModelPart("wake_elements_model_part")

    def _CalculateMetricWake(self):

        find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(self.body_model_part.GetRootModelPart())
        find_nodal_h.Execute()

        for node in self.body_model_part.GetRootModelPart().Nodes:
            this_h = max(node.GetValue(KratosMultiphysics.NODAL_H), 1.0)
            tensor = KratosMultiphysics.Vector(6, 0.0)
            for i in range(3):
                tensor[i]=1 / (this_h*this_h)
            node.SetValue(KratosMultiphysics.MeshingApplication.METRIC_TENSOR_3D, tensor)

        tensor = KratosMultiphysics.Vector(6, 0.0)
        for i in range(3):
            tensor[i]=1 / (self.target_h_wake*self.target_h_wake)
        for elem in self.body_model_part.GetRootModelPart().Elements:
            if elem.GetValue(CPFApp.WAKE):
                for node in elem.GetNodes():
                    node.SetValue(MeshingApplication.METRIC_TENSOR_3D, tensor)

        self.body_model_part.GetRootModelPart().RemoveSubModelPart("trailing_edge_elements_model_part")
        self.body_model_part.GetRootModelPart().RemoveSubModelPart("wake_elements_model_part")

    def _CallMMG(self):

        ini_time=time.time()

        mmg_parameters = KratosMultiphysics.Parameters("""
        {
            "discretization_type"              : "STANDARD",
            "save_external_files"              : false,
            "initialize_entities"              : false,
            "preserve_flags"                   : false,
            "save_mdpa_file"                       : false,
            "echo_level"                       : 3
        }
        """)

        mmg_process =MeshingApplication.MmgProcess3D(self.body_model_part.GetRootModelPart(), mmg_parameters)
        mmg_process.Execute()

        KratosMultiphysics.Logger.PrintInfo('DefineWakeProcess','Remesh time: ',time.time()-ini_time)