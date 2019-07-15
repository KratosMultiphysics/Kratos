import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp
import math


def DotProduct(A,B):
    result = 0
    for i,j in zip(A,B):
        result += i*j
    return result

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
            "wing_tips_model_part_name": "",
            "wake_stl_file_name" : "",
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

        wing_tips_model_part_name = settings["wing_tips_model_part_name"].GetString()
        if wing_tips_model_part_name == "":
            err_msg = "Empty model_part_name in DefineWakeProcess3D\n"
            err_msg += "Please specify the model part that contains the wing tips nodes"
            raise Exception(err_msg)
        self.wing_tips_model_part = Model[wing_tips_model_part_name]

        self.wake_stl_file_name = settings["wake_stl_file_name"].GetString()
        if self.wake_stl_file_name == "":
            err_msg = "Empty wake_stl_file_name in DefineWakeProcess3D\n"
            err_msg += "Please specify the stl file name that contains the wake surface nodes"
            raise Exception(err_msg)

        self.epsilon = settings["epsilon"].GetDouble()

        self.fluid_model_part = self.trailing_edge_model_part.GetRootModelPart()

    def ExecuteInitialize(self):

        self.__SetWakeDirectionAndNormal()
        self.__MarkTrailingEdgeAndWingTipNodes()
        self.__CreateWakeModelPart()
        self.__MarkWakeElements()
        # Mark the elements touching the trailing edge from below as kutta
        self.__MarkKuttaElements()
        self.__VisualizeWake()

    def __SetWakeDirectionAndNormal(self):
        free_stream_velocity = self.fluid_model_part.ProcessInfo.GetValue(CPFApp.FREE_STREAM_VELOCITY)
        if(free_stream_velocity.Size() != 3):
            raise Exception('The free stream velocity should be a vector with 3 components!')
        self.wake_direction = KratosMultiphysics.Vector(3)
        vnorm = math.sqrt(
            free_stream_velocity[0]**2 + free_stream_velocity[1]**2 + free_stream_velocity[2]**2)
        self.wake_direction[0] = free_stream_velocity[0]/vnorm
        self.wake_direction[1] = free_stream_velocity[1]/vnorm
        self.wake_direction[2] = free_stream_velocity[2]/vnorm

        # For now plane wake surfaces are considered.
        # TODO: Generalize this to curved wake surfaces
        self.wake_normal = KratosMultiphysics.Vector(3)
        self.wake_normal[0] = 0.0
        self.wake_normal[1] = 0.0
        self.wake_normal[2] = 1.0

        self.span_direction = KratosMultiphysics.Vector(3)
        self.span_direction[0] = self.wake_normal[1] * self.wake_direction[2] - self.wake_normal[2] * self.wake_direction[1]
        self.span_direction[1] = self.wake_normal[2] * self.wake_direction[0] - self.wake_normal[0] * self.wake_direction[2]
        self.span_direction[0] = self.wake_normal[0] * self.wake_direction[1] - self.wake_normal[1] * self.wake_direction[0]

    def __MarkTrailingEdgeAndWingTipNodes(self):
        # Mark trailing edge nodes
        for node in self.trailing_edge_model_part.Nodes:
            node.SetValue(CPFApp.TRAILING_EDGE, True)

        # Mark wing tip nodes
        for node in self.wing_tips_model_part.Nodes:
            node.SetValue(CPFApp.WING_TIP, True)
            # Compute the local span direction for later
            # TODO: Generalize this for more than one lifting surface
            for other_node in self.wing_tips_model_part.Nodes:
                wing_span_direction = other_node - node
                if(DotProduct(wing_span_direction,wing_span_direction) > self.epsilon):
                    node.SetValue(CPFApp.WING_SPAN_DIRECTION, wing_span_direction)

    # This function imports the stl file containing the wake and creates the wake model part out of it.
    # TODO: implement an automatic generation of the wake
    def __CreateWakeModelPart(self):
        from stl import mesh #this requires numpy-stl
        wake_stl_mesh = mesh.Mesh.from_multi_file(self.wake_stl_file_name)
        self.wake_model_part = self.model.CreateModelPart("wake_model_part")

        dummy_property = self.wake_model_part.Properties[0]
        node_id = 1
        elem_id = 1

        # Looping over stl meshes
        for stl_mesh in wake_stl_mesh:
            for vertex in stl_mesh.points:
                node1 = self.wake_model_part.CreateNewNode(node_id, float(vertex[0]), float(vertex[1]), float(vertex[2]))
                node_id+=1
                node2 = self.wake_model_part.CreateNewNode(node_id, float(vertex[3]), float(vertex[4]), float(vertex[5]))
                node_id+=1
                node3 = self.wake_model_part.CreateNewNode(node_id, float(vertex[6]), float(vertex[7]), float(vertex[8]))
                node_id+=1

                self.wake_model_part.CreateNewElement("Element3D3N", elem_id,  [
                                              node1.Id, node2.Id, node3.Id], dummy_property)
                elem_id += 1

    # Check which elements are cut and mark them as wake
    def __MarkWakeElements(self):
        KratosMultiphysics.Logger.PrintInfo('...Selecting wake elements...')
        distance_calculator = KratosMultiphysics.CalculateSignedDistanceTo3DSkinProcess(
            self.wake_model_part, self.fluid_model_part)
        distance_calculator.Execute()

        if(self.fluid_model_part.HasSubModelPart("trailing_edge_model_part")):
            for elem in self.trailing_edge_model_part.Elements:
                elem.Set(KratosMultiphysics.TO_ERASE)
            self.trailing_edge_model_part.RemoveElements(KratosMultiphysics.TO_ERASE)
        else:
            self.trailing_edge_model_part = self.fluid_model_part.CreateSubModelPart("trailing_edge_model_part")
        #List to store trailing edge elements id
        self.trailing_edge_element_id_list = []

        for elem in self.fluid_model_part.Elements:
            # Mark and save the elements touching the trailing edge
            self.__MarkTrailingEdgeElement(elem)
            if(elem.Is(KratosMultiphysics.TO_SPLIT)):
                #elem.SetValue(CPFApp.WAKE, True)
                elem.SetValue(CPFApp.DECOUPLED_TRAILING_EDGE_ELEMENT, True)
                wake_elemental_distances = elem.GetValue(KratosMultiphysics.ELEMENTAL_DISTANCES)
                for i in range(len(wake_elemental_distances)):
                    if(abs(wake_elemental_distances[i]) < self.epsilon ):
                        if(wake_elemental_distances[i] < 0.0):
                            wake_elemental_distances[i] = -self.epsilon
                        else:
                            wake_elemental_distances[i] = self.epsilon
                elem.SetValue(CPFApp.WAKE_ELEMENTAL_DISTANCES,wake_elemental_distances)
                counter = 0
                for node in elem.GetNodes():
                    node.SetSolutionStepValue(CPFApp.WAKE_DISTANCE, wake_elemental_distances[counter])
                    if(wake_elemental_distances[counter] > 0.0):
                        node.SetValue(KratosMultiphysics.WATER_PRESSURE, 1.0)
                    else:
                        node.SetValue(KratosMultiphysics.WATER_PRESSURE, -1.0)#
                    counter +=1

        self.__SaveTrailingEdgeElements()
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
        # This function stores the trailing edge element
        # to its submodelpart.
        self.trailing_edge_model_part.AddElements(self.trailing_edge_element_id_list)

    def __MarkKuttaElements(self):
        # This function selects the kutta elements. Kutta elements
        # are touching the trailing edge from below.
        for elem in self.trailing_edge_model_part.Elements:
            # Check if the element is touching the wing tip
            wing_tip, trailing_edge_node, number_of_non_te_nodes = self.__CheckIfWingTip(elem)
            if wing_tip:
                # Elements touching the wing tip are set to normal
                elem.SetValue(CPFApp.WAKE, False)
                elem.SetValue(CPFApp.DECOUPLED_TRAILING_EDGE_ELEMENT, False)
            else:
                distances_to_te = self.__ComputeDistancesToTrailingEdgeNode(elem, trailing_edge_node, number_of_non_te_nodes)
                self.__CheckIfKuttaElement(elem, distances_to_te, number_of_non_te_nodes)

    def __CheckIfWingTip(self, elem):
        wing_tip = False
        number_of_te_nodes = 0

        for elnode in elem.GetNodes():
            if(elnode.GetValue(CPFApp.TRAILING_EDGE)):
                    trailing_edge_node = elnode
                    number_of_te_nodes += 1
            if(elnode.GetValue(CPFApp.WING_TIP)):
                wing_tip = True
                wing_tip_node = elnode
                wing_span_direction = elnode.GetValue(CPFApp.WING_SPAN_DIRECTION)

        if(number_of_te_nodes > 1):
            wing_tip = False
        number_of_non_te_nodes = 4 - number_of_te_nodes

        if(wing_tip):
            for elnode in elem.GetNodes():
                if not (elnode.GetValue(CPFApp.WING_TIP)):
                    distance = elnode - wing_tip_node
                    wing_span_projection = DotProduct(distance, wing_span_direction)
                    free_stream_projection = DotProduct(distance, self.wake_direction)
                    if(wing_span_projection > 0.0 and free_stream_projection > 0.0):
                        wing_tip = False

        return wing_tip, trailing_edge_node, number_of_non_te_nodes

    def __ComputeDistancesToTrailingEdgeNode(self, elem, trailing_edge_node, number_of_non_te_nodes):
        # This function computes the distance of the element nodes
        # to the wake
        nodal_distances_to_te = KratosMultiphysics.Vector(number_of_non_te_nodes)
        counter = 0
        for elnode in elem.GetNodes():
            # Looping only over non traling edge nodes
            if not (elnode.GetValue(CPFApp.TRAILING_EDGE)):
                # Compute the distance from the node to the trailing edge
                distance = elnode - trailing_edge_node

                # Compute the projection of the distance vector in the wake normal direction
                distance_to_wake = DotProduct(distance, self.wake_normal)

                # Nodes laying on the wake have a negative distance
                if(abs(distance_to_wake) < self.epsilon):
                    distance_to_wake = -self.epsilon

                nodal_distances_to_te[counter] = distance_to_wake
                counter += 1

        return nodal_distances_to_te

    def __CheckIfKuttaElement(self, elem, distances_to_te, number_of_non_te_nodes):
        # This function checks whether the element is kutta
        # Initialize counters
        number_of_nodes_with_positive_distance = 0
        number_of_nodes_with_negative_distance = 0

        # Count how many element nodes are above and below the wake
        for nodal_distance_to_wake in distances_to_te:
            if(nodal_distance_to_wake < 0.0):
                number_of_nodes_with_negative_distance += 1
            else:
                number_of_nodes_with_positive_distance += 1

        if(number_of_nodes_with_negative_distance > 0 and number_of_nodes_with_positive_distance > 0):
            if(elem.GetValue(CPFApp.DECOUPLED_TRAILING_EDGE_ELEMENT)):
                # Elements with nodes above and below the wake and that were already marked as wake, are wake elements
                #elem.SetValue(CPFApp.WAKE, True)

                elem.SetValue(CPFApp.DECOUPLED_TRAILING_EDGE_ELEMENT, True)
        elif(number_of_nodes_with_negative_distance > number_of_non_te_nodes - 1):
            # Elements with all non trailing edge nodes below the wake are kutta
            elem.SetValue(CPFApp.KUTTA, True)
        else:
            # Elements with all non trailing edge nodes above the wake are normal
            elem.SetValue(CPFApp.WAKE, False)
            elem.SetValue(CPFApp.DECOUPLED_TRAILING_EDGE_ELEMENT, False)
            elem.SetValue(CPFApp.ZERO_VELOCITY_CONDITION, True)

    def __VisualizeWake(self):


        for elem in self.trailing_edge_model_part.Elements:
            if(elem.GetValue(CPFApp.DECOUPLED_TRAILING_EDGE_ELEMENT)):
                print(elem.Id)
                pass
            elif(elem.GetValue(CPFApp.KUTTA)):
                #print(elem.Id)
                pass
            else:
                #print(elem.Id)
                pass
                if not (elem.GetValue(CPFApp.ZERO_VELOCITY_CONDITION)):
                    print(elem.Id)
                    pass




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

        from gid_output_process import GiDOutputProcess
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
    #'''
