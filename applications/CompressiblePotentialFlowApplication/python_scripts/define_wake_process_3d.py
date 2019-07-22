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
            "body_model_part_name": "",
            "wake_stl_file_name" : "",
            "wake_normal": [0.0,0.0,1.0],
            "output_wake": true,
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

    def ExecuteInitialize(self):

        self.__SetWakeAndSpanDirections()
        # Save the trailing edge and wing tip nodes for further computations
        self.__MarkTrailingEdgeAndWingTipNodes()
        # Save the lower surface normals to help mark kutta elements later on
        self.__ComputeLowerSurfaceNormals()
        # Read wake from stl and create the wake model part
        self.__CreateWakeModelPart()
        # Check which elements are cut and mark them as wake
        self.__MarkWakeElements()
        # Mark the elements touching the trailing edge from below as kutta
        self.__MarkKuttaElements()
        # Output the wake in GiD for visualization
        if(self.output_wake):
            self.__VisualizeWake()
        # Output element ids in the terminal
        self.__TerminalPrint()

        # # # Mark wing tip nodes
        # for node in self.wing_tips_model_part.Nodes:
        #     node.SetValue(CPFApp.TRAILING_EDGE, False)

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

    def __MarkTrailingEdgeAndWingTipNodes(self):
        # Mark trailing edge nodes
        for node in self.trailing_edge_model_part.Nodes:
            node.SetValue(CPFApp.TRAILING_EDGE, True)
            # Marking only one node
            # if(node.Id == 56):
            #     node.SetValue(CPFApp.ZERO_VELOCITY_CONDITION, True)

        # Mark wing tip nodes
        for node in self.wing_tips_model_part.Nodes:
            node.SetValue(CPFApp.WING_TIP, True)
            # Compute the local span direction for later
            # TODO: Generalize this for more than one lifting surface
            for other_node in self.wing_tips_model_part.Nodes:
                wing_span_direction = other_node - node
                if(DotProduct(wing_span_direction,wing_span_direction) > self.epsilon):
                    node.SetValue(CPFApp.WING_SPAN_DIRECTION, wing_span_direction)

        # Selecting only one wing tip
        counter = 0
        for node in self.wing_tips_model_part.Nodes:
            if counter < 1:
                node.SetValue(CPFApp.WING_TIP, False)
                print(node.Id)
                counter += 1
            else:
                print('aaaa')
                print(node.Id)
                counter += 1

    def __ComputeLowerSurfaceNormals(self):
        for cond in self.body_model_part.Conditions:
            # The surface normal points outisde the domain
            surface_normal = cond.GetGeometry().Normal()
            projection = DotProduct(surface_normal, self.wake_normal)
            # The surface normal in the same direction as the wake normal belongs to the lower_surface
            if(projection > 0.0):
                for node in cond.GetNodes():
                    node.SetValue(KratosMultiphysics.NORMAL,surface_normal)

    # This function imports the stl file containing the wake and creates the wake model part out of it.
    # TODO: implement an automatic generation of the wake
    def __CreateWakeModelPart(self):
        from stl import mesh #this requires numpy-stl
        wake_stl_mesh = mesh.Mesh.from_multi_file(self.wake_stl_file_name)
        self.wake_model_part = self.model.CreateModelPart("wake_model_part")

        dummy_property = self.wake_model_part.Properties[0]
        node_id = 1
        elem_id = 1
        # Vertical translation of the wake
        z =  0#self.epsilon

        # Looping over stl meshes
        for stl_mesh in wake_stl_mesh:
            for vertex in stl_mesh.points:
                node1 = self.wake_model_part.CreateNewNode(node_id, float(vertex[0]), float(vertex[1]), float(vertex[2])+z)
                node_id+=1
                node2 = self.wake_model_part.CreateNewNode(node_id, float(vertex[3]), float(vertex[4]), float(vertex[5])+z)
                node_id+=1
                node3 = self.wake_model_part.CreateNewNode(node_id, float(vertex[6]), float(vertex[7]), float(vertex[8])+z)
                node_id+=1

                self.wake_model_part.CreateNewElement("Element3D3N", elem_id,  [
                                              node1.Id, node2.Id, node3.Id], dummy_property)
                elem_id += 1

    # Check which elements are cut and mark them as wake
    def __MarkWakeElements(self):
        KratosMultiphysics.Logger.PrintInfo('...Selecting wake elements...')

        # Mark cut elements and compute elemental distances
        # Attention: Note that in this process a negative distance is assigned to nodes
        # laying on the wake. In 2D it is done viceversa.
        use_discontinuous_distance_process = False
        if(use_discontinuous_distance_process):
            distance_calculator = KratosMultiphysics.CalculateDiscontinuousDistanceToSkinProcess3D(
            self.fluid_model_part, self.wake_model_part)
        else:
            distance_calculator = KratosMultiphysics.CalculateDistanceToSkinProcess3D(
            self.fluid_model_part, self.wake_model_part)
        distance_calculator.Execute()

        #List to store trailing edge elements id
        self.trailing_edge_element_id_list = []

        if not self.fluid_model_part.HasSubModelPart("wake__elements"):
            self.wake_sub_model_part = self.fluid_model_part.CreateSubModelPart("wake__elements")
        else: self.wake_sub_model_part = self.fluid_model_part.GetSubModelPart("wake__elements")
        self.wake_element_id_list = []

        for node in self.fluid_model_part.Nodes:
            node.SetValue(CPFApp.WAKE_DISTANCE, 100)
            #node.SetValue(KratosMultiphysics.WATER_PRESSURE, 100)

        for elem in self.fluid_model_part.Elements:
            # Mark and save the elements touching the trailing edge
            self.__MarkTrailingEdgeElement(elem)


            # wake_elemental_distances = elem.GetValue(KratosMultiphysics.ELEMENTAL_DISTANCES)
            # # Save wake nodal distances
            # counter = 0
            # for node in elem.GetNodes():
            #     node.SetValue(CPFApp.WAKE_DISTANCE,wake_elemental_distances[counter])
            #     counter += 1

            # counter = 0
            # npos = 0
            # nneg = 0
            # for counter in range(0,4):
            #     if wake_elemental_distances[counter] > 0.0:
            #         npos += 1
            #     else:
            #         nneg += 1

            # if (npos>0 and nneg >0):

            # Cut elements are wake
            if(elem.Is(KratosMultiphysics.TO_SPLIT)):
                # Mark wake element
                elem.SetValue(CPFApp.WAKE, True)
                self.wake_element_id_list.append(elem.Id)
                # Save wake elemental distances
                wake_elemental_distances = elem.GetValue(KratosMultiphysics.ELEMENTAL_DISTANCES)
                # if(elem.Id==78):
                #     print(wake_elemental_distances)
                # # Save wake nodal distances
                # counter = 0
                # for node in elem.GetNodes():
                #     node.SetValue(CPFApp.WAKE_DISTANCE,wake_elemental_distances[counter])
                #     counter += 1
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
                # Mark nodes above and below the wake with WATER_PRESSURE variable for visualization
                counter = 0
                for node in elem.GetNodes():
                    if(wake_elemental_distances[counter] > 0.0):
                        node.SetValue(KratosMultiphysics.WATER_PRESSURE, 1.0)
                    else:
                        node.SetValue(KratosMultiphysics.WATER_PRESSURE, -1.0)#
                    counter +=1

        self.wake_sub_model_part.AddElements(self.wake_element_id_list)
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
        # This function stores the trailing edge elements
        # in its submodelpart.
        if(self.fluid_model_part.HasSubModelPart("trailing_edge_model_part")):
            for elem in self.trailing_edge_model_part.Elements:
                elem.Set(KratosMultiphysics.TO_ERASE)
            self.trailing_edge_model_part.RemoveElements(KratosMultiphysics.TO_ERASE)
        else:
            self.trailing_edge_model_part = self.fluid_model_part.CreateSubModelPart("trailing_edge_model_part")
        self.trailing_edge_model_part.AddElements(self.trailing_edge_element_id_list)

    def __MarkKuttaElements(self):
        # This function selects the kutta elements. Kutta elements
        # are touching the trailing edge from below.

        # Loop over elements touching the trailing edge
        for elem in self.trailing_edge_model_part.Elements:
            # Check if it is a wing tip element
            wing_tip = self.__CheckIfWingTipElement(elem)
            if wing_tip:
                # Wing tip elements are set to normal
                # TODO: Check what to do with wing tip elements
                # and tip vortice elements in general
                elem.SetValue(CPFApp.WAKE, False)
                self.wake_sub_model_part.RemoveElement(elem)
            else:
                trailing_edge_node, number_of_non_te_nodes = self.__GetATrailingEdgeNodeAndNumberOfNonTENodes(elem)
                nodal_distances = self.__ComputeNodalDistancesToWakeAndLowerSurface(elem, trailing_edge_node, number_of_non_te_nodes)
                self.__CheckIfKuttaElement(elem, nodal_distances, number_of_non_te_nodes)

    def __CheckIfWingTipElement(self, elem):
        # Wing tip elements are elements with one node at the wing tip
        # and with the rest of the nodes on the side of the wing and wake.
        wing_tip = False

        # Checking if the element has a node at the wing tip
        for elnode in elem.GetNodes():
            if(elnode.GetValue(CPFApp.WING_TIP)):
                wing_tip = True
                wing_tip_node = elnode
                wing_span_direction = elnode.GetValue(CPFApp.WING_SPAN_DIRECTION)

        if(wing_tip):
            # Setting only elements touching the wing tip as structure
            elem.Set(KratosMultiphysics.STRUCTURE)
            for elnode in elem.GetNodes():
                if not (elnode.GetValue(CPFApp.WING_TIP)):
                    # Checking if the rest of the nodes are on the side
                    distance = elnode - wing_tip_node
                    wing_span_projection = DotProduct(distance, wing_span_direction)
                    # A positive wing_span_projection means that the node is not on the side
                    # but actually laying somewhere above or below the wing and/or wake
                    if(wing_span_projection > 0.0):
                        wing_tip = False

        return wing_tip

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

                # Node laying either above or below the lower surface
                if(free_stream_direction_distance < 0.0):
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
                elnode.SetValue(CPFApp.WAKE_DISTANCE,distance)
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
            #elem.Set(KratosMultiphysics.STRUCTURE)
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

    def __TerminalPrint(self):
        # Print trailing_edge_model_part elements
        for elem in self.trailing_edge_model_part.Elements:
            if(elem.GetValue(CPFApp.WAKE)):
                #print(elem.Id)
                pass
            elif(elem.GetValue(CPFApp.KUTTA)):
                #print(elem.Id)
                pass
            else:
                #print(elem.Id)
                pass

        # Print fluid_model_part elements
        counter_structure = 0
        counter_wake = 0
        for elem in self.fluid_model_part.Elements:
            #print(elem.Id)
            if(elem.Is(KratosMultiphysics.STRUCTURE) and elem.GetValue(CPFApp.WAKE)):
                #print('aa')
                #print(elem.Id)
                counter_structure += 1
                pass
            if(elem.GetValue(CPFApp.WAKE)):
                #print(elem.Id)
                counter_wake += 1
                pass
            elif(elem.GetValue(CPFApp.KUTTA)):
                #print(elem.Id)
                pass
            else:
                #print(elem.Id)
                pass

        # is_te = KratosMultiphysics.Vector(4)
        # for i in range(len(is_te)):
        #     is_te[i] = 0
        # for elem in self.wake_sub_model_part.Elements:
        #     number_of_te_nodes = 0
        #     for elnode in elem.GetNodes():
        #         if(elnode.GetValue(CPFApp.TRAILING_EDGE)):
        #             number_of_te_nodes += 1
        #     if(number_of_te_nodes > 1):
        #         counter = 0
        #         to_be_marked = True
        #         for elnode in elem.GetNodes():
        #             if(elnode.GetValue(CPFApp.TRAILING_EDGE) and to_be_marked):
        #                 is_te[counter] = True
        #                 counter += 1
        #                 to_be_marked = True
        #                 #elnode.SetValue(CPFApp.TRAILING_EDGE, False)

        # # Selecting only one wing tip
        # counter = 0
        # for node in self.wing_tips_model_part.Nodes:
        #     if counter < 1:
        #         node.SetValue(CPFApp.WING_TIP, False)
        #         counter += 1
        #         pass

        for node in self.fluid_model_part.Nodes:
            if(node.Id == 1151):
                #node.SetValue(CPFApp.WING_TIP, True)
                pass



        counter_wake_elements = 0
        for elem in self.wake_sub_model_part.Elements:
            #print(elem.Id)
            counter_wake_elements += 1

        print('counter_wake = ', counter_wake)
        print('counter_structure = ', counter_structure)
        print('counter_wake_elements = ', counter_wake_elements)

    def ExecuteFinalize(self):
        CPFApp.CheckWakeConditionProcess3D(self.wake_sub_model_part, 1e-13, 1).Execute()
