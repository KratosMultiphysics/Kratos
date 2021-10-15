import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp
import time
import math


def Factory(settings, Model):
    if( not isinstance(settings,KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    return DefineEmbeddedWakeProcess(Model, settings["Parameters"])


class DefineEmbeddedWakeProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings):
        # Call the base Kratos process constructor
        KratosMultiphysics.Process.__init__(self)

        # Check default settings
        default_settings = KratosMultiphysics.Parameters(r'''{
            "model_part_name": "",
            "epsilon": 1e-9
        }''')
        settings.ValidateAndAssignDefaults(default_settings)

        self.main_model_part = Model[settings["model_part_name"].GetString()].GetRootModelPart()
        if Model.HasModelPart("wake"):
            Model.DeleteModelPart("wake")
        self.wake_model_part=Model.CreateModelPart("wake")
        self.model=Model

        self.epsilon = settings["epsilon"].GetDouble()

    # def ExecuteInitialize(self):
    def ExecuteInitializeSolutionStep(self):
        KratosMultiphysics.FindGlobalNodalElementalNeighboursProcess(self.main_model_part).Execute()

        ini_time = time.time()


        # self.list_of_failed_te_nodes =[]
        # max_x_coordinate = -1e30
        # restarted_search_maximum = 1e30
        # trailing_edge_candidate = self.FindNode(max_x_coordinate, restarted_search_maximum)
        # is_valid = self.CheckIfValid(trailing_edge_candidate)
        # while(not is_valid):
        #     trailing_edge_candidate.SetSolutionStepValue(CPFApp.GEOMETRY_DISTANCE, 1e-9)
        #     self.list_of_failed_te_nodes.append(trailing_edge_candidate.Id)
        #     trailing_edge_candidate = self.FindNode(max_x_coordinate, trailing_edge_candidate.X)
        #     is_valid = self.CheckIfValid(trailing_edge_candidate)
        #     if not is_valid:
        #         self.list_of_failed_te_nodes.append(trailing_edge_candidate.Id)
        #         trailing_edge_candidate.SetSolutionStepValue(CPFApp.GEOMETRY_DISTANCE, 1e-9)
        # stop
        self._DefineWakeModelPart()
        # self._MoveAndRotateWake()
        # Executing define wake process

        cpp_ini_time = time.time()
        CPFApp.DefineEmbeddedWakeProcess(self.main_model_part, self.wake_model_part).Execute()
        KratosMultiphysics.Logger.PrintInfo('EmbeddedWake','CPP Wake computation time: ',time.
        time()-cpp_ini_time)
        # angle = self.main_model_part.ProcessInfo.GetValue(CPFApp.ROTATION_ANGLE)

        # for node in self.main_model_part.GetSubModelPart("trailing_edge_sub_model_part").Nodes:
        #     if node.GetValue(CPFApp.AIRFOIL):
        #         self.structure_node = node

        # self._RedefineWake()

        # self.main_model_part.GetElement(144354).Set(KratosMultiphysics.STRUCTURE, True)
        # self.main_model_part.GetElement(73449).SetValue(KratosMultiphysics.CompressiblePotentialFlowApplication.WAKE, True)
        # self.main_model_part.GetElement(83215).Set(KratosMultiphysics.ACTIVE, False)
        # self.main_model_part.GetElement(97893).Set(KratosMultiphysics.ACTIVE, False)
        # self.main_model_part.GetElement(83215).Set(KratosMultiphysics.STRUCTURE, False)
        # self.main_model_part.GetElement(97893).Set(KratosMultiphysics.STRUCTURE, False)
        # for node in self.main_model_part.GetElement(83215).GetNodes():
        #     node.SetValue(KratosMultiphysics.CompressiblePotentialFlowApplication.WING_TIP, False)
        # for node in self.main_model_part.GetElement(97893).GetNodes():
        #     node.SetValue(KratosMultiphysics.CompressiblePotentialFlowApplication.WING_TIP, False)
        # self.main_model_part.GetNode(50251).SetValue(KratosMultiphysics.CompressiblePotentialFlowApplication.TRAILING_EDGE, True)
        # self.main_model_part.GetElement(83215).SetValue(KratosMultiphysics.CompressiblePotentialFlowApplication.WAKE, False)
        # self.main_model_part.GetElement(97893).SetValue(KratosMultiphysics.CompressiblePotentialFlowApplication.WAKE, False)
        # self.main_model_part.GetElement(77001).SetValue(KratosMultiphysics.CompressiblePotentialFlowApplication.KUTTA, False)
        # self.main_model_part.GetElement(123199).SetValue(KratosMultiphysics.CompressiblePotentialFlowApplication.KUTTA, False)

        # self.main_model_part.GetElement(80156).SetValue(KratosMultiphysics.CompressiblePotentialFlowApplication.KUTTA, True)
        # self.main_model_part.GetElement(106586).SetValue(KratosMultiphysics.CompressiblePotentialFlowApplication.KUTTA, True)
        # self.main_model_part.GetElement(126747).SetValue(KratosMultiphysics.CompressiblePotentialFlowApplication.KUTTA, True)


        # self.main_model_part.GetElement(51798).Set(KratosMultiphysics.STRUCTURE, True)
        # self.main_model_part.GetNode(13407).SetValue(KratosMultiphysics.CompressiblePotentialFlowApplication.WING_TIP, True)
        # self.main_model_part.GetNode(37316).SetValue(KratosMultiphysics.CompressiblePotentialFlowApplication.WING_TIP, True)
        # self.main_model_part.GetNode(49212).SetValue(KratosMultiphysics.CompressiblePotentialFlowApplication.WING_TIP, True)
        # self.main_model_part.GetNode(53939).SetValue(KratosMultiphysics.CompressiblePotentialFlowApplication.WING_TIP, True)
        # self.main_model_part.GetNode(22240).SetValue(KratosMultiphysics.CompressiblePotentialFlowApplication.WING_TIP, True)
        # self.main_model_part.GetNode(40990).SetValue(KratosMultiphysics.CompressiblePotentialFlowApplication.WING_TIP, True)
        # self.main_model_part.GetNode(15936).SetValue(KratosMultiphysics.CompressiblePotentialFlowApplication.WING_TIP, False)
        # self.main_model_part.GetNode(18238).SetValue(KratosMultiphysics.CompressiblePotentialFlowApplication.WING_TIP, False)
        # self.main_model_part.GetNode(22386).SetValue(KratosMultiphysics.CompressiblePotentialFlowApplication.WING_TIP, False)
        # self.main_model_part.GetNode(22790).SetValue(KratosMultiphysics.CompressiblePotentialFlowApplication.WING_TIP, False)
        # self.main_model_part.GetNode(42142).SetValue(KratosMultiphysics.CompressiblePotentialFlowApplication.WING_TIP, False)
        # self.main_model_part.GetNode(44935).SetValue(KratosMultiphysics.CompressiblePotentialFlowApplication.WING_TIP, False)
        # self.main_model_part.GetNode(39287).SetValue(KratosMultiphysics.CompressiblePotentialFlowApplication.WING_TIP, True)
        # self.main_model_part.GetNode(42146).SetValue(KratosMultiphysics.CompressiblePotentialFlowApplication.WING_TIP, True)
        # self.main_model_part.GetNode(984).SetValue(KratosMultiphysics.CompressiblePotentialFlowApplication.WING_TIP, False)
        # self.main_model_part.GetNode(2654).SetValue(KratosMultiphysics.CompressiblePotentialFlowApplication.WING_TIP, False)


        # self.main_model_part.GetNode(64620).SetValue(KratosMultiphysics.CompressiblePotentialFlowApplication.WING_TIP, True)
        # self.main_model_part.GetNode(69404).SetValue(KratosMultiphysics.CompressiblePotentialFlowApplication.WING_TIP, True)





        KratosMultiphysics.Logger.PrintInfo('EmbeddedWake','Wake computation time: ',time.time()-ini_time)
    def _RedefineWake(self):
        ini_time = time.time()

        # lower_surface_eslem_ids=[]
        # upper_surface_elem_ids=[]
        # for elem in self.main_model_part.Elements:
        #     if elem.IsNot(KratosMultiphysics.ACTIVE):
        #         if max_inactive_x < elem.GetGeometry().Center().X:
        #             max_inactive_x = elem.GetGeometry().Center().X
        #     if elem.Is(KratosMultiphysics.TO_SPLIT) and elem.Is(KratosMultiphysics.ACTIVE):
        #         pos_nodes=[]
        #         neg_nodes=[]
        #         for node in elem.GetNodes():
        #             distance = node.GetSolutionStepValue(CPFApp.GEOMETRY_DISTANCE)
        #             # print(elem.Id, distance)
        #             if distance > 0.0:
        #                 pos_nodes.append(node)
        #             else:
        #                 neg_nodes.append(node)
        #         max_distance=0
        #         final_distance = 0
        #         for pos_node in pos_nodes:
        #             for neg_node in neg_nodes:
        #                 node_distance = pos_node.Y - neg_node.Y
        #                 if abs(node_distance) > max_distance:
        #                     max_distance = abs(node_distance)
        #                     final_distance = node_distance

        #         if final_distance > 0.0:
        #             elem.SetValue(CPFApp.KUTTA, False)
        #             for node in elem.GetNodes():
        #                 node.SetValue(CPFApp.UPPER_SURFACE, True)
        #         else:
        #             for node in elem.GetNodes():
        #                 node.SetValue(CPFApp.LOWER_SURFACE, True)
        # for elem in self.main_model_part.Elements:
        #     counter_lower = 0
        #     for node in elem.GetNodes():
        #         if node.GetValue(CPFApp.LOWER_SURFACE):
        #             counter_lower += 1
        #     if counter_lower == 3:
        #         lower_surface_elem_ids.append(elem.Id)
        #     counter_upper = 0
        #     for node in elem.GetNodes():
        #         if node.GetValue(CPFApp.UPPER_SURFACE):
        #             counter_upper += 1
        #     if counter_upper == 3:
        #         upper_surface_elem_ids.append(elem.Id)


        # if self.main_model_part.HasSubModelPart("lower_surface_sub_model_part"):
        #     self.main_model_part.RemoveSubModelPart("lower_surface_sub_model_part")
        # self.lower_surface_element_sub_model_part = self.main_model_part.CreateSubModelPart("lower_surface_sub_model_part")
        # self.lower_surface_element_sub_model_part.AddElements(lower_surface_elem_ids)


        # if self.main_model_part.HasSubModelPart("upper_surface_sub_model_part"):
        #     self.main_model_part.RemoveSubModelPart("upper_surface_sub_model_part")
        # self.upper_surface_element_sub_model_part = self.main_model_part.CreateSubModelPart("upper_surface_sub_model_part")
        # self.upper_surface_element_sub_model_part.AddElements(upper_surface_elem_ids)
        # self.lower_surface_element_sub_model_part = self.main_model_part.GetSubModelPart("lower_surface_sub_model_part")
        # self.upper_surface_element_sub_model_part = self.main_model_part.GetSubModelPart("upper_surface_sub_model_part")

        # for node in self.main_model_part.Nodes:
        # #     # if node.GetValue(CPFApp.UPPER_SURFACE) and node.GetValue(CPFApp.LOWER_SURFACE) and node.X > max_inactive_x-0.1 and node.GetSolutionStepValue(CPFApp.GEOMETRY_DISTANCE) < 0.0:
        #     if node.GetValue(CPFApp.UPPER_SURFACE) and node.GetValue(CPFApp.LOWER_SURFACE) and node.X > self.max_inactive_x-0.5 and node.X < self.max_inactive_x and node.GetSolutionStepValue(CPFApp.GEOMETRY_DISTANCE) < 0.0:
        # #     # if node.GetValue(CPFApp.UPPER_SURFACE) and node.GetValue(CPFApp.LOWER_SURFACE) and node.X < self.trailing_edge_node.X  and node.X > max_inactive_x-0.1 and node.GetSolutionStepValue(CPFApp.GEOMETRY_DISTANCE) < 0.0:
        #         node.SetValue(CPFApp.TRAILING_EDGE, True)
        # #         print("SETTING NODE TO TRAILING EDGE", node.Id)
        # stop

        # # for elem in self.main_model_part.Elements:
        # if self.structure_node.GetValue(CPFApp.WAKE_DISTANCE) > 0.0:
        #     sub_model_part_to_mark_as_kutta=self.lower_surface_element_sub_model_part
        # else:
        #     sub_model_part_to_mark_as_kutta=self.upper_surface_element_sub_model_part
        # for elem in sub_model_part_to_mark_as_kutta.Elements:
        #     if not elem.GetValue(CPFApp.WAKE) and not elem.GetValue(CPFApp.KUTTA) and elem.Is(KratosMultiphysics.ACTIVE):
        #         for node in elem.GetNodes():
        #             if node.GetValue(CPFApp.TRAILING_EDGE) and node.GetSolutionStepValue(CPFApp.GEOMETRY_DISTANCE) < 0.0:
        #                 elem.SetValue(CPFApp.KUTTA, True)
        #                 # pass

        # for node in self.main_model_part.Nodes:
        #     geometry_distance = node.GetSolutionStepValue(CPFApp.GEOMETRY_DISTANCE)
        #     node.SetValue(KratosMultiphysics.TEMPERATURE, geometry_distance)
        #     if node.GetValue(CPFApp.TRAILING_EDGE) and geometry_distance< 0.0:
        #         node.Set(KratosMultiphysics.INSIDE) # is negative?
        #         node.SetSolutionStepValue(CPFApp.GEOMETRY_DISTANCE, -1e30)
        #         node.SetValue(KratosMultiphysics.TEMPERATURE, -1e30)
        # for elem in self.main_model_part.Elements:
        #     is_trailing_edge = False
        #     for node in elem.GetNodes():
        #         if node.GetValue(CPFApp.TRAILING_EDGE):
        #             is_trailing_edge = True
        #     if is_trailing_edge:
        #         elem_distances = elem.GetValue(KratosMultiphysics.ELEMENTAL_DISTANCES)
        #         if elem.GetValue(CPFApp.KUTTA):
        #             for node, elem_distance in zip(elem.GetNodes(), elem_distances):
        #                 if node.GetValue(CPFApp.TRAILING_EDGE):
        #                     old_distance = node.GetValue(KratosMultiphysics.TEMPERATURE)
        #                     # if old_distance*elem_distance > 0.0:
        #                     #     if old_distance*old_distance < 0.0 or elem_distance*elem_distance < 0.0:
        #                     #         print("WARNING, DIFFERENT SIGNS:",elem.Id, node.Id, elem_distance, old_distance)
        #                     if old_distance < 0.0:
        #                         new_distance = max(-abs(old_distance),-abs(elem_distance))
        #                     else:
        #                         new_distance = min(abs(old_distance),abs(elem_distance))
        #                     node.SetValue(KratosMultiphysics.TEMPERATURE, new_distance)
        #                     print(elem.Id, node.Id, elem_distance, node.Is(KratosMultiphysics.INSIDE), new_distance)
        #         else:
        #             for node, elem_distance in zip(elem.GetNodes(), elem_distances):
        #                 if node.GetValue(CPFApp.TRAILING_EDGE):
        #                     old_distance = node.GetSolutionStepValue(CPFApp.GEOMETRY_DISTANCE)
        #                     if old_distance < 0.0:
        #                         new_distance = max(-abs(old_distance),-abs(elem_distance))
        #                     else:
        #                         new_distance = min(abs(old_distance),abs(elem_distance))
        #                     node.SetSolutionStepValue(CPFApp.GEOMETRY_DISTANCE, new_distance)
        #                     print(elem.Id, node.Id, elem_distance, node.Is(KratosMultiphysics.INSIDE), new_distance)
        # print("List of failed te nodes", self.list_of_failed_te_nodes)
        # for node_id in self.list_of_failed_te_nodes:
        #     self.main_model_part.GetNode(node_id).SetSolutionStepValue(CPFApp.GEOMETRY_DISTANCE, 1e-9)


        KratosMultiphysics.Logger.PrintInfo('EmbeddedWake','Redefine time: ',time.time()-ini_time)

    def _DefineWakeModelPart(self):
        ''' This function generates the modelpart of the wake. TODO: make end of the domain user-definable.
        '''
        # stop
        # self.wake_model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        # self.wake_model_part.CreateNewNode(2, 200.0, 0.0, 0.0)
        # self.wake_model_part.CreateNewElement("Element2D2N", 1, [1,2], KratosMultiphysics.Properties(0))

        # skin_model_part = self.model["skin"]
        # te_node = -1
        # for node in skin_model_part.GetSubModelPart("TrailingEdgeNode").Nodes:
        #     te_node = node
        # le_node = -1
        # for node in skin_model_part.GetSubModelPart("LeadingEdgeNode").Nodes:
        #     le_node = node
        # te_weight=0.95
        # print("TE NODE", te_node.X, te_node.Y)
        # self.wake_model_part.CreateNewNode(1, (1-te_weight)*le_node.X+te_weight*te_node.X, (1-te_weight)*le_node.Y+te_weight*te_node.Y, 0.0)
        # self.wake_model_part.CreateNewNode(2, te_node.X, te_node.Y, 0.0)
        # self.wake_model_part.CreateNewNode(3, 200.0, te_node.Y, 0.0)
        # self.wake_model_part.CreateNewElement("Element2D2N", 1, [1,2], KratosMultiphysics.Properties(0))
        # self.wake_model_part.CreateNewElement("Element2D2N", 2, [2,3], KratosMultiphysics.Properties(0))
        skin_model_part = self.model["skin"]
        te_node = -1
        for node in skin_model_part.GetSubModelPart("TrailingEdgeNode").Nodes:
            te_node = node
        le_node = -1
        for node in skin_model_part.GetSubModelPart("LeadingEdgeNode").Nodes:
            le_node = node
        te_weight=0.95
        self.wake_model_part.CreateNewNode(1, (1-te_weight)*le_node.X+te_weight*te_node.X, (1-te_weight)*le_node.Y+te_weight*te_node.Y, 0.0)
        self.wake_model_part.CreateNewNode(2, 200.0, te_node.Y, 0.0)
        self.wake_model_part.CreateNewElement("Element2D2N", 1, [1,2], KratosMultiphysics.Properties(0))

        self.moving_parameters = KratosMultiphysics.Parameters()
        self.moving_parameters.AddEmptyValue("rotation_point")
        self.moving_parameters["rotation_point"].SetVector([(1-te_weight)*le_node.X+te_weight*te_node.X, (1-te_weight)*le_node.Y+te_weight*te_node.Y, 0.0])
        self.moving_parameters.AddEmptyValue("rotation_angle")
        angle=math.radians(-self.main_model_part.ProcessInfo.GetValue(CPFApp.ROTATION_ANGLE))
        self.moving_parameters["rotation_angle"].SetDouble(angle)
        CPFApp.MoveModelPartProcess(self.wake_model_part, self.moving_parameters).Execute()
        print(self.wake_model_part.GetNode(1).X, self.wake_model_part.GetNode(1).Y)
        print(self.wake_model_part.GetNode(2).X, self.wake_model_part.GetNode(2).Y)

    def _MoveAndRotateWake(self):
        ''' This function moves and rotates the wake with the same parameters as the geometry.
        '''
        # self.moving_parameters = KratosMultiphysics.Parameters()
        # self.moving_parameters.AddEmptyValue("origin")
        # self.moving_parameters["origin"].SetVector(self.main_model_part.ProcessInfo.GetValue(CPFApp.WAKE_ORIGIN))
        # self.moving_parameters.AddEmptyValue("rotation_angle")
        # angle=math.radians(-self.main_model_part.ProcessInfo.GetValue(CPFApp.ROTATION_ANGLE))
        # self.moving_parameters["rotation_angle"].SetDouble(angle)
        # CPFApp.MoveModelPartProcess(self.wake_model_part, self.moving_parameters).Execute()
        # print(self.main_model_part.ProcessInfo.GetValue(CPFApp.WAKE_ORIGIN))
        # print(angle)

        # inital_point = self.main_model_part.ProcessInfo.GetValue(CPFApp.WAKE_ORIGIN)
        # vector = KratosMultiphysics.Vector(3, 0.0)
        # vector[0] =  math.cos(angle)
        # vector[1] = -math.sin(-angle)

        # max_projection = -1e30
        # for element in self.main_model_part.Elements:
        #     if element.IsNot(KratosMultiphysics.ACTIVE):
        #         x_center = element.GetGeometry().Center().X
        #         y_center = element.GetGeometry().Center().X
        #         center_vector = KratosMultiphysics.Vector(3, 0.0)
        #         center_vector[0] = x_center - inital_point[0]
        #         center_vector[1] = y_center - inital_point[1]
        #         product = center_vector[0]*vector[0] + center_vector[1]*vector[1]
        #         if product > max_projection:
        #             max_projection = product
        #             max_y = element.GetGeometry().Center().Y
        #             max_x = element.GetGeometry().Center().X
        #             print(element.Id, max_projection, max_x, max_y)
        # wake_origin = KratosMultiphysics.Vector(3, 0.0)
        # wake_origin[0] = max_x - 0.00
        # wake_origin[1] = max_y


        max_x = -1e30
        for element in self.main_model_part.Elements:
            if element.IsNot(KratosMultiphysics.ACTIVE):
                x_center = element.GetGeometry().Center().X
                if x_center > max_x:
                    max_x = x_center
                    max_y = element.GetGeometry().Center().Y
                    max_elem = element
        wake_origin = KratosMultiphysics.Vector(3, 0.0)

        # for node in max_elem.GetNodes():
        #     node.SetValue(CPFApp.WING_TIP, True)
        # max_x_negative = -1e30
        # for node in self.main_model_part.Nodes:
        #     distance = node.GetSolutionStepValue(CPFApp.GEOMETRY_DISTANCE)
        #     if distance < 0.0 and node.X > max_x_negative:
        #         max_x_negative = node.X
        #         max_y_negative = node.Y + 1e-7
        # wake_origin = KratosMultiphysics.Vector(3, 0.0)
        self.max_inactive_x = max_x
        wake_origin[0] = max_x
        wake_origin[1] = max_y
        # self.wake_model_part.CreateNewNode(1, max_x, max_y, 0.0)
        # self.wake_model_part.CreateNewNode(2, self.model.GetModelPart("skin").GetNode(1).X, self.model.GetModelPart("skin").GetNode(1).Y, 0.0)
        # self.wake_model_part.CreateNewNode(3, 200.0, self.model.GetModelPart("skin").GetNode(1).Y, 0.0)
        # self.wake_model_part.CreateNewElement("Element2D2N", 1, [1,2], KratosMultiphysics.Properties(0))
        # self.wake_model_part.CreateNewElement("Element2D2N", 2, [2,3], KratosMultiphysics.Properties(0))

        # max_y +=  0.00001
        self.wake_model_part.CreateNewNode(1, max_x, max_y, 0.0)
        self.wake_model_part.CreateNewNode(2, 200.0, max_y, 0.0)
        # self.wake_model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        # self.wake_model_part.CreateNewNode(2, 200.0, 0.0, 0.0)
        self.wake_model_part.CreateNewElement("Element2D2N", 1, [1,2], KratosMultiphysics.Properties(0))
        angle = self.main_model_part.ProcessInfo.GetValue(CPFApp.ROTATION_ANGLE)
        self.main_model_part.ProcessInfo.SetValue(CPFApp.ROTATION_ANGLE, angle)
        print(angle)
        stop
        # print(self.model.GetModelPart("skin").GetNode(1).X, self.model.GetModelPart("skin").GetNode(1).Y)
        # self.moving_parameters["origin"].SetVector(wake_origin)
        # self.main_model_part.ProcessInfo.SetValue(CPFApp.WAKE_ORIGIN, wake_origin)
        # print("WAKE_ORIGIN", wake_origin)

        # self.moving_parameters.AddEmptyValue("rotation_angle")
        # self.moving_parameters["rotation_angle"].SetDouble(0.0)

        # print(self.moving_parameters)
        # CPFApp.MoveModelPartProcess(self.wake_model_part, self.moving_parameters).Execute()

    def ExecuteFinalizeSolutionStep(self):
        self.wake_sub_model_part = self.main_model_part.CreateSubModelPart("wake_sub_model_part")
        for elem in self.main_model_part.Elements:
            if elem.GetValue(CPFApp.WAKE) and elem.Is(KratosMultiphysics.ACTIVE) and elem.IsNot(KratosMultiphysics.MARKER):
                self.wake_sub_model_part.Elements.append(elem)

        absolute_tolerance = 1e-6
        CPFApp.PotentialFlowUtilities.CheckIfWakeConditionsAreFulfilled2D(self.wake_sub_model_part, absolute_tolerance, 0)
        self.main_model_part.RemoveSubModelPart("wake_sub_model_part")

    def CheckIfValid(self, trailing_edge_candidate):
        for elem in self.main_model_part.Elements:
            is_neighbour = False
            for node in elem.GetNodes():
                if (node.Id == trailing_edge_candidate.Id):
                    is_neighbour = True
            if is_neighbour:
                if elem.IsNot(KratosMultiphysics.ACTIVE):
                    return True
                for node in elem.GetNodes():
                    if node.GetSolutionStepValue(CPFApp.GEOMETRY_DISTANCE) < 0.0 and node.X < trailing_edge_candidate.X :
                        is_valid = self.CheckIfValid(node)
                        if is_valid:
                            return True
        return False

    def FindNode(self, max_x_coordinate, restarted_search_maximum):
        print(max_x_coordinate, restarted_search_maximum)
        for elem in self.main_model_part.Elements:
            if elem.Is(KratosMultiphysics.TO_SPLIT) and elem.Is(KratosMultiphysics.ACTIVE):
                for node in elem.GetNodes():
                    if(node.X > max_x_coordinate) and (node.X<restarted_search_maximum) and node.GetSolutionStepValue(CPFApp.GEOMETRY_DISTANCE) < 0.0:
                        max_x_coordinate = node.X
                        trailing_edge_node = node
        return trailing_edge_node
