import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp
import KratosMultiphysics.MeshingApplication as MeshingApplication
from KratosMultiphysics.CompressiblePotentialFlowApplication.define_wake_process_2d import DefineWakeProcess2D
import time


def Factory(settings, Model):
    if( not isinstance(settings,KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    return DefineEmbeddedWakeProcess(Model, settings["Parameters"])

class DefineEmbeddedWakeProcess(DefineWakeProcess2D):

    def ExecuteInitialize(self):
        ini_time = time.time()

        self.wake_model_part=self.model.CreateModelPart("wake")

        self.wake_model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        self.wake_model_part.CreateNewNode(2, 200.0, 0.0, 0.0)
        self.wake_model_part.CreateNewElement("Element2D2N", 1, [1,2], KratosMultiphysics.Properties(0))

        self.moving_parameters = KratosMultiphysics.Parameters()
        self.moving_parameters.AddEmptyValue("origin")
        self.moving_parameters["origin"].SetVector(self.fluid_model_part.ProcessInfo.GetValue(CPFApp.WAKE_ORIGIN))
        self.moving_parameters.AddEmptyValue("rotation_angle")
        self.moving_parameters["rotation_angle"].SetDouble(self.fluid_model_part.ProcessInfo.GetValue(CPFApp.ROTATION_ANGLE))
        CPFApp.MoveModelPartProcess(self.wake_model_part, self.moving_parameters).Execute()

        for elem in self.fluid_model_part.Elements:
            boolean = elem.Is(KratosMultiphysics.TO_SPLIT)
            elem.Set(KratosMultiphysics.BOUNDARY,boolean)

        CPFApp.DefineEmbeddedWakeProcess(self.fluid_model_part, self.wake_model_part).Execute()

        for elem in self.fluid_model_part.Elements:
            boolean = elem.Is(KratosMultiphysics.BOUNDARY)
            elem.Set(KratosMultiphysics.TO_SPLIT,boolean)

        KratosMultiphysics.Logger.PrintInfo('EmbeddedWake','Wake computation time: ',time.time()-ini_time)

    def _SaveTrailingEdgeNode(self):
        '''
        This function finds and saves the trailing edge for further computations and deactivates the neighbour elements
        '''
        import time
        ini_time=time.time()
        self.__FindMaximumInactiveNode()
        KratosMultiphysics.Logger.PrintInfo('EmbeddedWake','Find Max. Inactive Node time: ',time.time()-ini_time)
        ini_time=time.time()
        self.__CreateInactiveTrailingEdgeElementBall()
        KratosMultiphysics.Logger.PrintInfo('EmbeddedWake','Creating Inactive Element Ball time: ',time.time()-ini_time)

        import math



        # x0 =  -0.75#self.origin[0]
        # y0 =  0.0#self.origin[1]
        # direction = KratosMultiphysics.Vector(3)
        # angle=math.radians(5.0)
        # direction[0]=-math.cos(angle)
        # direction[1]=math.sin(angle)
        # n = KratosMultiphysics.Vector(3)
        # xn = KratosMultiphysics.Vector(3)
        # n[0] = -direction[1]
        # n[1] = direction[0]
        # n[2] = 0.0
        # max_x_center = -1e10
        # for elem in self.boundary_sub_model_part.Elements:
        #     distances = KratosMultiphysics.Vector(len(elem.GetNodes()))
        #     '''with positive epsilon'''
        #     npos = 0
        #     nneg = 0
        #     for elnode in elem.GetNodes():
        #         xn[0] = elnode.X - x0
        #         xn[1] = elnode.Y - y0
        #         xn[2] = 0.0
        #         d =  xn[0]*n[0] + xn[1]*n[1]
        #         if(d < 0):
        #             nneg += 1
        #         else:
        #             npos += 1
        #         elnode.SetValue(KratosMultiphysics.TEMPERATURE,d)
        #     if(nneg>0 and npos>0):
        #         center_X=elem.GetGeometry().Center().X
        #         if (center_X>x0):
        #             elem.Set(KratosMultiphysics.THERMAL,True)
        #             elem.Set(KratosMultiphysics.ACTIVE,False)
        #             elem.Set(KratosMultiphysics.TO_SPLIT,False)
        #             if (center_X > max_x_center):
        #                 max_x_center = center_X
        #                 max_elem = elem

        # max_elem.Set(KratosMultiphysics.ACTIVE,False)
        # max_elem.Set(KratosMultiphysics.TO_SPLIT,False)
        # max_x_node=-1e10
        # for elnode in max_elem.GetNodes():
        #     if elnode.X>max_x_node:
        #         max_node=elnode
        #         max_x_node=elnode.X
        # self.trailing_edge_node=max_node
##################################################################
####################################################################

        # self.wake_model_part=self.model.GetModelPart("wake")

        # ini_time = time.time()
        # for elem in self.fluid_model_part.Elements:
        #     boolean = elem.Is(KratosMultiphysics.TO_SPLIT)
        #     elem.Set(KratosMultiphysics.BOUNDARY,boolean)

        # distance_calculator = KratosMultiphysics.CalculateDiscontinuousDistanceToSkinProcess2D(self.fluid_model_part, self.wake_model_part)
        # distance_calculator.Execute()

        # for elem in self.fluid_model_part.Elements:
        #     # print(elem.GetValue(KratosMultiphysics.ELEMENTAL_DISTANCES))
        #     boolean = elem.Is(KratosMultiphysics.BOUNDARY)
        #     elem.Set(KratosMultiphysics.TO_SPLIT,boolean)
        # KratosMultiphysics.Logger.PrintInfo('EmbeddedWake','Computing Wake Distance: ',time.time()-ini_time)

        # # stop
        # x0 =  0.0#self.origin[0]
        # y0 =  0.0#self.origin[1]

        # max_x_center = -1e10
        # for elem in self.boundary_sub_model_part.Elements:
        #     distances = elem.GetValue(KratosMultiphysics.ELEMENTAL_DISTANCES)
        #     '''with positive epsilon'''
        #     npos = 0
        #     nneg = 0
        #     counter = -1
        #     for elnode in elem.GetNodes():
        #         counter += 1
        #         d =  distances[counter]
        #         if(d < 0):
        #             nneg += 1
        #         else:
        #             npos += 1
        #         elnode.SetValue(KratosMultiphysics.TEMPERATURE,d)
        #     if(nneg>0 and npos>0):
        #         center_X=elem.GetGeometry().Center().X
        #         if (center_X>x0):
        #             elem.Set(KratosMultiphysics.THERMAL,True)
        #             elem.Set(KratosMultiphysics.ACTIVE,False)
        #             # elem.Set(KratosMultiphysics.TO_SPLIT,False)
        #             if (center_X > max_x_center):
        #                 max_x_center = center_X
        #                 max_elem = elem

        # max_elem.Set(KratosMultiphysics.ACTIVE,False)
        # # max_elem.Set(KratosMultiphysics.TO_SPLIT,False)
        # max_x_node=-1e10
        # for elnode in max_elem.GetNodes():
        #     if elnode.X>max_x_node:
        #         max_node=elnode
        #         max_x_node=elnode.X
        # self.trailing_edge_node=max_node


        # self.__DeactivateActive(self.fluid_model_part.GetElement(86292))

        # self.trailing_edge_node = self.fluid_model_part.GetNode(79576)

        self.trailing_edge_node.SetValue(CPFApp.TRAILING_EDGE, True)

    def _MarkKuttaElements(self):
        # This function selects the kutta elements. Kutta elements
        # are touching the trailing edge from below.
        for elem in self.trailing_edge_model_part.Elements:
            elem.SetValue(CPFApp.KUTTA, True)
            for node in elem.GetNodes():
                if not node.GetValue(KratosMultiphysics.CompressiblePotentialFlowApplication.TRAILING_EDGE):
                    # # Compute the distance from the element center to the trailing edge
                    # x_distance_to_te = elem.GetGeometry().Center().X - self.trailing_edge_node.X
                    # y_distance_to_te = elem.GetGeometry().Center().Y - self.trailing_edge_node.Y
                    # Compute the distance from the element center to the trailing edge
                    x_distance_to_te = node.X - self.trailing_edge_node.X
                    y_distance_to_te = node.Y - self.trailing_edge_node.Y

                    # Compute the projection of the distance vector in the wake normal direction
                    distance_to_wake = x_distance_to_te*self.wake_normal[0] + \
                        y_distance_to_te*self.wake_normal[1]

                    print(elem.Id, distance_to_wake)
                    if(distance_to_wake > 0.0):
                        print(elem.Id, distance_to_wake)

                        elem.SetValue(CPFApp.KUTTA, False)

        # self.fluid_model_part.GetElement(75420).SetValue(CPFApp.KUTTA,True)

    def __FindMaximumInactiveNode(self):
        '''
        This function finds the maximum-x node that is contained in the negative part of the embedded body
        '''
        self.boundary_sub_model_part=self.fluid_model_part.CreateSubModelPart("boundary")
        max_x=-1e10
        boundary_element_id_list=[]
        for element in self.fluid_model_part.Elements:
            if element.Is(KratosMultiphysics.TO_SPLIT):
                boundary_element_id_list.append(element.Id)
            if element.IsNot(KratosMultiphysics.ACTIVE):
                for node in element.GetNodes():
                    if node.X>max_x:
                        max_x=node.X
                        self.max_inactive_node=node
        self.boundary_sub_model_part.AddElements(boundary_element_id_list)
        self.trailing_edge_node = self.max_inactive_node
    def __CreateInactiveTrailingEdgeElementBall(self):
        '''
        This function deactivates all the elements near the trailing edge node to consistenly attach the wake.
        '''
        self.deactivated_model_part=self.fluid_model_part.CreateSubModelPart("deactivated")
        deactivated_element_id_list=[]
        for element in self.boundary_sub_model_part.Elements:
            for elnode in element.GetNodes():
                if elnode.Id == self.max_inactive_node.Id:
                    self.__DeactivateActive(element)
                    deactivated_element_id_list.append(element.Id)
        self.deactivated_model_part.AddElements(deactivated_element_id_list)
        self.__FindTrailingEdgeNodeFromInactiveBall()

    def __FindTrailingEdgeNodeFromInactiveBall(self):
        '''
        This function find the furthest inactive node that will be defined as the trailing edge
        '''
        max_x_node=-1e10
        max_x_center=-1e10
        for element in self.deactivated_model_part.Elements:
            n_center = 0
            center_X=element.GetGeometry().Center().X
            if center_X>max_x_center:
                max_x_center=center_X
                max_elem=element

        for elnode in max_elem.GetNodes():
            if elnode.X>max_x_node:
                max_x_node=elnode.X
                max_node=elnode

        for element in self.boundary_sub_model_part.Elements:
            center_X=element.GetGeometry().Center().X
            if center_X>self.max_inactive_node.X:
                self.__DeactivateBoundary(element)

        self.trailing_edge_node = max_node

    def __DeactivateBoundary(self,elem):
        elem.Set(KratosMultiphysics.TO_SPLIT,False)

    def __DeactivateActive(self,elem):
        elem.Set(KratosMultiphysics.ACTIVE,False)
        elem.Set(KratosMultiphysics.TO_SPLIT,False)

