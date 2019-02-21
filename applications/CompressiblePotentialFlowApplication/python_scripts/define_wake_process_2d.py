import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp
import math


def Factory(settings, Model):
    if(not isinstance(settings, KratosMultiphysics.Parameters)):
        raise Exception(
            "expected input shall be a Parameters object, encapsulating a json string")

    return DefineWakeProcess(Model, settings["Parameters"])


class DefineWakeProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
            {
                "mesh_id"                   : 0,
                "model_part_name"           : "please specify the model part that contains the kutta nodes",
                "upper_surface_model_part_name" : "please specify the model part that contains the upper surface nodes",
                "lower_surface_model_part_name" : "please specify the model part that contains the lower surface nodes",
                "fluid_part_name"           : "MainModelPart",
                "wake_direction"                 : [1.0,0.0,0.0],
                "epsilon"    : 1e-9
            }
            """)

        settings.ValidateAndAssignDefaults(default_settings)
        # TODO Implement this process in C++ and make it open mp parallel to save time selecting the wake elements

        self.wake_direction = settings["wake_direction"].GetVector()
        if(self.wake_direction.Size() != 3):
            raise Exception('The wake direction should be a vector with 3 components!')

        dnorm = math.sqrt(
            self.wake_direction[0]**2 + self.wake_direction[1]**2 + self.wake_direction[2]**2)
        self.wake_direction[0] /= dnorm
        self.wake_direction[1] /= dnorm
        self.wake_direction[2] /= dnorm

        self.wake_normal = KratosMultiphysics.Vector(3)
        self.wake_normal[0] = -self.wake_direction[1]
        self.wake_normal[1] = self.wake_direction[0]
        self.wake_normal[2] = 0.0

        self.epsilon = settings["epsilon"].GetDouble()

        self.upper_surface_model_part = Model[settings["upper_surface_model_part_name"].GetString(
        )]
        self.lower_surface_model_part = Model[settings["lower_surface_model_part_name"].GetString(
        )]

        self.fluid_model_part = Model[settings["fluid_part_name"].GetString()].GetRootModelPart()
        self.trailing_edge_model_part = self.fluid_model_part.CreateSubModelPart(
            "trailing_edge_model_part")

        KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplex(self.fluid_model_part,
                                                                       self.fluid_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE])

        # Neigbour search tool instance
        AvgElemNum = 10
        AvgNodeNum = 10
        nodal_neighbour_search = KratosMultiphysics.FindNodalNeighboursProcess(
            self.fluid_model_part, AvgElemNum, AvgNodeNum)
        # Find neighbours
        nodal_neighbour_search.Execute()

    def ExecuteInitialize(self):
        # Save the trailing edge for further computations
        self.SaveTrailingEdgeNode()
        # Check which elements are cut and mark them as wake
        self.MarkWakeElements()
        # Mark the elements touching the trailing edge from below as kutta
        self.MarkKuttaElements()
        # Mark the trailing edge element that is further downstream as wake
        self.MarkWakeTEElement()

    def SaveTrailingEdgeNode(self):
        # This function finds and saves the trailing edge for further computations
        max_x_coordinate = -1e30
        for node in self.upper_surface_model_part.Nodes:
            if(node.X > max_x_coordinate):
                max_x_coordinate = node.X
                self.te = node

        self.te.SetValue(CPFApp.TRAILING_EDGE, True)

    def MarkWakeElements(self):
        # This function checks which elements are cut by the wake
        # and marks them as wake elements
        KratosMultiphysics.Logger.PrintInfo('...Selecting wake elements...')

        for elem in self.fluid_model_part.Elements:
            # Mark and save the elements touching the trailing edge
            self.MarkTrailingEdgeElements(elem)

            # Elements downstream the trailing edge can be wake elements
            potentially_wake = self.SelectPotentiallyWakeElements(elem)

            if(potentially_wake):
                # Compute the nodal distances of the element to the wake
                distances_to_wake = self.ComputeDistancesToWake(elem)

                # Selecting the cut (wake) elements
                wake_element = self.SelectWakeElements(distances_to_wake)

                if(wake_element):
                    elem.SetValue(CPFApp.WAKE, True)
                    elem.SetValue(
                        KratosMultiphysics.ELEMENTAL_DISTANCES, distances_to_wake)

        KratosMultiphysics.Logger.PrintInfo('...Selecting wake elements finished...')

    def MarkTrailingEdgeElements(self, elem):
        # This function marks the elements touching the trailing
        # edge and saves them in the trailing_edge_model_part for
        # further computations
        for elnode in elem.GetNodes():
            if(elnode.GetValue(CPFApp.TRAILING_EDGE)):
                elem.SetValue(CPFApp.TRAILING_EDGE, True)
                self.trailing_edge_model_part.Elements.append(elem)
                break

    def SelectPotentiallyWakeElements(self, elem):
        # This function selects the elements downstream the
        # trailing edge as potentially wake elements

        # Compute the distance from the element's center to
        # the trailing edge
        x_distance_to_te = elem.GetGeometry().Center().X - self.te.X
        y_distance_to_te = elem.GetGeometry().Center().Y - self.te.Y

        # Compute the projection of the distance in the wake direction
        projection_on_wake = x_distance_to_te*self.wake_direction[0] + \
            y_distance_to_te*self.wake_direction[1]

        # Elements downstream the trailing edge can be wake elements
        if(projection_on_wake > 0):
            return True
        else:
            return False

    def ComputeDistancesToWake(self, elem):
        # This function computes the distance of the element nodes
        # to the wake
        nodal_distances_to_wake = KratosMultiphysics.Vector(3)
        counter = 0
        for elnode in elem.GetNodes():
            # Compute the distance from the node to the trailing edge
            x_distance_to_te = elnode.X - self.te.X
            y_distance_to_te = elnode.Y - self.te.Y

            # Compute the projection of the distance vector in the wake normal direction
            distance_to_wake = x_distance_to_te*self.wake_normal[0] + \
                y_distance_to_te*self.wake_normal[1]

            # Nodes laying on the wake have a positive distance
            if(abs(distance_to_wake) < self.epsilon):
                distance_to_wake = self.epsilon

            nodal_distances_to_wake[counter] = distance_to_wake
            counter += 1

        return nodal_distances_to_wake

    @staticmethod
    def SelectWakeElements(distances_to_wake):
        # This function checks whether the element is cut by the wake

        # Initialize counters
        number_of_nodes_with_positive_distance = 0
        number_of_nodes_with_negative_distance = 0

        # Count how many element nodes are above and below the wake
        for nodal_distance_to_wake in distances_to_wake:
            if(nodal_distance_to_wake < 0.0):
                number_of_nodes_with_negative_distance += 1
            else:
                number_of_nodes_with_positive_distance += 1

        # Elements with nodes above and below the wake are wake elements
        return(number_of_nodes_with_negative_distance > 0 and number_of_nodes_with_positive_distance > 0)

    def MarkKuttaElements(self):
        # This function selects the kutta elements. Kutta elements
        # are touching the trailing edge from below.
        for elem in self.trailing_edge_model_part.Elements:
            # Compute the distance from the element center to the trailing edge
            x_distance_to_te = elem.GetGeometry().Center().X - self.te.X
            y_distance_to_te = elem.GetGeometry().Center().Y - self.te.Y

            # Compute the projection of the distance vector in the wake normal direction
            distance_to_wake = x_distance_to_te*self.wake_normal[0] + \
                y_distance_to_te*self.wake_normal[1]

            # Marking the elements under the trailing edge as kutta
            if(distance_to_wake < 0.0):
                elem.SetValue(CPFApp.KUTTA, True)

    def ComputeMaximumX(self):
        # This function computes the maximum x coordinate within
        # all the trailing edge elements' centers
        maximum_x_coordinate = -1e30

        # Loop over the elements touching the trailing edge
        for elem in self.trailing_edge_model_part.Elements:
            # Find the element touching the trailing edge with maximum x coordinate
            if(elem.GetGeometry().Center().X > maximum_x_coordinate):
                maximum_x_coordinate = elem.GetGeometry().Center().X

        return maximum_x_coordinate

    def MarkWakeTEElement(self):
        # This function finds the trailing edge element that is further downstream
        # and marks it as wake trailing edge element

        # Compute the maximum x coordinate within all the trailing edge
        # elements' centers
        maximum_x_coordinate = self.ComputeMaximumX()

        for elem in self.trailing_edge_model_part.Elements:
            if(elem.GetGeometry().Center().X < maximum_x_coordinate):
                elem.SetValue(CPFApp.WAKE, False)
            else:
                elem.Set(KratosMultiphysics.STRUCTURE)
                elem.SetValue(CPFApp.KUTTA, False)
