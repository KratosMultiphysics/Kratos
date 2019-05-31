import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp
import KratosMultiphysics.MeshingApplication as MeshingApplication
from KratosMultiphysics.CompressiblePotentialFlowApplication.define_wake_process_2d import DefineWakeProcess2D
import math

def Factory(settings, Model):
    if( not isinstance(settings,KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    return DefineEmbeddedWakeProcess(Model, settings["Parameters"])

class DefineEmbeddedWakeProcess(DefineWakeProcess2D):

    def _SaveTrailingEdgeNode(self):
        # This function finds and saves the trailing edge for further computations
        max_x=-1e10
        print("SAVING TRAILING EDGE NODE")
        boundary_sub_model_part=self.fluid_model_part.CreateSubModelPart("boundary")
        for element in self.fluid_model_part.Elements:
            if element.Is(KratosMultiphysics.TO_SPLIT):
                boundary_sub_model_part.Elements.append(element)
            if element.IsNot(KratosMultiphysics.ACTIVE):
                for node in element.GetNodes():
                    if node.X>max_x:
                        max_x=node.X
                        max_inactive_node=node
        for element in boundary_sub_model_part.Elements:
            for node in element.GetNodes():
                if node.X > max_x:
                    element.Set(KratosMultiphysics.TO_SPLIT,False)




        max_x_node=-1e10
        deactivated_model_part=self.fluid_model_part.CreateSubModelPart('deactivated')
        for element in boundary_sub_model_part.Elements:
            for elnode in element.GetNodes():
                if elnode.Id == max_inactive_node.Id:
                    self.DeactivateActive(element)
                    deactivated_model_part.Elements.append(element)
        max_x_center=-1e10
        for element in deactivated_model_part.Elements:
            n_center = 0
            center_X=element.GetGeometry().Center().X
            if center_X>max_x_center:
                max_x_center=center_X
                max_elem=element

        for elnode in max_elem.GetNodes():
            if elnode.X>max_x_node:
                max_x_node=elnode.X
                max_node=elnode

        self.trailing_edge_node = max_node

        for element in boundary_sub_model_part.Elements:
            center_X=element.GetGeometry().Center().X
            if center_X>max_inactive_node.X:
                self.DeactivateBoundary(element)

        self.trailing_edge_node.SetValue(CPFApp.TRAILING_EDGE, True)

    def DeactivateElement(self,elem):
        self.fluid_model_part.GetElement(elem).Set(KratosMultiphysics.ACTIVE,False)

    def DeactivateBoundary(self,elem):
        elem.Set(KratosMultiphysics.TO_SPLIT,False)
    def DeactivateActive(self,elem):
        elem.Set(KratosMultiphysics.ACTIVE,False)
        elem.Set(KratosMultiphysics.TO_SPLIT,False)
