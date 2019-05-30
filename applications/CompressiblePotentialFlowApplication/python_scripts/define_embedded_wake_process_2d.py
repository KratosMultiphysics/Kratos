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
        for element in self.fluid_model_part.Elements:
            if element.IsNot(KratosMultiphysics.ACTIVE):
                for node in element.GetNodes():
                    if node.X>max_x:
                        max_x=node.X
                        self.trailing_edge_node = node
        self.trailing_edge_node.SetValue(CPFApp.TRAILING_EDGE, True)

