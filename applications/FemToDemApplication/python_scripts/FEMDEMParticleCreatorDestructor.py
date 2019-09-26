from KratosMultiphysics.DEMApplication import *

def Wait():
    input("Press Something")


class FemDemParticleCreatorDestructor:

    def __init__(self, SpheresModelPart, Properties, Parameters):
        self.SpheresModelPart = SpheresModelPart
        self.Properties = Properties
        self.creator_destructor = ParticleCreatorDestructor()
        self.Parameters = Parameters

    def FEMDEM_CreateSphericParticle(self, Coordinates, Radius, Id):

        element_name = self.Parameters["ElementType"].GetString()

        # only work with this name WHYYYY
        element_name = "SphericParticle3D"

        created_element = self.creator_destructor.CreateSphericParticle(self.SpheresModelPart,
                                                                        Coordinates,
                                                                        self.Properties,
                                                                        Radius,
                                                                        element_name)
        created_node = created_element.GetNodes()[0]
        created_node.Id = Id

    def GetSpheresModelPart(self):
        return self.SpheresModelPart

    def ClearElementsAndNodes(self):
        self.SpheresModelPart.Elements.clear()
        self.SpheresModelPart.Nodes.clear()
