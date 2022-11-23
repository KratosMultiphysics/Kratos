# Importing the Kratos Library
import KratosMultiphysics

# Import the functions I need
from Functions import *

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return SBMBoundaryConditionsProcess(Model, settings["Parameters"])


## All the processes python should be derived from "Process"
class SBMBoundaryConditionsProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings):
        KratosMultiphysics.Process.__init__(self)
        # default_settings = self.GetDefaultParameters()
        # settings.ValidateAndAssignDefaults(default_settings)
        
        # Ottengo la stringhe con i nomi dei model part 
        fluid_part_name = settings["model_part_name"].GetString()
        skin_part_name = settings["skin_model_part_name"].GetString()
        self.Model_Fluid = Model.GetModelPart(fluid_part_name)
        self.Model_Skin = Model.GetModelPart(skin_part_name)
        

    # @classmethod
    # def GetDefaultParameters(cls):
    #     default_parameters = KratosMultiphysics.Parameters("""
    #     {
    #     "model_part_name" : "fluid_model_part",
    #     "skin_model_part_name" : "skin_model_part"
    #     }
    #     """
    #     )
    #     return default_parameters

 

    def Execute(self):
        
        print('ciao, sto applicando sbm process')
        # Compute the sign distance from the skin_model_part
        KratosMultiphysics.CalculateDistanceToSkinProcess2D(self.Model_Fluid, self.Model_Skin).Execute()

        # Find the surrogate boundary nodes
        tot_sur_nodes = Find_surrogate_nodes(self.Model_Fluid)

        # Total number of skin elements
        tot_skin_el = len(self.Model_Skin.Conditions)
        print('Number of skin elements: ', tot_skin_el)


        # Find the closest skin element for each surr node
        closest_element = Find_closest_skin_element(self.Model_Fluid,self.Model_Skin,tot_sur_nodes,tot_skin_el)


        # Find the projection onto the skin elements for each surr node
        projection_surr_nodes = Find_projections(self.Model_Fluid,self.Model_Skin,tot_sur_nodes,tot_skin_el,closest_element)

        # Then we create a sub_model part with just the elements & the nodes "outside" the surrogate boundary
        sub_model_part_fluid = self.Model_Fluid.CreateSubModelPart("sub_model_part")
        for node in self.Model_Fluid.Nodes :
            if node.GetSolutionStepValue(KratosMultiphysics.DISTANCE) > 0 :
                sub_model_part_fluid.AddNode(node,0)
        for elem in self.Model_Fluid.Elements :
            if elem.Is(MARKER):
                sub_model_part_fluid.AddElement(elem,0)

        # Set a scalar field -> VELOCITY_X
        for node in sub_model_part_fluid.Nodes:
                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X, math.sqrt((node.X)**2+(node.Y)**2) )
                    node.SetValue(KratosMultiphysics.NODAL_AREA, 0.0)


        # Compute the gradient with the function ComputeNodalGradientProcess
        KratosMultiphysics.ComputeNodalGradientProcess(
        sub_model_part_fluid,
        KratosMultiphysics.VELOCITY_X,
        KratosMultiphysics.VELOCITY_X_GRADIENT,
        KratosMultiphysics.NODAL_AREA).Execute()

        # Compute the Dirichlet BC at the surr nodes
        surr_BC = Dirichlet_BC (self.Model_Fluid,self.Model_Skin,tot_sur_nodes,closest_element,projection_surr_nodes)
        
        return self.Model_Fluid, surr_BC, sub_model_part_fluid


