from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.ShallowWaterApplication as Shallow

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()


def CreateSolver(main_model_part, custom_settings):
    return HelmholtzSolver(main_model_part, custom_settings)

class HelmholtzSolver(object):

    def __init__(self, main_model_part, custom_settings):

        #TODO: shall obtain the compute_model_part from the MODEL once the object is implemented
        self.main_model_part = main_model_part

        ##settings string in json format
        #TODO: default settings
        #~ base_settings = KratosMultiphysics.Parameters("""
        #~ {
        #~ }""")

        ## Overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        #~ self.settings.ValidateAndAssignDefaults(base_settings)
        
        ## Construct the linear solvers
        import linear_solver_factory
        self.linear_solver = linear_solver_factory.ConstructSolver(self.settings["solver_settings"])

        ## Set the element replace settings
        self.settings.AddEmptyValue("element_replace_settings")
        if(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
            self.settings["element_replace_settings"] = KratosMultiphysics.Parameters("""
                {
                "element_name":"Element2D3N",
                "condition_name": "LineCondition2D2N"
                }
                """)
        else:
            raise Exception("domain size is not 2")

        print("Construction of HelmholtzSolver finished")

    def AddVariables(self):
        ## Add class variables
        self.main_model_part.AddNodalSolutionStepVariable(Shallow.ELEVATION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(Shallow.BATHYMETRY)
        self.main_model_part.AddNodalSolutionStepVariable(Shallow.DEPTH)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.GRAVITY)

        print("Shallow water solver variables added correctly")

    def ImportModelPart(self):
        ## Read model part
        self._ModelPartReading()
        ## Replace elements and conditions
        self._ExecuteAfterReading()
        ## Set buffer size
        self._SetBufferSize()

        print ("Model reading finished.")

    def AddDofs(self):
        ## Adding dofs
        for node in self.main_model_part.Nodes:
            node.AddDof(Shallow.ELEVATION)

        print("Shallow water solver DOFs added correctly.")

    def GetComputingModelPart(self):
        pass
        #~ return self.main_model_part.GetSubModelPart("fluid_computational_model_part")
