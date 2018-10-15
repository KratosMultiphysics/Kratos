from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.DelaunayMeshingApplication as KratosDelaunay
import KratosMultiphysics.ContactMechanicsApplication as KratosContact
KratosMultiphysics.CheckForPreviousImport()


def Factory(settings, Model):
    if( not isinstance(settings,KratosMultiphysics.Parameters) ):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return RigidBodiesProcess(Model, settings["Parameters"])


class RigidBodiesProcess(KratosMultiphysics.Process):
    #
    def __init__(self, Model, custom_settings ):

        KratosMultiphysics.Process.__init__(self)

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "model_part_name" : "Solid Domain",
            "rigid_bodies"    : []
        }
        """)

        ##overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        self.echo_level        = 1

        self.Model = Model

    #
    def ExecuteInitialize(self):

        self.main_model_part = self.Model[self.settings["model_part_name"].GetString()]
        self.dimension       = self.main_model_part.ProcessInfo[KratosMultiphysics.SPACE_DIMENSION]

        # construct rigid body domains
        self.rigid_bodies = []
        bodies_list = self.settings["rigid_bodies"]
        self.number_of_bodies = bodies_list.size()
        for i in range(0,self.number_of_bodies):
            item = bodies_list[i]
            rigid_body_module = __import__(item["python_module"].GetString())
            body = rigid_body_module.CreateRigidBody( self.main_model_part, item )
            self.rigid_bodies.append(body)

        # initialize rigid body domains
        import domain_utilities
        domain_utils = domain_utilities.DomainUtilities()
        domain_utils.InitializeDomains(self.main_model_part,self.echo_level)

        for body in self.rigid_bodies:
            body.Initialize();

        print(self._class_prefix()+" Ready")

    ###

    #
    def ExecuteInitializeSolutionStep(self):
        pass


    #
    def ExecuteFinalizeSolutionStep(self):
        pass


    ###

    #
    @classmethod
    def GetVariables(self):
        nodal_variables = ['RIGID_WALL']
        return nodal_variables

    #
    @classmethod
    def _class_prefix(self):
        header = "::[----Rigid Bodies---]::"
        return header
