from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.DelaunayMeshingApplication as KratosDelaunay
KratosMultiphysics.CheckForPreviousImport()

#from multiprocessing import Pool

def Factory(settings, Model):
    if( not isinstance(settings,KratosMultiphysics.Parameters) ):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return RemeshDomainsProcess(Model, settings["Parameters"])


class RemeshDomainsProcess(KratosMultiphysics.Process):
    #
    def __init__(self, Model, custom_settings ):

        KratosMultiphysics.Process.__init__(self)

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "echo_level"            : 0,
            "model_part_name"       : "Meshing Domain",
            "meshing_control_type"  : "step",
            "meshing_frequency"     : 1.0,
            "meshing_before_output" : true,
            "meshing_domains"       : []
        }
        """)

        ##overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        self.echo_level        = self.settings["echo_level"].GetInt()
        self.meshing_frequency = self.settings["meshing_frequency"].GetDouble()

        self.meshing_control_is_time = False
        meshing_control_type   = self.settings["meshing_control_type"].GetString()
        if(meshing_control_type == "time"):
            self.meshing_control_is_time = True
        elif(meshing_control_type == "step"):
            self.meshing_control_is_time = False

        self.step_count   = 1
        self.counter      = 1
        self.next_meshing = 0.0
        self.meshing_before_output = self.settings["meshing_before_output"].GetBool()

        self.model = Model

        #construct meshing domains
        self.meshing_domains = []
        domains_list = self.settings["meshing_domains"]
        self.number_of_domains = domains_list.size()
        for i in range(0,self.number_of_domains):
            item = domains_list[i]
            domain_module = __import__(item["python_module"].GetString())
            domain = domain_module.CreateMeshingDomain(Model,item)
            self.meshing_domains.append(domain)

    #
    def ExecuteInitialize(self):

        self.main_model_part = self.model[self.settings["model_part_name"].GetString()]
        self.dimension = self.main_model_part.ProcessInfo[KratosMultiphysics.SPACE_DIMENSION]

        # mesh mesher initial values
        self.remesh_domains_active = False
        for domain in self.meshing_domains:
            if( domain.Active() ):
                self.remesh_domains_active = True

        # check restart
        self.restart = False
        if( self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] == True ):
            self.restart = True
            self.step_count = self.main_model_part.ProcessInfo[KratosMultiphysics.STEP]

            if self.meshing_control_is_time:
                self.next_meshing  = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME] + self.meshing_frequency
            else:
                self.next_meshing = self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] + self.meshing_frequency
        else:
            #if commented in the first step meshing is applied
            self.next_meshing = self.meshing_frequency
            # it must be initialized if restart is called//no

        self.main_model_part.ProcessInfo.SetValue(KratosDelaunay.INITIALIZED_DOMAINS, False);

        # initialize all meshing domains
        if( self.remesh_domains_active ):

            self.InitializeDomains()

            for domain in self.meshing_domains:
                domain.SetEchoLevel(self.echo_level)
                domain.Initialize()
                #domain.Check()

        print(self._class_prefix()+" Ready")


    def InitializeDomains(self):

        print(self._class_prefix()+" Initialize Domains")
        import domain_utilities
        domain_utils = domain_utilities.DomainUtilities()
        domain_utils.InitializeDomains(self.main_model_part,self.echo_level)

    ###

    #
    def ExecuteInitializeSolutionStep(self):

        self.step_count += 1

    #
    def ExecuteBeforeOutputStep(self):

        if(self.remesh_domains_active):
            if( self.meshing_before_output ):
                if(self.IsMeshingStep()):
                    self.RemeshDomains()

    #
    def ExecuteAfterOutputStep(self):

        if(self.remesh_domains_active):
            if( not self.meshing_before_output ):
                if(self.IsMeshingStep()):
                    self.RemeshDomains()

    ###

    #
    def ExecuteMeshing(domain):
        domain.ExecuteMeshing()

    #
    def RemeshDomains(self):

        if( self.echo_level >= 0 ):
            print(self._class_prefix()+" [ Mesh Generation (call:"+str(self.counter)+") ]")

        self.model_manager = self.GetModelManager()

        self.model_manager.ExecuteInitialize()

        #serial
        for domain in self.meshing_domains:
            domain.ExecuteMeshing()

        if(self.echo_level>1):
            print("")
            print(self.main_model_part)

        self.model_manager.ExecuteFinalize()

        self.counter += 1

        # set meshing step time
        self.SetMeshingStepTime()

        # schedule next meshing
        if(self.meshing_frequency > 0.0): # note: if == 0 always active
            if(self.meshing_control_is_time):
                time = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]
                while(self.next_meshing <= time):
                    self.next_meshing += self.meshing_frequency
            else:
                while(self.next_meshing <= self.step_count):
                    self.next_meshing += self.meshing_frequency


    #
    def GetModelManager(self):
        meshing_options = KratosMultiphysics.Flags()
        return KratosDelaunay.ModelStructure(self.main_model_part, meshing_options, self.echo_level)

    #
    def GetMeshingStep(self):
        return self.counter

    #
    def SetMeshingStepTime(self):
        current_time = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]
        self.main_model_part.ProcessInfo.SetValue(KratosDelaunay.MESHING_STEP_TIME, current_time)

    #
    def IsMeshingStep(self):

        if(self.meshing_control_is_time):
            #print( str(self.main_model_part.ProcessInfo[KratosMultiphysics.TIME])+">"+ str(self.next_meshing) )
            return ( self.main_model_part.ProcessInfo[KratosMultiphysics.TIME] >= self.next_meshing )
        else:
            return ( self.step_count >= self.next_meshing )

    #
    def GetVariables(self):
        import domain_utilities
        nodal_variables = domain_utilities.DomainUtilities().GetVariables()
        nodal_variables = nodal_variables + ['DETERMINANT_F'] # variables smoothing
        nodal_variables = nodal_variables + ['MEAN_ERROR'] # removing nodes

        #nodal_variables = nodal_variables + ['CAUCHY_STRESS_VECTOR', 'DEFORMATION_GRADIENT'] # transfer variables
        for domain in self.meshing_domains:
            nodal_variables = nodal_variables + domain.GetVariables()

        # print(self._class_prefix()+" Variables added")

        return nodal_variables

    #
    @classmethod
    def _class_prefix(self):
        header = "::[--Meshing_Process--]::"
        return header
