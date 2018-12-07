from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.DelaunayMeshingApplication as KratosDelaunay
import KratosMultiphysics.ContactMechanicsApplication as KratosContact
KratosMultiphysics.CheckForPreviousImport()

import remesh_domains_process

def Factory(settings, Model):
    if( not isinstance(settings,KratosMultiphysics.Parameters) ):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return ContactDomainProcess(Model, settings["Parameters"])


class ContactDomainProcess(remesh_domains_process.RemeshDomainsProcess):
    #
    def __init__(self, Model, custom_settings):

        super(ContactDomainProcess, self).__init__(Model, custom_settings)
    #
    def ExecuteInitialize(self):


        self.main_model_part = self.model[self.settings["model_part_name"].GetString()]
        self.dimension = self.main_model_part.ProcessInfo[KratosMultiphysics.SPACE_DIMENSION]


        # check restart
        self.restart = False
        if( self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] == True ):
            self.restart = True
            self.step_count = self.main_model_part.ProcessInfo[KratosMultiphysics.STEP]

            if self.meshing_control_is_time:
                self.next_meshing  = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]
            else:
                self.next_meshing = self.main_model_part.ProcessInfo[KratosMultiphysics.STEP]


        # execute initialize base class
        if( self.main_model_part.ProcessInfo[KratosDelaunay.INITIALIZED_DOMAINS] == False ):
            # self.main_model_part.ProcessInfo[KratosDelaunay.INITIALIZED_DOMAINS] == False
            self.InitializeDomains()
            print(" initialize domains ")

        # initialize contact domains
        for domain in self.meshing_domains:
            domain.SetEchoLevel(self.echo_level)
            domain.Initialize()

        if self.restart:
            self.RemeshDomains()

        print(self._class_prefix()+" Ready")

    #

    def ExecuteBeforeOutputStep(self):

        if(self._domain_parts_updated() or self.IsMeshingStep()):
            if self.meshing_before_output:
                self.RemeshDomains()
    #
    def ExecuteAfterOutputStep(self):

        if(self._domain_parts_updated() or self.IsMeshingStep()):
            if not self.meshing_before_output:
                self.RemeshDomains()

    ###
    def _domain_parts_updated(self):

        process_info = self.main_model_part.ProcessInfo
        if process_info.Has(KratosDelaunay.MESHING_STEP_TIME):
            current_time = process_info[KratosMultiphysics.TIME]
            delta_time = process_info[KratosMultiphysics.DELTA_TIME]
            previous_time = current_time - delta_time

            #arithmetic floating point tolerance
            tolerance = delta_time * 0.001

            meshing_step_time = process_info[KratosDelaunay.MESHING_STEP_TIME]

            if meshing_step_time > previous_time-tolerance and meshing_step_time < previous_time+tolerance:
                return True

        return False

    #
    def GetModelManager(self):
        meshing_options = KratosMultiphysics.Flags()
        return KratosContact.ContactModelStructure(self.main_model_part, meshing_options, self.echo_level)
    #
    def SetMeshingStepTime(self):
        current_time = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]
        self.main_model_part.ProcessInfo.SetValue(KratosContact.CONTACT_STEP_TIME, current_time)
    #
    def GetVariables(self):

        nodal_variables = super(ContactDomainProcess, self).GetVariables()

        nodal_variables = nodal_variables + ['OFFSET']
        nodal_variables = nodal_variables + ['CONTACT_NORMAL', 'CONTACT_FORCE']
        nodal_variables = nodal_variables + ['CONTACT_STRESS', 'EFFECTIVE_CONTACT_STRESS', 'EFFECTIVE_CONTACT_FORCE']

        return nodal_variables

    #
    @classmethod
    def _class_prefix(self):
        header = "::[--Meshing Contact--]::"
        return header
