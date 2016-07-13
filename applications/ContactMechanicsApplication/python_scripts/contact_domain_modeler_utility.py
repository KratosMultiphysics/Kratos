from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics 
import KratosMultiphysics.PfemBaseApplication as KratosPfemBase
import KratosMultiphysics.ContactMechanicsApplication as KratosContact
KratosMultiphysics.CheckForPreviousImport()

import meshing_domains_modeler_utility as modeler_utility

class ContactModelerUtility(modeler_utility.ModelerUtility):
    #

                      
    #
    def SetMeshingDomains(self, meshing_domains ):

        # set mesing domains
        self.meshing_domains = meshing_domains

        # set modeler utilities
        self.modeler_utils = KratosPfemBase.ModelerUtilities()

 
        # check if some of them is active:
        if( self.meshing_domains.size() ):
            self.remesh_domains_active = True                  
            
    ###
    #
    def RemeshDomains(self):
        
        if(self.remesh_domains_active):
            if( self.echo_level > 0 ):
                print("::[Modeler_Utility]:: CONTACT SEARCH : ", self.contact_condition)

            meshing_options = KratosMultiphysics.Flags()
            self.model_meshing =  KratosContact.ContactModelMeshing(self.model_part, meshing_options, self.echo_level)

            self.model_meshing.ExecuteInitialize()

            for domain in self.meshing_domains:

                domain.ExecuteMeshing();

                self.remesh_executed = True

 
            self.model_meshing.ExecuteFinalize()

            self.counter += 1 
            
            self.model_meshing.ExecuteFinalize()

