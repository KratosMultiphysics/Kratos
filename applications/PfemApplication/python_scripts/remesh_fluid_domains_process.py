from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.PfemApplication as KratosPfem
import KratosMultiphysics.DelaunayMeshingApplication as KratosDelaunay
KratosMultiphysics.CheckForPreviousImport()

import remesh_domains_process

def Factory(settings, Model):
    if( not isinstance(settings,KratosMultiphysics.Parameters) ):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return RemeshFluidDomainsProcess(Model, settings["Parameters"])


class RemeshFluidDomainsProcess(remesh_domains_process.RemeshDomainsProcess):
    #
    def __init__(self, Model, custom_settings ):

        super(RemeshFluidDomainsProcess,self).__init__(Model,custom_settings)

    #
    def ExecuteInitialize(self):

        remesh_domains_process.RemeshDomainsProcess.ExecuteInitialize(self)

        if(self.remesh_domains_active):
            self.counter = 0
            self.RemeshDomains()

    #
    def GetModelManager(self):
        meshing_options = KratosMultiphysics.Flags()
        self.mesher_utils = KratosDelaunay.MesherUtilities()
        meshing_options.Set(self.mesher_utils.KEEP_ISOLATED_NODES, True)

        return KratosDelaunay.ModelStructure(self.main_model_part, meshing_options, self.echo_level)

    #
    @classmethod
    def _class_prefix(self):
        header = "::[---Meshing_Fluid---]::"
        return header
