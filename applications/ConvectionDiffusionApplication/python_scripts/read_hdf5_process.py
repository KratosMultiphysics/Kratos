from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import bisect as bi
import numpy as np
import h5py
import KratosMultiphysics as Kratos
import KratosMultiphysics.ConvectionDiffusionApplication as CD
import KratosMultiphysics.ConvectionDiffusionApplication.CDProjectionModule as CDProjectionModule
import KratosMultiphysics.ConvectionDiffusionApplication.WriteMdpaToHdf5 as WriteMdpaToHdf5
from KratosMultiphysics import *
import sys
import os
import random

def Factory(settings, Model):
    if(type(settings) != Kratos.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ReadHdf5Process(Model, settings["Parameters"])

class ReadHdf5Process(Kratos.Process):
    def __init__(self, Model, settings):
        Kratos.Process.__init__(self)

        default_settings = Kratos.Parameters("""
            {
                "file_name" : "PLEASE_SPECIFY_HDF5_FILENAME",
                "model_part_name" : "please_specify_model_part_name",
                "reference_model_part_name" : "specify_model_part_name",
                "list_of_variables" : ["VELOCITY", "TEMPERATURE", "COORDINATES"],
                "error_projection_parameters" : {
                	"u_characteristic"  : 1.0
            	},
                "domain_size" : 3
            }
            """
            )

        settings.ValidateAndAssignDefaults(default_settings)
        self.parameters = settings
        self.problem_name = self.parameters["file_name"].GetString()
        self.reference_mdpa_name = str(os.getcwd() + '/' + self.parameters["reference_model_part_name"].GetString())
        self.destination_model_part = Model[settings["model_part_name"].GetString()]
        self.u_characteristic = settings["error_projection_parameters"]["u_characteristic"].GetDouble()
        self.domain_size = self.parameters["domain_size"].GetInt()
        for element in self.destination_model_part.Elements:
            rho = element.Properties.GetValue(Kratos.DENSITY)
            break

        self.p_characteristic = (1/2)*rho*self.u_characteristic**2
        self.nodal_velocity_error = []
        self.nodal_concentration_error = []
        self.step = []
    def ExecuteFinalizeSolutionStep(self):
        if self.destination_model_part.ProcessInfo[Kratos.STEP] % 10 == 0:
            self.reference_model_part = self.test_HDF5NodalSolutionStepDataIO(self.destination_model_part)
            self.fluid_variables = []
            self.projected_variables = []
            self.fluid_variables += [Kratos.VELOCITY]
            self.fluid_variables += [Kratos.TEMPERATURE]
            self.projected_variables += [CD.VELOCITY_PROJECTED]
            self.projected_variables += [CD.CONCENTRATION_PROJECTED]
            #Instead of defining variables in convection_diffusion_application:
            # self.element_name = "Element3D4N"
            # self.model = KratosMultiphysics.Model()
            # self.destination_model_part_new = self.model.CreateModelPart("NewDestinationModelPart")
            # self.destination_model_part_new.AddNodalSolutionStepVariable(CD.VELOCITY_PROJECTED)
            # self.destination_model_part_new.AddNodalSolutionStepVariable(CD.CONCENTRATION_PROJECTED)
            # model_part_cloner = KratosMultiphysics.ConnectivityPreserveModeler()
            # model_part_cloner.GenerateModelPart(self.destination_model_part, self.destination_model_part_new, self.element_name)
            CDProjectionModule.ProjectionModuleConvectionDiffusion(self.reference_model_part, self.destination_model_part, self.parameters, self.fluid_variables, self.projected_variables, self.domain_size)
            velocity_error, concentration_error = self.CalculateErrorNorm()
            WriteMdpaToHdf5.WriteConvergenceNodalErrorToHdf5(self.destination_model_part, velocity_error, concentration_error, self.problem_name)
    

    def ReadModelPartFile(self): 
        self.model_exact = Kratos.Model()
        model_part = self.model_exact.CreateModelPart('reference_model_part')
        model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        model_part.AddNodalSolutionStepVariable(Kratos.TEMPERATURE)

        # mdpa_name = '/home/aitor/Escritorio/Norouzi_mesh_1.gid/Norouzi_mesh_1'
        mdpa_name = self.reference_mdpa_name
        model_part_io_exactMesh = Kratos.ModelPartIO(mdpa_name)
        model_part_io_exactMesh.ReadModelPart(model_part)
        model_part.SetBufferSize(2)

        return model_part

    def GetFieldDataFile(self, model):
        import h5py
        # file_name = '/home/aitor/Escritorio/Norouzi_mesh_1.gid/Norouzi_mesh_1'
        file_name = self.reference_mdpa_name
        self.hf = h5py.File(file_name+".hdf5", 'r')
        group_names = list(self.hf.keys())
        self.velocities = self.hf[str(self.destination_model_part.ProcessInfo[Kratos.STEP])+'/VELOCITIES']
        self.temperatures = self.hf[str(self.destination_model_part.ProcessInfo[Kratos.STEP])+'/TEMPERATURES']
        self.coordinates = self.hf[str(self.destination_model_part.ProcessInfo[Kratos.STEP])+'/COORDINATES']

    def WriteDataInModelPart(self, exact_model_part, velocities, temperatures, coordinates):
        for node_i, node in enumerate(self.exact_model_part.Nodes):
            node.SetSolutionStepValue(Kratos.VELOCITY_X, self.velocities[node_i][0])
            node.SetSolutionStepValue(Kratos.VELOCITY_Y, self.velocities[node_i][1])
            node.SetSolutionStepValue(Kratos.VELOCITY_Z, self.velocities[node_i][2])
            node.SetSolutionStepValue(Kratos.TEMPERATURE, self.temperatures[node_i])
        self.hf.close()
    
    def test_HDF5NodalSolutionStepDataIO(self, model):
        self.exact_model_part = self.ReadModelPartFile()
        self.GetFieldDataFile(self.destination_model_part)
        self.WriteDataInModelPart(self.exact_model_part, self.velocities, self.temperatures, self.coordinates)
    
        return self.exact_model_part
                
    def CalculateErrorNorm(self):
        self.velocity_error_norm = self.VectorL2ErrorNorm(self.destination_model_part)
        self.concentration_error_norm = self.ScalarL2ErrorNorm(self.destination_model_part)
        return self.velocity_error_norm/self.u_characteristic, self.concentration_error_norm/self.p_characteristic

    def VectorL2ErrorNorm(self, destination_model_part):
        return CD.MeshConvergenceErrorCalculator().GetL2VectorErrorNorm(self.destination_model_part)

    def ScalarL2ErrorNorm(self, destination_model_part):
        return CD.MeshConvergenceErrorCalculator().GetL2ScalarErrorNorm(self.destination_model_part)

