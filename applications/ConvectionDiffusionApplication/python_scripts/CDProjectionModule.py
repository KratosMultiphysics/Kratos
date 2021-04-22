import KratosMultiphysics as Kratos
from KratosMultiphysics import Parameters
import KratosMultiphysics.ConvectionDiffusionApplication as CD
import sys
import numpy as np

class ProjectionModuleConvectionDiffusion:

    def __init__(self,
                reference_model_part,
                destination_model_part,
                project_parameters,
                fluid_variables,
                projected_variables,
                domain_size):

        self.reference_model_part = reference_model_part
        self.destination_model_part = destination_model_part
        self.project_parameters = project_parameters
        self.fluid_variables = fluid_variables
        self.projected_variables = projected_variables
        self.domain_size = domain_size
        if self.domain_size == 2:
            self.projector = CD.BinBasedInterpolationUtility2D()
            self.bin_of_objects_fluid = Kratos.BinBasedFastPointLocator2D(self.reference_model_part)
        elif self.domain_size == 3:
            self.projector = CD.BinBasedInterpolationUtility3D()
            self.bin_of_objects_fluid = Kratos.BinBasedFastPointLocator3D(self.reference_model_part)
        HMin = 0.00001
        self.bin_of_objects_fluid.UpdateSearchDatabaseAssignedSize(HMin)
        for var in self.fluid_variables:
            self.projector.AddFluidVariable(var)
        for var in self.projected_variables:
            self.projector.AddVariablesToImposeProjection(var)
        self.ProjectFromFluid()
        
    def ProjectFromFluid(self):
        self.projector.InterpolateFromFluidMeshCD(self.reference_model_part,
                                                self.destination_model_part,
                                                self.project_parameters,
                                                self.bin_of_objects_fluid)

    
