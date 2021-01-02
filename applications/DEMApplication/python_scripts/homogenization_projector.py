import KratosMultiphysics as Kratos
from KratosMultiphysics import Parameters
import KratosMultiphysics.DEMApplication as DEM
import sys
import numpy as np

class ProjectionModule:

    def __init__(self,
                homogenization_model_part,
                balls_model_part,
                project_parameters,
                coupling_homogenization_vars,
                homogenization_type,
                time_filtered_vars,
                domain_size=3):

        self.homogenization_model_part = homogenization_model_part
        self.particles_model_part = balls_model_part
        self.project_parameters = project_parameters
        self.dimension = domain_size
        self.homogenization_type = project_parameters["homogenization"]["homogenization_type"].GetInt()
        #self.shape_factor = self.project_parameters["homogenization"]["shape_factor"].GetDouble()
        #self.meso_scale_length = self.project_parameters["homogenization"]["meso_scale_length"].GetDouble()
        self.averaging_time_interval = self.project_parameters["homogenization"]["averaging_time_interval"].GetInt()


        # Create projector_parameters
        self.projector_parameters = Parameters("{}")
        self.projector_parameters.AddValue("homogenization_type", project_parameters["homogenization"]["homogenization_type"])
        if self.dimension == 3:
            self.projector = DEM.BinBasedDEMHomogenizationMapper3D(self.projector_parameters)
            self.bin_of_objects_homogenization = Kratos.BinBasedFastPointLocator3D(homogenization_model_part)

        else:
            self.projector = DEM.BinBasedDEMHomogenizationMapper2D(self.projector_parameters)
            self.bin_of_objects_homogenization = Kratos.BinBasedFastPointLocator2D(homogenization_model_part)

        # telling the projector which variables we are interested in modifying

        for var in coupling_homogenization_vars:
        self.projector.AddHomogenizationCouplingVariable(var)

        for var in time_filtered_vars:
            averaging_time_interval = self.averaging_time_interval
            if self.averaging_time_interval < 15:
                alpha = 1 - np.exp(- averaging_time_interval)
            else:
                alpha = 1.0 / averaging_time_interval
            print('A'*100, alpha)
            self.projector.AddHomogenizationVariableToBeTimeFiltered(var, alpha)

    def UpdateDatabase(self, HMin):

        if self.dimension == 3:
            self.bin_of_objects_homogenization.UpdateSearchDatabase()

        else:
            self.bin_of_objects_homogenization.UpdateSearchDatabaseAssignedSize(HMin)

    def ProjectFromParticles(self, recalculate_neigh = True):

        if self.coupling_type != 3:
            self.projector.InterpolateFromDEMMesh(self.particles_model_part, self.homogenization_model_part, self.bin_of_objects_homogenization)

        else:
            pass
            #self.projector.HomogenizeFromDEMMesh(self.particles_model_part, self.homogenization_model_part, self.meso_scale_length, self.shape_factor, recalculate_neigh)

