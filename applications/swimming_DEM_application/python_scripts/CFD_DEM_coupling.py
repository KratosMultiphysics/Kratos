from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import math
from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.SwimmingDEMApplication import *

class ProjectionModule:

    def __init__(self, fluid_model_part, balls_model_part, FEM_DEM_model_part, dimension, pp):

        self.fluid_model_part            = fluid_model_part
        self.particles_model_part        = balls_model_part
        self.FEM_DEM_model_part          = FEM_DEM_model_part
        self.dimension                   = dimension
        self.min_fluid_fraction          = pp.CFD_DEM.min_fluid_fraction
        self.coupling_type               = pp.CFD_DEM.coupling_weighing_type
        self.viscosity_modification_type = pp.CFD_DEM.viscosity_modification_type
        self.n_particles_in_depth        = pp.CFD_DEM.n_particles_in_depth
        self.meso_scale_length           = pp.CFD_DEM.meso_scale_length
        self.shape_factor                = pp.CFD_DEM.shape_factor

        if (self.dimension == 3):

            if pp.CFD_DEM.ElementType == "SwimmingNanoParticle":
                self.projector = BinBasedNanoDEMFluidCoupledMapping3D(self.min_fluid_fraction, self.coupling_type, self.viscosity_modification_type)

            else:
                self.projector = BinBasedDEMFluidCoupledMapping3D(self.min_fluid_fraction, self.coupling_type, self.viscosity_modification_type)
            self.bin_of_objects_fluid = BinBasedFastPointLocator3D(fluid_model_part)

        else:
            if pp.CFD_DEM.ElementType == "SwimmingNanoParticle":
                self.projector = BinBasedNanoDEMFluidCoupledMapping2D(self.min_fluid_fraction, self.coupling_type, self.viscosity_modification_type, self.n_particles_in_depth)

            else:
                self.projector = BinBasedDEMFluidCoupledMapping2D(self.min_fluid_fraction, self.coupling_type, self.viscosity_modification_type, self.n_particles_in_depth)
            self.bin_of_objects_fluid = BinBasedFastPointLocator2D(fluid_model_part)

        # telling the projector which variables we are interested in modifying

        for var in pp.coupling_dem_vars:
            self.projector.AddDEMCouplingVariable(var)

        for var in pp.coupling_fluid_vars:
            self.projector.AddFluidCouplingVariable(var)

        # calculating the fluid nodal areas that are needed for the coupling

        self.area_calculator = CalculateNodalAreaProcess(self.fluid_model_part, self.dimension)
        self.area_calculator.Execute()

    def UpdateDatabase(self, HMin):

        if (self.dimension == 3):
            self.bin_of_objects_fluid.UpdateSearchDatabase()

        else:
            self.bin_of_objects_fluid.UpdateSearchDatabaseAssignedSize(HMin)

    def ProjectFromFluid(self, alpha):
        self.projector.InterpolateFromFluidMesh(self.fluid_model_part, self.particles_model_part, self.bin_of_objects_fluid, alpha)

    def ProjectFromNewestFluid(self):
        self.projector.InterpolateFromNewestFluidMesh(self.fluid_model_part, self.particles_model_part, self.bin_of_objects_fluid)

    def ProjectFromParticles(self, recalculate_neigh = True):

        if (self.coupling_type != 3):
            self.projector.InterpolateFromDEMMesh(self.particles_model_part, self.fluid_model_part, self.bin_of_objects_fluid)

        else:
            self.projector.HomogenizeFromDEMMesh(self.particles_model_part, self.fluid_model_part, self.meso_scale_length, self.shape_factor, recalculate_neigh)

    def ComputePostProcessResults(self, particles_process_info):
        self.projector.ComputePostProcessResults(self.particles_model_part, self.fluid_model_part, self.FEM_DEM_model_part, self.bin_of_objects_fluid, particles_process_info)


