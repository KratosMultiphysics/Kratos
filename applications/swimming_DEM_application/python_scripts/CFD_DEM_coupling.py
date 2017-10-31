from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import math
from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
#from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.SwimmingDEMApplication import *
import sys

class ProjectionModule:

    def __init__(self, fluid_model_part, balls_model_part, FEM_DEM_model_part, dimension, pp, flow_field = None):

        self.fluid_model_part            = fluid_model_part
        self.particles_model_part        = balls_model_part
        self.FEM_DEM_model_part          = FEM_DEM_model_part
        self.dimension                   = dimension
        self.min_fluid_fraction          = pp.CFD_DEM["min_fluid_fraction"].GetDouble()
        self.coupling_type               = pp.CFD_DEM["coupling_weighing_type"].GetInt()
        self.time_averaging_type         = pp.CFD_DEM["time_averaging_type"].GetInt()
        self.viscosity_modification_type = pp.CFD_DEM["viscosity_modification_type"].GetInt()
        self.n_particles_in_depth        = pp.CFD_DEM.n_particles_in_depth
        self.meso_scale_length           = pp.CFD_DEM["meso_scale_length"].GetDouble()
        self.shape_factor                = pp.CFD_DEM["shape_factor"].GetDouble()
        self.do_impose_flow_from_field   = pp.CFD_DEM["do_impose_flow_from_field_option"].GetBool()
        self.flow_field                  = flow_field

        if self.dimension == 3:

            if pp.CFD_DEM["ElementType"].GetString() == "SwimmingNanoParticle":
                self.projector = BinBasedNanoDEMFluidCoupledMapping3D(self.min_fluid_fraction, self.coupling_type, self.time_averaging_type , self.viscosity_modification_type)

            else:
                self.projector = BinBasedDEMFluidCoupledMapping3D(self.min_fluid_fraction, self.coupling_type, self.time_averaging_type , self.viscosity_modification_type)
            self.bin_of_objects_fluid = BinBasedFastPointLocator3D(fluid_model_part)

        else:
            if pp.CFD_DEM["ElementType"].GetString() == "SwimmingNanoParticle":
                self.projector = BinBasedNanoDEMFluidCoupledMapping2D(self.min_fluid_fraction, self.coupling_type, self.time_averaging_type , self.viscosity_modification_type, self.n_particles_in_depth)

            else:
                self.projector = BinBasedDEMFluidCoupledMapping2D(self.min_fluid_fraction, self.coupling_type, self.time_averaging_type , self.viscosity_modification_type, self.n_particles_in_depth)
            self.bin_of_objects_fluid = BinBasedFastPointLocator2D(fluid_model_part)

        # telling the projector which variables we are interested in modifying

        for var in pp.coupling_dem_vars:
            self.projector.AddDEMCouplingVariable(var)

        for var in pp.coupling_fluid_vars:
            self.projector.AddFluidCouplingVariable(var)

        for var in pp.coupling_dem_vars:
            if var in {FLUID_VEL_PROJECTED, FLUID_ACCEL_PROJECTED, FLUID_VEL_LAPL_PROJECTED, FLUID_ACCEL_FOLLOWING_PARTICLE_PROJECTED}:
                self.projector.AddDEMVariablesToImpose(var)
                pp.coupling_dem_vars.remove(var)
            self.projector.AddDEMVariablesToImpose(SLIP_VELOCITY)

        for var in pp.time_filtered_vars:
            self.projector.AddFluidVariableToBeTimeFiltered(var, 0.004)

        # calculating the fluid nodal areas that are needed for the coupling

        self.area_calculator = CalculateNodalAreaProcess(self.fluid_model_part, self.dimension)
        self.area_calculator.Execute()

    def UpdateDatabase(self, HMin):

        if self.dimension == 3:
            self.bin_of_objects_fluid.UpdateSearchDatabase()

        else:
            self.bin_of_objects_fluid.UpdateSearchDatabaseAssignedSize(HMin)

    def ApplyForwardCoupling(self, alpha = None):
        if self.do_impose_flow_from_field:
            self.ImposeFluidFlowOnParticles()
        else:
            if alpha == None:
                self.ProjectFromNewestFluid()
            else:
                self.ProjectFromFluid(alpha)

    def ApplyForwardCouplingOfVelocityOnly(self):
        if self.do_impose_flow_from_field:
            self.ImposeVelocityOnDEMFromField()
        else:
            self.InterpolateVelocity()

    def ProjectFromFluid(self, alpha):
        self.projector.InterpolateFromFluidMesh(self.fluid_model_part, self.particles_model_part, self.bin_of_objects_fluid, alpha)

    def ProjectFromNewestFluid(self):
        self.projector.InterpolateFromNewestFluidMesh(self.fluid_model_part, self.particles_model_part, self.bin_of_objects_fluid)

    def InterpolateVelocity(self):
        self.projector.InterpolateVelocity(self.fluid_model_part, self.particles_model_part, self.bin_of_objects_fluid)

    def ImposeFluidFlowOnParticles(self):
        self.projector.ImposeFlowOnDEMFromField(self.flow_field, self.particles_model_part)

    def ImposeVelocityOnDEMFromField(self):
        self.projector.ImposeVelocityOnDEMFromField(self.flow_field, self.particles_model_part)

    def ProjectFromParticles(self, recalculate_neigh = True):
        #print("\nProjecting from particles to the fluid...")
        #sys.stdout.flush()

        if self.coupling_type != 3:
            self.projector.InterpolateFromDEMMesh(self.particles_model_part, self.fluid_model_part, self.bin_of_objects_fluid)

        else:
            self.projector.HomogenizeFromDEMMesh(self.particles_model_part, self.fluid_model_part, self.meso_scale_length, self.shape_factor, recalculate_neigh)

        #print("\nFinished projecting from particles to the fluid...")
        #sys.stdout.flush()

    def ComputePostProcessResults(self, particles_process_info):
        self.projector.ComputePostProcessResults(self.particles_model_part, self.fluid_model_part, self.FEM_DEM_model_part, self.bin_of_objects_fluid, particles_process_info)
