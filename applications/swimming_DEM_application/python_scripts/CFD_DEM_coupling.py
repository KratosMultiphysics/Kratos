from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import math
from KratosMultiphysics import *
#from KratosMultiphysics.IncompressibleFluidApplication import *
#from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.SwimmingDEMApplication import *
import sys

class ProjectionModule:

    def __init__(self,
                fluid_model_part,
                balls_model_part,
                FEM_DEM_model_part,
                pp,
                flow_field = None):

        self.fluid_model_part            = fluid_model_part
        self.particles_model_part        = balls_model_part
        self.FEM_DEM_model_part          = FEM_DEM_model_part
        self.pp                          = pp
        self.dimension                   = pp.domain_size
        self.coupling_type               = pp.CFD_DEM["coupling_weighing_type"].GetInt()
        self.meso_scale_length           = pp.CFD_DEM["meso_scale_length"].GetDouble()
        self.shape_factor                = pp.CFD_DEM["shape_factor"].GetDouble()
        self.do_impose_flow_from_field   = pp.CFD_DEM["do_impose_flow_from_field_option"].GetBool()
        self.flow_field                  = flow_field

        # Create projector_parameters
        self.projector_parameters = Parameters("{}")
        self.projector_parameters.AddValue("min_fluid_fraction",pp.CFD_DEM["min_fluid_fraction"])
        self.projector_parameters.AddValue("coupling_type",pp.CFD_DEM["coupling_weighing_type"])
        self.projector_parameters.AddValue("time_averaging_type",pp.CFD_DEM["time_averaging_type"])
        self.projector_parameters.AddValue("viscosity_modification_type",pp.CFD_DEM["viscosity_modification_type"])
        self.projector_parameters.AddValue("n_particles_per_depth_distance",pp.CFD_DEM["n_particles_in_depth"])
        self.projector_parameters.AddValue("body_force_per_unit_mass_variable_name",pp.CFD_DEM["body_force_per_unit_mass_variable_name"])

        if self.dimension == 3:

            if pp.CFD_DEM["ElementType"].GetString() == "SwimmingNanoParticle":
                self.projector = BinBasedNanoDEMFluidCoupledMapping3D(self.projector_parameters)

            else:
                self.projector = BinBasedDEMFluidCoupledMapping3D(self.projector_parameters)
            self.bin_of_objects_fluid = BinBasedFastPointLocator3D(fluid_model_part)

        else:
            if pp.CFD_DEM["ElementType"].GetString() == "SwimmingNanoParticle":
                self.projector = BinBasedNanoDEMFluidCoupledMapping2D(self.projector_parameters)

            else:
                self.projector = BinBasedDEMFluidCoupledMapping2D(self.projector_parameters)
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

    def ApplyForwardCouplingOfVelocityToSlipVelocityOnly(self):
        if self.do_impose_flow_from_field:
            self.ImposeVelocityOnDEMFromFieldToSlipVelocity()
        else:
            self.InterpolateVelocityOnSlipVelocity()

    def ProjectFromFluid(self, alpha):
        self.projector.InterpolateFromFluidMesh(self.fluid_model_part,
                                                self.particles_model_part,
                                                self.pp.CFD_DEM,
                                                self.bin_of_objects_fluid,
                                                alpha)

    def ProjectFromNewestFluid(self):
        self.projector.InterpolateFromNewestFluidMesh(self.fluid_model_part, self.particles_model_part, self.bin_of_objects_fluid)

    def InterpolateVelocityOnSlipVelocity(self):
        self.projector.InterpolateVelocityOnSlipVelocity(self.fluid_model_part, self.particles_model_part, self.bin_of_objects_fluid)

    def ImposeFluidFlowOnParticles(self):
        self.projector.ImposeFlowOnDEMFromField(self.flow_field, self.particles_model_part)

    def ImposeVelocityOnDEMFromFieldToSlipVelocity(self):
        self.projector.ImposeVelocityOnDEMFromFieldToSlipVelocity(self.flow_field, self.particles_model_part)

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
