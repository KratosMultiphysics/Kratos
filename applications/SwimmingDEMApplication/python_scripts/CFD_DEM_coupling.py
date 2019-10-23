from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import KratosMultiphysics as Kratos
from KratosMultiphysics import Parameters
import KratosMultiphysics.SwimmingDEMApplication as SDEM
import sys

class ProjectionModule:

    def __init__(self,
                fluid_model_part,
                balls_model_part,
                FEM_DEM_model_part,
                project_parameters,
                coupling_dem_vars,
                coupling_fluid_vars,
                time_filtered_vars,
                flow_field=None,
                domain_size=3):

        self.fluid_model_part = fluid_model_part
        self.particles_model_part = balls_model_part
        self.FEM_DEM_model_part = FEM_DEM_model_part
        self.project_parameters = project_parameters
        self.dimension = domain_size
        self.coupling_type = project_parameters["coupling"]["coupling_weighing_type"].GetInt()
        self.backward_coupling_parameters = project_parameters["coupling"]["backward_coupling"]
        self.meso_scale_length = self.backward_coupling_parameters["meso_scale_length"].GetDouble()
        self.shape_factor = self.backward_coupling_parameters["shape_factor"].GetDouble()
        self.do_impose_flow_from_field = project_parameters["custom_fluid"]["do_impose_flow_from_field_option"].GetBool()
        self.flow_field = flow_field

        # Create projector_parameters
        self.projector_parameters = Parameters("{}")
        self.projector_parameters.AddValue("backward_coupling", project_parameters["coupling"]["backward_coupling"])
        self.projector_parameters.AddValue("coupling_type", project_parameters["coupling"]["coupling_weighing_type"])
        self.projector_parameters.AddValue("forward_coupling", project_parameters["coupling"]["forward_coupling"])
        self.projector_parameters.AddValue("viscosity_modification_type", project_parameters["coupling"]["backward_coupling"]["viscosity_modification_type"])
        self.projector_parameters.AddValue("n_particles_per_depth_distance", project_parameters["n_particles_in_depth"])
        self.projector_parameters.AddValue("body_force_per_unit_mass_variable_name", project_parameters["body_force_per_unit_mass_variable_name"])

        if self.dimension == 3:

            if project_parameters["ElementType"].GetString() == "SwimmingNanoParticle":
                self.projector = SDEM.BinBasedNanoDEMFluidCoupledMapping3D(self.projector_parameters)

            else:
                self.projector = SDEM.BinBasedDEMFluidCoupledMapping3D(self.projector_parameters)
            self.bin_of_objects_fluid = Kratos.BinBasedFastPointLocator3D(fluid_model_part)

        else:
            if project_parameters["ElementType"].GetString() == "SwimmingNanoParticle":
                self.projector = SDEM.BinBasedNanoDEMFluidCoupledMapping2D(self.projector_parameters)

            else:
                self.projector = SDEM.BinBasedDEMFluidCoupledMapping2D(self.projector_parameters)
            self.bin_of_objects_fluid = Kratos.BinBasedFastPointLocator2D(fluid_model_part)

        # telling the projector which variables we are interested in modifying

        for var in coupling_dem_vars:
            self.projector.AddDEMCouplingVariable(var)

        for var in coupling_fluid_vars:
            self.projector.AddFluidCouplingVariable(var)

        for var in coupling_dem_vars:
            if var in {Kratos.FLUID_VEL_PROJECTED, Kratos.FLUID_ACCEL_PROJECTED, Kratos.FLUID_VEL_LAPL_PROJECTED, Kratos.FLUID_ACCEL_FOLLOWING_PARTICLE_PROJECTED}:
                self.projector.AddDEMVariablesToImpose(var)
                coupling_dem_vars.remove(var)
            self.projector.AddDEMVariablesToImpose(Kratos.AUX_VEL)

        for var in time_filtered_vars:
            self.projector.AddFluidVariableToBeTimeFiltered(var, 0.004)

    def UpdateDatabase(self, HMin):

        if self.dimension == 3:
            self.bin_of_objects_fluid.UpdateSearchDatabase()

        else:
            self.bin_of_objects_fluid.UpdateSearchDatabaseAssignedSize(HMin)

    def ApplyForwardCoupling(self, alpha=None):
        if self.do_impose_flow_from_field:
            self.ImposeFluidFlowOnParticles()
        else:
            if alpha == None:
                alpha = 1.0
            self.ProjectFromFluid(alpha)

    def ApplyForwardCouplingOfVelocityToAuxVelocityOnly(self, alpha=None):
        if self.do_impose_flow_from_field:
            self.ImposeVelocityOnDEMFromFieldToAuxVelocity()
        elif alpha == None:
            alpha = 1.0
        self.InterpolateVelocityOnAuxVelocity(alpha)

    def ProjectFromFluid(self, alpha):

        self.projector.InterpolateFromFluidMesh(self.fluid_model_part,
                                                self.particles_model_part,
                                                self.project_parameters,
                                                self.bin_of_objects_fluid,
                                                alpha)

    def InterpolateVelocityOnAuxVelocity(self, alpha):
        self.projector.InterpolateVelocityOnAuxVelocity(self.fluid_model_part, self.particles_model_part, self.bin_of_objects_fluid, alpha)

    def ImposeFluidFlowOnParticles(self):
        self.projector.ImposeFlowOnDEMFromField(self.flow_field, self.particles_model_part)

    def ImposeVelocityOnDEMFromFieldToAuxVelocity(self):
        self.projector.ImposeVelocityOnDEMFromFieldToAuxVelocity(self.flow_field, self.particles_model_part)

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
