import KratosMultiphysics as Kratos
from KratosMultiphysics import Parameters
import KratosMultiphysics.SwimmingDEMApplication as SDEM
import KratosMultiphysics.DEMApplication as DEM
import sys
import numpy as np

class ProjectionModule:

    def __init__(self,
                fluid_model_part,
                balls_model_part,
                FEM_DEM_model_part,
                project_parameters,
                coupling_dem_vars,
                coupling_fluid_vars,
                time_filtered_vars,
                fluid_model_type,
                flow_field=None,
                domain_size=3):

        self.fluid_model_part = fluid_model_part
        self.particles_model_part = balls_model_part
        self.FEM_DEM_model_part = FEM_DEM_model_part
        self.project_parameters = project_parameters
        self.DEM_parameters = self.project_parameters["dem_parameters"]
        self.dimension = domain_size
        self.coupling_type = project_parameters["coupling"]["coupling_weighing_type"].GetInt()
        self.backward_coupling_parameters = project_parameters["coupling"]["backward_coupling"]
        self.meso_scale_length = self.backward_coupling_parameters["meso_scale_length"].GetDouble()
        self.shape_factor = self.backward_coupling_parameters["shape_factor"].GetDouble()
        self.do_impose_flow_from_field = project_parameters["custom_fluid"]["do_impose_flow_from_field_option"].GetBool()
        self.flow_field = flow_field
        self.use_drew_model = False
        self.use_new_implementation = True
        if fluid_model_type == "advmsDEM" or fluid_model_type == "aqsvmsDEM":
            self.use_drew_model = True

        if self.backward_coupling_parameters.Has("averaging_time_interval"):
            self.averaging_time_interval = self.backward_coupling_parameters["averaging_time_interval"].GetDouble()

         # Create projector_parameters
        self.projector_parameters = Parameters("{}")
        self.projector_parameters.AddValue("backward_coupling", project_parameters["coupling"]["backward_coupling"])
        self.projector_parameters.AddValue("coupling_type", project_parameters["coupling"]["coupling_weighing_type"])
        self.projector_parameters.AddValue("forward_coupling", project_parameters["coupling"]["forward_coupling"])
        self.projector_parameters.AddValue("viscosity_modification_type", project_parameters["coupling"]["backward_coupling"]["viscosity_modification_type"])
        self.projector_parameters.AddValue("n_particles_per_depth_distance", project_parameters["n_particles_in_depth"])
        self.projector_parameters.AddValue("body_force_per_unit_mass_variable_name", project_parameters["body_force_per_unit_mass_variable_name"])
        self.projector_parameters.AddValue("gentle_coupling_initiation", project_parameters["coupling"]["gentle_coupling_initiation"])
        self.search_strategy = DEM.OMP_DEMSearch()
        if "PeriodicDomainOption" in self.DEM_parameters.keys():
            if self.DEM_parameters["PeriodicDomainOption"].GetBool():
                self.search_strategy = DEM.OMP_DEMSearch(self.DEM_parameters["BoundingBoxMinX"].GetDouble(),
                                                         self.DEM_parameters["BoundingBoxMinY"].GetDouble(),
                                                         self.DEM_parameters["BoundingBoxMinZ"].GetDouble(),
                                                         self.DEM_parameters["BoundingBoxMaxX"].GetDouble(),
                                                         self.DEM_parameters["BoundingBoxMaxY"].GetDouble(),
                                                         self.DEM_parameters["BoundingBoxMaxZ"].GetDouble())


        if self.dimension == 3:

            if project_parameters["ElementType"].GetString() == "SwimmingNanoParticle":
                self.projector = SDEM.BinBasedNanoDEMFluidCoupledMapping3D(self.projector_parameters, self.search_strategy)

            else:
                self.projector = SDEM.BinBasedDEMFluidCoupledMapping3D(self.projector_parameters, self.search_strategy)
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
            averaging_time_interval = self.averaging_time_interval
            if self.averaging_time_interval < 15:
                alpha = 1 - np.exp(- averaging_time_interval)
            else:
                alpha = 1.0 / averaging_time_interval
            self.projector.AddFluidVariableToBeTimeFiltered(var, alpha)

    def UpdateHydrodynamicForces(self):
        self.projector.UpdateHydrodynamicForces(self.particles_model_part)

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

    def ProjectFromParticles(self, recalculate_neigh = True, update_impulse = True):
        if self.coupling_type != 3:
            self.projector.InterpolateFromDEMMesh(self.particles_model_part, self.fluid_model_part, self.bin_of_objects_fluid)
        else:
            self.projector.HomogenizeFromDEMMesh(self.particles_model_part, self.fluid_model_part, self.meso_scale_length, self.bin_of_objects_fluid, self.use_drew_model, self.use_new_implementation, recalculate_neigh, update_impulse)
            # self.projector.HomogenizeFromDEMMeshToElements(self.particles_model_part, self.fluid_model_part, self.meso_scale_length, self.bin_of_objects_fluid, self.use_drew_model, self.use_new_implementation, recalculate_neigh)

    def ComputePostProcessResults(self, particles_process_info):
        self.projector.ComputePostProcessResults(self.particles_model_part, self.fluid_model_part, self.FEM_DEM_model_part, self.bin_of_objects_fluid, particles_process_info)