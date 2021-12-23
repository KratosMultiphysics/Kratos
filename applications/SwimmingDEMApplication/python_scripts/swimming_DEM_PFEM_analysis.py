import KratosMultiphysics as Kratos
import KratosMultiphysics.DEMApplication as DEM
import KratosMultiphysics.SwimmingDEMApplication as SDEM
import KratosMultiphysics.SwimmingDEMApplication.swimming_DEM_procedures as SDP
import KratosMultiphysics.SwimmingDEMApplication.swimming_DEM_analysis as swimming_DEM_analysis
BaseAnalysis = swimming_DEM_analysis.SwimmingDEMAnalysis

class SDEMPFEMAnalysis(BaseAnalysis):

    def SetBetaParameters(self):
        super().SetBetaParameters()
        self.project_parameters["body_force_per_unit_mass_variable_name"].SetString('VOLUME_ACCELERATION')
        self.project_parameters["material_acceleration_calculation_type"].SetInt(8)

    def SetAllModelParts(self):
        self.all_model_parts = self._GetDEMAnalysis().all_model_parts

        # defining a fluid model
        self.all_model_parts.Add(self.fluid_model_part, "FluidPart")

        # defining a model part for the mixed part
        self.all_model_parts.Add(self.model.CreateModelPart("MixedPart"))

        self.mixed_model_part = self.all_model_parts.Get('MixedPart')

    def Initialize(self):
        super().Initialize()
        self.TransferWallsFromPfemToDem()

    def TransferWallsFromPfemToDem(self):
        destination_model_part = self._GetDEMAnalysis().rigid_face_model_part
        bodies_parts_list = self._GetFluidAnalysis().project_parameters["solver_settings"]["bodies_list"]
        for i in range(bodies_parts_list.size()):
            body_model_part_type = bodies_parts_list[i]["body_type"].GetString()
            if body_model_part_type == "Rigid":
                body_model_part_name = bodies_parts_list[i]["body_name"].GetString()
                source_model_part = self._GetFluidAnalysis().main_model_part.GetSubModelPart(body_model_part_name)
                SDEM.SwimmingDemInPfemUtils().TransferWalls(source_model_part, destination_model_part)

    def TransferGravityFromDisperseToFluid(self):
        # setting fluid's body force to the same as DEM's
        if self.project_parameters["custom_fluid"]["body_force_on_fluid_option"].GetBool():
            gravity = self.project_parameters["gravity_parameters"]["direction"].GetVector()
            gravity *= self.project_parameters["gravity_parameters"]["modulus"].GetDouble()
            Kratos.VariableUtils().SetVectorVar(Kratos.VOLUME_ACCELERATION, gravity, self.fluid_model_part.Nodes)

    def AssignKinematicViscosityFromDynamicViscosity(self):
        for node in self.fluid_model_part.Nodes:
            kinematic_viscosity = node.GetSolutionStepValue(Kratos.DYNAMIC_VISCOSITY,0) / node.GetSolutionStepValue(Kratos.DENSITY,0)
            node.SetSolutionStepValue(Kratos.VISCOSITY, 0, kinematic_viscosity)

    def FluidInitialize(self):
        self._GetFluidAnalysis().vars_man=self.vars_man
        self._GetFluidAnalysis().Initialize()
        bodies_parts_list = self._GetFluidAnalysis().project_parameters["solver_settings"]["bodies_list"]
        for i in range(bodies_parts_list.size()):
            body_model_part_type = bodies_parts_list[i]["body_type"].GetString()
            if body_model_part_type == "Fluid":
                body_model_part_name = bodies_parts_list[i]["body_name"].GetString()
                self.fluid_model_part = self._GetFluidAnalysis().main_model_part.GetSubModelPart(body_model_part_name)
                break

    def TransferTimeToFluidSolver(self):
        if self.step < self.GetFirstStepForFluidComputation() or self.stationarity:
            self._GetFluidAnalysis().time = self.time
            #self._GetFluidAnalysis().step = self.step #DO NOT INCREASE STEP IN PFEM, IT CRASHES (PROBABLY IT MUST DO SPECIAL THINGS FOR STEP=0)
            #self._GetFluidAnalysis().main_model_part.ProcessInfo[STEP] = self.step #DO NOT INCREASE STEP IN PFEM, IT CRASHES (PROBABLY IT MUST DO SPECIAL THINGS FOR STEP=0)
            Kratos.Logger.PrintInfo("SwimmingDEM", " [STEP:",self.step," TIME:",self.time,"]")

    def CloneTimeStep(self):
        self.TransferTimeToFluidSolver()

        if self.step < self.GetFirstStepForFluidComputation() or self.stationarity:
            self._GetFluidAnalysis().main_model_part.CloneTimeStep(self.time)

    def FluidSolve(self, time = 'None', solve_system = True):

        self._GetFluidAnalysis().InitializeSolutionStep()
        self._GetFluidAnalysis().SolveSolutionStep()
        self._GetFluidAnalysis().FinalizeSolutionStep()
        self.projection_module.UpdateDatabase(self.h_min)


    def GetFirstStepForFluidComputation(self):
        return 1

    def SetDragOutput(self):
        pass

    def FinalizeDragOutput(self):
        pass

    def SetPointGraphPrinter(self):
        pass

    def SetFluidSolverParameters(self):
        self.time           = self.pp.Start_time
        self.Dt             = self.pp.Dt
        self.out            = self.Dt
        if "REACTION" in self.pp.nodal_results:
            self.fluid_model_part.AddNodalSolutionStepVariable(Kratos.REACTION)
        if "DISTANCE" in self.pp.nodal_results:
            self.fluid_model_part.AddNodalSolutionStepVariable(Kratos.DISTANCE)
        self.vars_man.AddNodalVariables(self.fluid_model_part, self.vars_man.fluid_vars)
        if self.pp.type_of_inlet == 'ForceImposed':
            self.DEM_inlet = DEM.DEM_Force_Based_Inlet(self.dem_inlet_model_part, self.pp.force)

    def SetPostUtils(self):
        general_model_part = self._GetFluidAnalysis().main_model_part
        self.post_utils = SDP.PostUtils(self.swimming_DEM_gid_io,
                                        self.project_parameters,
                                        self.vars_man,
                                        general_model_part,
                                        self._GetDEMAnalysis().spheres_model_part,
                                        self._GetDEMAnalysis().cluster_model_part,
                                        self._GetDEMAnalysis().rigid_face_model_part,
                                        self.mixed_model_part)
    def SetEmbeddedTools(self):
        pass

    def PerformEmbeddedOperations(self):
        pass

    def _GetFluidAnalysis(self):
        if not hasattr(self, '_fluid_phase_analysis'):
            import KratosMultiphysics.SwimmingDEMApplication.DEM_coupled_pfem_fluid_dynamics_analysis as fluid_analysis
            self._fluid_phase_analysis = fluid_analysis.DEMCoupledPFEMFluidDynamicsAnalysis(self.model, self.project_parameters, self.vars_man)
            self._fluid_phase_analysis.main_path = self.main_path
        return self._fluid_phase_analysis
