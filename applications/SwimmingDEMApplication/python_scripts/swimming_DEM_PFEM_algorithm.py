from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import os
import sys
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.SwimmingDEMApplication import *
import swimming_DEM_analysis
import swimming_DEM_procedures as SDP
import variables_management as vars_man

sys.path.insert(0,'')
#import DEM_explicit_solver_var as DEM_parameters
BaseAlgorithm = swimming_DEM_analysis.Algorithm

class Algorithm(BaseAlgorithm):

    def SetFluidAlgorithm(self):
        import pfem_fluid_ready_for_dem_coupling as fluid_solution
        self.fluid_solution = fluid_solution.Solution(self.model)
        self.fluid_solution.main_path = self.main_path

    def SetCouplingParameters(self, varying_parameters):

        super(Algorithm,self).SetCouplingParameters(varying_parameters)
        self.pp.domain_size = self.fluid_solution.ProjectParameters["problem_data"]["dimension"].GetInt()

    def SetBetaParameters(self):

        self.pp.Dt = self.fluid_solution.GetDeltaTimeFromParameters()
        super(Algorithm,self).SetBetaParameters()
        self.pp.CFD_DEM["body_force_per_unit_mass_variable_name"].SetString('VOLUME_ACCELERATION')

    #def SetCouplingParameters(self, varying_parameters):
        #parameters_file = open("ProjectParametersDEM.json",'r')
        #self.pp.CFD_DEM = Parameters(parameters_file.read())
        #self.SetDoSolveDEMVariable()
        #self.pp.Dt = self.fluid_solution.GetDeltaTimeFromParameters()
        #self.SetBetaParameters()
        #self.SetCustomBetaParameters(varying_parameters)
        #self.pp.domain_size = self.fluid_solution.ProjectParameters["problem_data"]["domain_size"].GetInt()
        #super(Algorithm,self).SetCouplingParameters(varying_parameters)

    def SetAllModelParts(self):
        self.all_model_parts = self.disperse_phase_solution.all_model_parts

        # defining a fluid model
        self.all_model_parts.Add(self.fluid_model_part, "FluidPart")

        # defining a model part for the mixed part
        self.all_model_parts.Add(self.model.CreateModelPart("MixedPart"))

        self.mixed_model_part = self.all_model_parts.Get('MixedPart')

    def Initialize(self):
        super(Algorithm,self).Initialize()
        self.TransferWallsFromPfemToDem()

    def TransferWallsFromPfemToDem(self):
        destination_model_part = self.disperse_phase_solution.rigid_face_model_part
        bodies_parts_list = self.fluid_solution.ProjectParameters["solver_settings"]["bodies_list"]
        for i in range(bodies_parts_list.size()):
            body_model_part_type = bodies_parts_list[i]["body_type"].GetString()
            if body_model_part_type == "Rigid":
                body_model_part_name = bodies_parts_list[i]["body_name"].GetString()
                source_model_part = self.fluid_solution.main_model_part.GetSubModelPart(body_model_part_name)
                SwimmingDemInPfemUtils().TransferWalls(source_model_part, destination_model_part)

    def TransferGravityFromDisperseToFluid(self):
        # setting fluid's body force to the same as DEM's
        if self.pp.CFD_DEM["body_force_on_fluid_option"].GetBool():

            for node in self.fluid_model_part.Nodes:
                node.SetSolutionStepValue(VOLUME_ACCELERATION_X, 0, self.pp.CFD_DEM["GravityX"].GetDouble())
                node.SetSolutionStepValue(VOLUME_ACCELERATION_Y, 0, self.pp.CFD_DEM["GravityY"].GetDouble())
                node.SetSolutionStepValue(VOLUME_ACCELERATION_Z, 0, self.pp.CFD_DEM["GravityZ"].GetDouble())

    def AssignKinematicViscosityFromDynamicViscosity(self):
        for node in self.fluid_model_part.Nodes:
            kinematic_viscosity = node.GetSolutionStepValue(DYNAMIC_VISCOSITY,0) / node.GetSolutionStepValue(DENSITY,0)
            node.SetSolutionStepValue(VISCOSITY, 0, kinematic_viscosity)

    def FluidInitialize(self):

        self.fluid_solution.vars_man=vars_man
        self.fluid_solution.Initialize()
        bodies_parts_list = self.fluid_solution.ProjectParameters["solver_settings"]["bodies_list"]
        for i in range(bodies_parts_list.size()):
            body_model_part_type = bodies_parts_list[i]["body_type"].GetString()
            if body_model_part_type == "Fluid":
                body_model_part_name = bodies_parts_list[i]["body_name"].GetString()
                self.fluid_model_part = self.fluid_solution.main_model_part.GetSubModelPart(body_model_part_name)
                break

    def TransferTimeToFluidSolver(self):
        if self.step < self.GetFirstStepForFluidComputation() or self.stationarity:
            self.fluid_solution.time = self.time
            #self.fluid_solution.step = self.step #DO NOT INCREASE STEP IN PFEM, IT CRASHES (PROBABLY IT MUST DO SPECIAL THINGS FOR STEP=0)
            #self.fluid_solution.main_model_part.ProcessInfo[STEP] = self.step #DO NOT INCREASE STEP IN PFEM, IT CRASHES (PROBABLY IT MUST DO SPECIAL THINGS FOR STEP=0)
            print(" [STEP:",self.step," TIME:",self.time,"]")

    def CloneTimeStep(self):
        self.TransferTimeToFluidSolver()

        if self.step < self.GetFirstStepForFluidComputation() or self.stationarity:
            self.fluid_solution.main_model_part.CloneTimeStep(self.time)

    def FluidSolve(self, time = 'None', solve_system = True):

        self.fluid_solution.InitializeSolutionStep()
        self.fluid_solution.SolveSolutionStep()
        self.fluid_solution.FinalizeSolutionStep()
        self.projection_module.UpdateDatabase(self.h_min)


    def GetFirstStepForFluidComputation(self):
        return 1;

    def SetCutsOutput(self):
        pass

    def SetDragOutput(self):
        pass

    def FinalizeDragOutput(self):
        pass

    def SetPointGraphPrinter(self):
        pass

    def SetFluidSolverParameters(self):

        self.pp.domain_size
        self.time           = self.pp.Start_time
        self.Dt             = self.pp.Dt
        self.out            = self.Dt
        Nsteps         = self.pp.nsteps
        if "REACTION" in self.pp.nodal_results:
            self.fluid_model_part.AddNodalSolutionStepVariable(REACTION)
        if "DISTANCE" in self.pp.nodal_results:
            self.fluid_model_part.AddNodalSolutionStepVariable(DISTANCE)
        self.vars_man.AddNodalVariables(self.fluid_model_part, self.pp.fluid_vars)
        self.pp.variables_to_print_in_file
        if self.pp.type_of_inlet == 'ForceImposed':
            self.DEM_inlet = DEM_Force_Based_Inlet(self.DEM_inlet_model_part, self.pp.force)

    def SetPostUtils(self):
        general_model_part = self.fluid_solution.main_model_part
        self.post_utils = SDP.PostUtils(self.swimming_DEM_gid_io,
                                        self.pp,
                                        general_model_part,
                                        self.disperse_phase_solution.spheres_model_part,
                                        self.disperse_phase_solution.cluster_model_part,
                                        self.disperse_phase_solution.rigid_face_model_part,
                                        self.mixed_model_part)
    def SetEmbeddedTools(self):
        pass

    def PerformEmbeddedOperations(self):
        pass
