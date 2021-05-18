import time as timer
import os
import KratosMultiphysics as Kratos
from Kratos import Logger
from  KratosMultiphysics.deprecation_management import DeprecationManager

import KratosMultiphysics.DEMApplication as Dem
Logger.Print("Running under OpenMP........", label="DEM")
import KratosMultiphysics.StructuralMechanicsApplication as Structural
import KratosMultiphysics.DemStructuresCouplingApplication as DemFem
from KratosMultiphysics.DemStructuresCouplingApplication.dem_fem_coupling_algorithm import Algorithm

class SPAlgorithm(Algorithm):

    def __init__(self):
        super(SPAlgorithm,self).__init__()

        sp_parameters_file_name = "SPParameters.json"

        with open(sp_parameters_file_name,'r') as parameter_file:
            self.sp_parameters = Kratos.Parameters(parameter_file.read())

        self.ValidateSettings()

        self.test_number = self.sp_parameters["problem_data"]["test_number"].GetInt()
        # Test types (4 different options):
        # Test number 0: no test simulation
        # Test number 1: CTW16 specimen
        # Test number 2: CTW10 specimen
        # Test number 3: Blind test specimen

        self.post_process_step_count = 0
        self.post_process_frequency = self.sp_parameters["post_process_tool"]["output_interval"].GetInt()
        self.use_post_process_tool = self.sp_parameters["post_process_tool"]["use_post_process_tool"].GetBool()
        if not self.use_post_process_tool:
            self.post_process_frequency = 0
        self.post_process_write_count = self.post_process_frequency

    @classmethod
    def GetDefaultParameters(cls):
        """This function returns the default-settings used by this class
        """
        return Kratos.Parameters("""{
            "problem_data"     : {
                "test_number"  : 1,
                "center" : [0.0,0.0,0.0],
                "axis"    : [0.0,0.0,1.0],
                "sp_radius": 0.042
            },
            "post_process_tool":{
                "use_post_process_tool": false,
                "output_interval": 0
            },
            "multiaxial_control_module_fem_dem_generalized_2d_utility" : {
                "Parameters"    : {
                    "control_module_delta_time": 2.0e-6,
                    "perturbation_tolerance": 1.0e-3,
                    "perturbation_period": 10,
                    "max_reaction_rate_factor": 10.0,
                    "stiffness_averaging_time_interval": 2.0e-6,
                    "velocity_averaging_time_interval": 2.0e-6,
                    "reaction_averaging_time_interval": 2.0e-8,
                    "output_interval": 0
                },
                "list_of_actuators" : []
            }
        }""")

    # This function can be extended with new deprecated variables as they are generated
    def TranslateLegacyVariablesAccordingToCurrentStandard(self, settings):
        # Defining a string to help the user understand where the warnings come from (in case any is thrown)
        context_string = type(self).__name__

        if settings.Has('post_process_tool'):
            sub_settings_where_var_is = settings['post_process_tool']
            old_name = 'output_frequency'
            new_name = 'output_interval'

            if DeprecationManager.HasDeprecatedVariable(context_string, sub_settings_where_var_is, old_name, new_name):
                DeprecationManager.ReplaceDeprecatedVariableName(sub_settings_where_var_is, old_name, new_name)


    def ValidateSettings(self):
        """This function validates the settings of the solver
        """
        default_settings = self.GetDefaultParameters()
        self.TranslateLegacyVariablesAccordingToCurrentStandard(self.sp_parameters)
        self.sp_parameters.ValidateAndAssignDefaults(default_settings)

    def Initialize(self):
        super(SPAlgorithm,self).Initialize()

        self.InitializeAdditionalProcessInfoVars()

        from KratosMultiphysics.DemStructuresCouplingApplication.multiaxial_control_module_fem_dem_generalized_2d_utility import MultiaxialControlModuleFEMDEMGeneralized2DUtility
        self.multiaxial_control_module = MultiaxialControlModuleFEMDEMGeneralized2DUtility(self.model, self.sp_parameters)
        self.multiaxial_control_module.ExecuteInitialize()
        # if self.test_number:
        #     from KratosMultiphysics.DemStructuresCouplingApplication.control_module_fem_dem_utility import ControlModuleFemDemUtility
        #     self.control_module_fem_dem_utility = ControlModuleFemDemUtility(self.model, self.dem_solution.spheres_model_part, self.test_number)
        #     self.control_module_fem_dem_utility.ExecuteInitialize()

        self.CreateSPMeasuringRingSubmodelpart()

        # Create Postprocess tool for SP
        if self.use_post_process_tool:
            from KratosMultiphysics.DemStructuresCouplingApplication.sand_production_post_process_tool import SandProductionPostProcessTool
            self.sp_post_process_tool = SandProductionPostProcessTool(self.structural_solution._GetSolver().GetComputingModelPart(),
                                                                        self.dem_solution.spheres_model_part,
                                                                        self.test_number)

        from KratosMultiphysics.DemStructuresCouplingApplication import stress_failure_check_utility
        self.stress_failure_check_utility = stress_failure_check_utility.StressFailureCheckUtility(self.dem_solution.spheres_model_part, self.test_number)

    def CreateSPMeasuringRingSubmodelpart(self):

        if not self.dem_solution.spheres_model_part.HasSubModelPart("RingSubmodelPart"):
            self.dem_solution.spheres_model_part.CreateSubModelPart('RingSubmodelPart')
        self.ring_submodelpart = self.dem_solution.spheres_model_part.GetSubModelPart('RingSubmodelPart')

        zone_radius_to_measure_2d_sp = self.sp_parameters["problem_data"]["sp_radius"].GetDouble()
        nodes_in_zone_radius_list = []
        elements_in_zone_radius_list = []

        for element in self.dem_solution.spheres_model_part.Elements:
            node = element.GetNode(0)
            x = node.X
            y = node.Y

            if (x * x + y * y) < zone_radius_to_measure_2d_sp * zone_radius_to_measure_2d_sp:
                nodes_in_zone_radius_list.append(node.Id)
                elements_in_zone_radius_list.append(element.Id)

        self.ring_submodelpart.AddNodes(nodes_in_zone_radius_list)
        self.ring_submodelpart.AddElements(elements_in_zone_radius_list)

    def InitializeAdditionalProcessInfoVars(self):
        self.dem_solution.spheres_model_part.ProcessInfo.SetValue(Dem.SIGMA_3_AVERAGE, 0.0)

    def _TransferStructuresSkinToDem(self):
        self.structural_mp = self.structural_solution._GetSolver().GetComputingModelPart()
        self.skin_mp = self.structural_mp.GetSubModelPart("DetectedByProcessSkinModelPart")
        # dem_walls_mp = self.dem_solution.rigid_face_model_part.CreateSubModelPart("SkinTransferredFromStructure")
        dem_walls_mp = self.dem_solution.rigid_face_model_part
        max_prop_id = 0
        for prop in dem_walls_mp.Properties:
            if prop.Id > max_prop_id:
                max_prop_id = prop.Id
        props = Kratos.Properties(max_prop_id + 1)
        # NOTE: this should be more general
        props[Dem.STATIC_FRICTION] = 0.2
        props[Dem.DYNAMIC_FRICTION] = 0.2
        props[Dem.WALL_COHESION] = 0.0
        props[Dem.COMPUTE_WEAR] = False
        props[Dem.SEVERITY_OF_WEAR] = 0.001
        props[Dem.IMPACT_WEAR_SEVERITY] = 0.001
        props[Dem.BRINELL_HARDNESS] = 200.0
        props[Kratos.YOUNG_MODULUS] = 7e9
        props[Kratos.POISSON_RATIO] = 0.16
        dem_walls_mp.AddProperties(props)
        DemFem.DemStructuresCouplingUtilities().TransferStructuresSkinToDem(self.skin_mp, dem_walls_mp, props)

    def RunSolutionLoop(self):

        self.dem_solution.step = 0
        self.dem_solution.time = 0.0
        self.dem_solution.time_old_print = 0.0
        self.time_dem = 0.0
        self.Dt_structural = self.structural_solution._GetSolver().settings["time_stepping"]["time_step"].GetDouble()

        while self.structural_solution.time < self.structural_solution.end_time:

            portion_of_the_force_which_is_new = 0.1
            DemFem.DemStructuresCouplingUtilities().SmoothLoadTrasferredToFem(self.dem_solution.rigid_face_model_part, portion_of_the_force_which_is_new)

            self.structural_solution.time = self.structural_solution._GetSolver().AdvanceInTime(self.structural_solution.time)

            # if self.test_number:
            #     self.control_module_fem_dem_utility.ExecuteInitializeSolutionStep()
            self.multiaxial_control_module.ExecuteInitializeSolutionStep()
            self.structural_solution.InitializeSolutionStep()
            self.structural_solution._GetSolver().Predict()
            self.structural_solution._GetSolver().SolveSolutionStep()

            time_final_DEM_substepping = self.structural_solution.time

            self.Dt_DEM = self.dem_solution.spheres_model_part.ProcessInfo.GetValue(Kratos.DELTA_TIME)

            DemFem.InterpolateStructuralSolutionForDEM().SaveStructuralSolution(self.structural_mp)

            DemFem.ComputeDEMFaceLoadUtility().ClearDEMFaceLoads(self.skin_mp)

            for self.dem_solution.time_dem in self.yield_DEM_time(self.dem_solution.time, time_final_DEM_substepping, self.Dt_DEM):
                self.dem_solution.time = self.dem_solution.time + self.dem_solution._GetSolver().dt

                self.dem_solution.step += 1

                self.dem_solution.DEMFEMProcedures.UpdateTimeInModelParts(self.dem_solution.all_model_parts, self.dem_solution.time, self.dem_solution._GetSolver().dt, self.dem_solution.step)

                self.dem_solution.InitializeSolutionStep()

                self.dem_solution._GetSolver().Predict()

                DemFem.InterpolateStructuralSolutionForDEM().InterpolateStructuralSolution(self.structural_mp, self.Dt_structural, self.structural_solution.time, self.dem_solution._GetSolver().dt, self.dem_solution.time)

                self.dem_solution.SolverSolve()

                self.dem_solution.FinalizeSolutionStep()

                DemFem.ComputeDEMFaceLoadUtility().CalculateDEMFaceLoads(self.skin_mp, self.dem_solution._GetSolver().dt, self.Dt_structural)

                #### PRINTING GRAPHS ####
                os.chdir(self.dem_solution.graphs_path)
                self.dem_solution.post_utils.ComputeMeanVelocitiesInTrap("Average_Velocity.txt", self.dem_solution.time, self.dem_solution.graphs_path)

                self.dem_solution.materialTest.MeasureForcesAndPressure()
                self.dem_solution.materialTest.PrintGraph(self.dem_solution.time)

                self.dem_solution.DEMFEMProcedures.PrintGraph(self.dem_solution.time)
                self.dem_solution.DEMFEMProcedures.PrintBallsGraph(self.dem_solution.time)

                self.dem_solution.DEMEnergyCalculator.CalculateEnergyAndPlot(self.dem_solution.time)

                self.dem_solution.BeforePrintingOperations(self.dem_solution.time)

                #### GiD IO ##########################################
                if self.dem_solution.IsTimeToPrintPostProcess():
                    self.dem_solution._GetSolver().PrepareElementsForPrinting()
                    if self.dem_solution.DEM_parameters["ContactMeshOption"].GetBool():
                        self.dem_solution._GetSolver().PrepareContactElementsForPrinting()
                    self.dem_solution.PrintResultsForGid(self.dem_solution.time)
                    self.SPPostProcessResults(self.dem_solution.time)
                    self.dem_solution.demio.PrintMultifileLists(self.dem_solution.time, self.dem_solution.post_path)
                    self.dem_solution.time_old_print = self.dem_solution.time

                    if self.test_number:
                        self.stress_failure_check_utility.ExecuteFinalizeSolutionStep()

            DemFem.InterpolateStructuralSolutionForDEM().RestoreStructuralSolution(self.structural_mp)

            # if self.test_number:
            #     self.control_module_fem_dem_utility.ExecuteFinalizeSolutionStep()
            self.multiaxial_control_module.ExecuteFinalizeSolutionStep()
            self.structural_solution.FinalizeSolutionStep()
            self.structural_solution.OutputSolutionStep()

            # Write SP data
            if self.IsPostProcessWriteStep():
                self.sp_post_process_tool.WriteData()

    def IsPostProcessWriteStep(self):
        self.post_process_step_count += 1
        if self.post_process_step_count == self.post_process_write_count:
            self.post_process_write_count += self.post_process_frequency
            return True
        else:
            return False

    def SPPostProcessResults(self,time):
        DemFem.DemStructuresCouplingUtilities().MarkBrokenSpheres(self.ring_submodelpart)

        center = Kratos.Array3()
        center[0] = self.sp_parameters["problem_data"]["center"][0].GetDouble()
        center[1] = self.sp_parameters["problem_data"]["center"][1].GetDouble()
        center[2] = self.sp_parameters["problem_data"]["center"][2].GetDouble()
        axis = Kratos.Array3()
        axis[0] = self.sp_parameters["problem_data"]["axis"][0].GetDouble()
        axis[1] = self.sp_parameters["problem_data"]["axis"][1].GetDouble()
        axis[2] = self.sp_parameters["problem_data"]["axis"][2].GetDouble()

        radius = 0
        if self.test_number == 1:
            radius = 0.0036195; #95% of the real hole. CTW16 specimen
        elif self.test_number == 2:
            radius = 0.012065; #95% of the real hole. CTW10 specimen
        elif self.test_number == 3:
            radius = 0.036195; #95% of the real hole. Blind Test

        self.dem_solution.creator_destructor.MarkParticlesForErasingGivenCylinder(self.ring_submodelpart, center, axis, radius)

        if self.test_number == 1 or self.test_number == 2:
            if self.structural_solution._GetSolver().main_model_part.ProcessInfo[Kratos.DOMAIN_SIZE] == 2:
                self.outer_walls_model_part = self.model["Structure.LinePressure2D_Outer_line"]
            else:
                self.outer_walls_model_part = self.model["Structure.SurfacePressure3D_lateral_pressure"]
            DemFem.DemStructuresCouplingUtilities().ComputeSandProductionWithDepthFirstSearchNonRecursiveImplementation(self.ring_submodelpart, self.outer_walls_model_part, time)
            DemFem.DemStructuresCouplingUtilities().ComputeSandProduction(self.ring_submodelpart, self.outer_walls_model_part, time)
        elif self.test_number == 3:
            if self.structural_solution._GetSolver().main_model_part.ProcessInfo[Kratos.DOMAIN_SIZE] == 2:
                self.outer_walls_model_part_1 = self.model["Structure.LinePressure2D_Left_line"]
                # self.outer_walls_model_part_2 = self.model["Structure.LinePressure2D_Bot_line"]
            else:
                self.outer_walls_model_part_1 = self.model["Structure.SurfacePressure3D_sigmaXpos"]
                # self.outer_walls_model_part_2 = self.model["Structure.SurfacePressure3D_sigmaYpos"]
            # NOTE: The stress printed in this case will also be the SigmaZ, but probably SigmaX is more appropriate
            DemFem.DemStructuresCouplingUtilities().ComputeSandProductionWithDepthFirstSearchNonRecursiveImplementation(self.ring_submodelpart, self.outer_walls_model_part_1, time)
            DemFem.DemStructuresCouplingUtilities().ComputeSandProduction(self.ring_submodelpart, self.outer_walls_model_part_1, time)
            # DemFem.DemStructuresCouplingUtilities().ComputeTriaxialSandProduction(self.ring_submodelpart, self.outer_walls_model_part_1, self.outer_walls_model_part_2, time)


if __name__ == "__main__":
    SPAlgorithm().Run()
