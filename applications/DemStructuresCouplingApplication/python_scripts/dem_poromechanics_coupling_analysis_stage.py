import KratosMultiphysics as Kratos
from KratosMultiphysics.python_solver import PythonSolver
import KratosMultiphysics.analysis_stage

import KratosMultiphysics.DEMApplication as DEM
import KratosMultiphysics.PoromechanicsApplication as Poromechanics
import KratosMultiphysics.PoromechanicsApplication.poromechanics_analysis
import KratosMultiphysics.DemStructuresCouplingApplication as DemStructuresCouplingApplication
from sp_2d_rigid_fem_algorithm import DEMAnalysisStage2DSpRigidFem
from KratosMultiphysics.DEMApplication.DEM_analysis_stage import DEMAnalysisStage
import os
import math

class PoroMechanicsCouplingWithDemRadialMultiDofsControlModuleAnalysisStage(Kratos.analysis_stage.AnalysisStage):
    def __init__(self, model, parameters):

        self.parameters = parameters
        self.poromechanics_solution = Poromechanics.poromechanics_analysis.PoromechanicsAnalysis(model, parameters["poromechanics_parameters"])
        self.dem_solution = DEMAnalysisStage2DSpRigidFem(model, parameters["dem_parameters"])
        self.dem_solution.ReadModelParts = super(DEMAnalysisStage2DSpRigidFem, self.dem_solution).ReadModelParts
        super().__init__(model, parameters)

        self._CheckCoherentInputs()
        self.creator_destructor = DEM.ParticleCreatorDestructor()
        self.case_dir = os.getcwd()
        graph_folder_name = self.parameters["dem_parameters"]["problem_name"].GetString() + '_Graphs'
        graph_folder_path = os.path.join(self.case_dir, graph_folder_name)
        images_folder_name = 'dem_poro_coupling_images'
        self.images_folder_path = os.path.join(graph_folder_path, images_folder_name)
        if not os.path.isdir(self.images_folder_path):
            os.mkdir(self.images_folder_path)

        # General analysis settings
        self.number_of_DEM_steps_between_steadiness_checks = 25
        self.stationarity_measuring_tolerance = 5e1 #5e1 converges with GD 0.8 and CR 0.3 (case poro_dem_64b8)
        self.ignore_isolated_particles = False
        self.remove_isolated_particles = True
        self.breakable_material_after_first_dem_solution_steadiness_achieved = True
        self.use_precompression = True
        self.create_stress_ratio_animation = False
        # Casing settings
        self.add_casing = False
        self.casing_radius = 0.1
        self.casing_element_size = 0.003
        if not self.dem_solution.rigid_face_model_part.HasSubModelPart("Casing"):
            self.dem_solution.rigid_face_model_part.CreateSubModelPart('Casing')
        self.dem_solution.casing_submodelpart = self.dem_solution.rigid_face_model_part.GetSubModelPart('Casing')
        self.dem_solution.casing_submodelpart.SetValue(DEM.RIGID_BODY_OPTION, 1)

        file1 = os.path.join(os.getcwd(), "poro_solution_time_vs_sp_chunks.txt")
        file2 = os.path.join(os.getcwd(), "poro_solution_time_vs_sp_standard.txt")
        file3 = os.path.join(os.getcwd(), "poro_solution_time_vs_dem_time_vs_dems_steps.txt")
        file4 = os.path.join(os.getcwd(), "radial_normal_and_smoothed_reaction_stresses_values.txt")
        if os.path.exists(file1):
            os.remove(file1)
        if os.path.exists(file2):
            os.remove(file2)
        if os.path.exists(file3):
            os.remove(file3)
        if os.path.exists(file4):
            os.remove(file4)

        self.CreateInitialStressTensorFile()

    def FindSubModelPartWithExternalWalls(self):

        modelpart = self.model.GetModelPart('RigidFacePart')        
        submodelpart_with_external_walls = modelpart.GetSubModelPart('1')

        return submodelpart_with_external_walls

    def CreateInitialStressTensorFile(self):

        total_sigma_X = self.parameters["case_parameters"]["total_sigma_X"].GetDouble()
        total_sigma_Y = self.parameters["case_parameters"]["total_sigma_Y"].GetDouble()
        pore_pressure = self.parameters["case_parameters"]["pore_pressure"].GetDouble()
        biot_coefficient = self.parameters["case_parameters"]["biot_coefficient"].GetDouble()
        porosity = self.parameters["case_parameters"]["porosity"].GetDouble()
        permeability = self.parameters["case_parameters"]["permeability"].GetDouble()
        self.initial_effective_stress_X = -total_sigma_X + biot_coefficient * pore_pressure
        self.initial_effective_stress_Y = -total_sigma_Y + biot_coefficient * pore_pressure
        case_filename = self.parameters["poromechanics_parameters"]["problem_data"]["problem_name"].GetString()
        case_filename_plus_extension = case_filename + ".mdpa"
        file = open(case_filename_plus_extension, 'r')
        Xmin = Ymin = Xmax = Ymax = 0.0
        for line in file:
            if line.startswith("Begin Nodes"):
                while True:
                    nextline=next(file)
                    if nextline.startswith("End Nodes"):
                        break
                    data = nextline.split()
                    info = [float(data[1]), float(data[2])]
                    Xmin = info[0] if info[0] < Xmin else Xmin
                    Ymin = info[1] if info[1] < Ymin else Ymin
                    Xmax = info[0] if info[0] > Xmax else Xmax
                    Ymax = info[1] if info[1] > Ymax else Ymax                    
        file.close()
        initial_effective_sigmas_filename = self.parameters["poromechanics_parameters"]["problem_data"]["initial_stress_utility_settings"]["initial_input_filename"].GetString()
        initial_effective_sigmas_filename_plus_extension = initial_effective_sigmas_filename + ".mdpa"
        #TODO: print stresses in engineering format?
        effective_sigmas_datafile = open(initial_effective_sigmas_filename_plus_extension, 'w')
        effective_sigmas_datafile.write("Begin Properties 0\n")
        effective_sigmas_datafile.write("End Properties\n\n")
        effective_sigmas_datafile.write("Begin Nodes\n")
        effective_sigmas_datafile.write("1 "  + str(Xmin) + " "  + str(Ymin) + " 0\n")
        effective_sigmas_datafile.write("2  " + str(Xmax) + " "  + str(Ymin) + " 0\n")
        effective_sigmas_datafile.write("3  " + str(Xmax) + "  " + str(Ymax) + " 0\n")
        effective_sigmas_datafile.write("4 "  + str(Xmin) + "  " + str(Ymax) + " 0\n")
        effective_sigmas_datafile.write("End Nodes\n\n")
        effective_sigmas_datafile.write("Begin Elements Element2D3N\n")
        effective_sigmas_datafile.write("1 0 1 2 3\n")
        effective_sigmas_datafile.write("2 0 1 3 4\n")
        effective_sigmas_datafile.write("End Elements\n\n")
        effective_sigmas_datafile.write("Begin NodalData INITIAL_STRESS_TENSOR\n")
        effective_sigmas_datafile.write("1 0 [2,2] ((" + str(self.initial_effective_stress_X) +  ",0),(0," + str(self.initial_effective_stress_Y) + "))\n")
        effective_sigmas_datafile.write("2 0 [2,2] ((" + str(self.initial_effective_stress_X) +  ",0),(0," + str(self.initial_effective_stress_Y) + "))\n")
        effective_sigmas_datafile.write("3 0 [2,2] ((" + str(self.initial_effective_stress_X) +  ",0),(0," + str(self.initial_effective_stress_Y) + "))\n")
        effective_sigmas_datafile.write("4 0 [2,2] ((" + str(self.initial_effective_stress_X) +  ",0),(0," + str(self.initial_effective_stress_Y) + "))\n")        
        effective_sigmas_datafile.write("End NodalData\n")
        effective_sigmas_datafile.close()

    def Initialize(self):
        self.poromechanics_solution.Initialize()
        self.dem_solution.Initialize()
        super().Initialize()
        self.submodelpart_with_external_walls = self.FindSubModelPartWithExternalWalls()
        self.effective_stresses_communicator = DemStructuresCouplingApplication.EffectiveStressesCommunicatorUtility(self.poromechanics_solution._GetSolver().main_model_part, self.submodelpart_with_external_walls)
        self.effective_stresses_communicator.Initialize()
        self.pore_pressure_communicator_utility = DemStructuresCouplingApplication.PorePressureCommunicatorUtility(self.poromechanics_solution._GetSolver().main_model_part, self.dem_solution.spheres_model_part)
        self.pore_pressure_communicator_utility.Initialize()

        if self.use_precompression:
            self.ApplyPreCompression()

    def ApplyPreCompression(self):

        initial_effective_stress_matrix = Kratos.Matrix()
        initial_effective_stress_matrix.Resize(2,2)
        initial_effective_stress_matrix[0, 0] = self.initial_effective_stress_X
        initial_effective_stress_matrix[1, 0] = 0.0
        initial_effective_stress_matrix[0, 1] = 0.0
        initial_effective_stress_matrix[1, 1] = self.initial_effective_stress_Y

        self.effective_stresses_communicator.CommunicateGivenRadialEffectiveStressesToDemWalls(initial_effective_stress_matrix)
        self.dem_solution._GetSolver().cplusplus_strategy.BreakAllBonds()

        self.stationarity_checking_is_activated = False
        print("\n************************************ Solving DEM for the first time...\n", flush=True)
        self.DEM_steps_counter = 0
        self.dem_time_steps_per_fem_time_step = 0
        self.dem_solution.time = self.poromechanics_solution.time

        while not self.DEMSolutionIsSteady():
            self.dem_solution.time = self.dem_solution._GetSolver().AdvanceInTime(self.dem_solution.time)
            self.dem_solution.InitializeSolutionStep()
            self.dem_solution._GetSolver().Predict()
            self.dem_solution._GetSolver().SolveSolutionStep()
            self.dem_solution.FinalizeSolutionStep()
            self.AdditionalDEMOperationsToFinalizeDEMSolutionStep()
            self.dem_solution.OutputSolutionStep()

        self.dem_solution._GetSolver().cplusplus_strategy.HealAllBonds()
        DEM.ParallelBondUtilities().SetCurrentIndentationAsAReferenceInParallelBonds(self.dem_solution.spheres_model_part)

        if self.remove_isolated_particles:
            self.creator_destructor.MarkIsolatedParticlesForErasing(self.dem_solution.spheres_model_part)
        self.dem_solution.AdjustHoleInnerRadius()
        DEM.PreUtilities().ResetSkinParticles(self.dem_solution.spheres_model_part)
        self.dem_solution._GetSolver().cplusplus_strategy.ComputeSkin(self.dem_solution.spheres_model_part, 1.5)
        if self.add_casing:
            self.AddCasing()

    def AddCasing(self):

        max_FEM_node_Id = self.creator_destructor.FindMaxNodeIdInModelPart(self.dem_solution.rigid_face_model_part)
        max_FEM_elem_Id = self.creator_destructor.FindMaxConditionIdInModelPart(self.dem_solution.rigid_face_model_part)
        casing_length = 2.0 * math.pi * self.casing_radius
        number_of_casing_elems = casing_length / self.casing_element_size
        number_of_casing_elems = round(number_of_casing_elems)
        delta_angle_in_radians = 2.0 * math.pi / number_of_casing_elems
        for cond in self.submodelpart_with_external_walls.Conditions:
            props = cond.Properties
            break

        for i in range(number_of_casing_elems):
            Xi = self.casing_radius * math.cos(delta_angle_in_radians * i)
            Yi = self.casing_radius * math.sin(delta_angle_in_radians * i)
            self.dem_solution.casing_submodelpart.CreateNewNode(max_FEM_node_Id + i + 1, Xi, Yi, 0.0)
        for i in range(number_of_casing_elems - 1):
            self.dem_solution.casing_submodelpart.CreateNewCondition("RigidEdge2D2N", max_FEM_elem_Id + i + 1, [max_FEM_node_Id + i + 1, max_FEM_node_Id + i + 2], props)
        self.dem_solution.casing_submodelpart.CreateNewCondition("RigidEdge2D2N", max_FEM_elem_Id + number_of_casing_elems, [max_FEM_node_Id + number_of_casing_elems, max_FEM_node_Id + 1], props)

    def RunSolutionLoop(self):

        self.stationarity_checking_is_activated = False
        first_dem_solution_steadiness_achieved = False

        while self.poromechanics_solution.KeepAdvancingSolutionLoop():
            print("\n************************************ Solving FEM...\n", flush=True)
            self.poromechanics_solution.time = self.poromechanics_solution._GetSolver().AdvanceInTime(self.poromechanics_solution.time)
            self.poromechanics_solution.InitializeSolutionStep()
            self.poromechanics_solution._GetSolver().Predict()
            self.poromechanics_solution._GetSolver().SolveSolutionStep()
            self.poromechanics_solution.FinalizeSolutionStep()
            self.poromechanics_solution.OutputSolutionStep()

            self.effective_stresses_communicator.CopyWallCurrentEffectiveStressesToOldEffectiveStresses()
            self.effective_stresses_communicator.CommunicateCurrentRadialEffectiveStressesToDemWalls()
            self.pore_pressure_communicator_utility.ComputeForceOnParticlesDueToPorePressureGradient()

            print("\n************************************ Now solving DEM...\n", flush=True)
            self.DEM_steps_counter = 0
            self.dem_time_steps_per_fem_time_step = 0
            self.dem_solution.time = self.poromechanics_solution.time

            if first_dem_solution_steadiness_achieved and self.breakable_material_after_first_dem_solution_steadiness_achieved:
                for elem in self.dem_solution.spheres_model_part.Elements:
                    elem.Properties.GetSubProperties(elem.Properties.Id)[DEM.IS_UNBREAKABLE] = False
                    break

            while not self.DEMSolutionIsSteady():
                self.dem_solution.time = self.dem_solution._GetSolver().AdvanceInTime(self.dem_solution.time)
                self.center = KratosMultiphysics.Array3()
                self.center[0] = 0; # self.sp_parameters["problem_data"]["center"][0].GetDouble()
                self.center[1] = 0; # self.sp_parameters["problem_data"]["center"][1].GetDouble()
                self.center[2] = 0; # self.sp_parameters["problem_data"]["center"][2].GetDouble()
                self.axis = KratosMultiphysics.Array3()
                self.axis[0] = 0; # self.sp_parameters["problem_data"]["axis"][0].GetDouble()
                self.axis[1] = 0; # self.sp_parameters["problem_data"]["axis"][1].GetDouble()
                self.axis[2] = 1; # self.sp_parameters["problem_data"]["axis"][2].GetDouble()
                self.radius_to_delete_sp = 0.075
                self.creator_destructor.MarkParticlesForErasingGivenCylinder(self.dem_solution.spheres_model_part, self.center, self.axis, self.radius_to_delete_sp)

                self.dem_solution.InitializeSolutionStep()
                self.dem_solution._GetSolver().Predict()
                self.dem_solution._GetSolver().SolveSolutionStep()
                self.dem_solution.FinalizeSolutionStep()
                self.AdditionalDEMOperationsToFinalizeDEMSolutionStep()
                self.dem_solution.OutputSolutionStep()

            first_dem_solution_steadiness_achieved = True

    def AdditionalDEMOperationsToFinalizeDEMSolutionStep(self):

        if self.dem_solution.IsTimeToPrintPostProcess():

            os.chdir(self.images_folder_path)

            filename = 'radial_normal_and_smoothed_reaction_stresses_values_' + str(self.dem_solution.time) + '.txt'

            try:
                os.remove(filename)
            except OSError:
                pass

            radial_normal_and_smoothed_reaction_stresses_file = open(filename, 'a')

            average_ratio = 0.0

            for node in self.dem_solution.rigid_face_model_part.Nodes:
                total_radial_normal_stress_value = node.GetValue(DEM.RADIAL_NORMAL_STRESS_COMPONENT)
                A =  node.GetValue(DEM.SMOOTHED_REACTION_STRESS_X)
                B =  node.GetValue(DEM.SMOOTHED_REACTION_STRESS_Y)
                C =  node.GetValue(DEM.SMOOTHED_REACTION_STRESS_Z)
                total_smoothed_reaction_stress_value = math.sqrt(A*A + B*B + C*C)
                if abs(total_radial_normal_stress_value) > 0.0:
                    this_node_applied_stress_ratio = -total_smoothed_reaction_stress_value / total_radial_normal_stress_value
                else:
                    this_node_applied_stress_ratio = 0.0
                average_ratio += this_node_applied_stress_ratio
                displ = node.GetSolutionStepValue(Kratos.DISPLACEMENT)
                displ_norm = math.sqrt(displ[0]*displ[0]+displ[1]*displ[1]+displ[2]*displ[2])
                vel = node.GetSolutionStepValue(Kratos.VELOCITY)
                vel_norm = math.sqrt(vel[0]*vel[0]+vel[1]*vel[1]+vel[2]*vel[2])
                radial_normal_and_smoothed_reaction_stresses_file.write(str(node.Id) + " " + str(node.X) + " " + str(node.Y) + " " + str(node.Z) + " " + str(total_radial_normal_stress_value) + " " + str(total_smoothed_reaction_stress_value) + " " + str(displ_norm) + " " + str(vel_norm) + '\n')

            radial_normal_and_smoothed_reaction_stresses_file.close()
            average_ratio = average_ratio / self.dem_solution.rigid_face_model_part.NumberOfNodes()
            self.dem_solution.KratosPrintInfo("-------- Averaged loading ratio (average quotient current stress/ target stress) is  " + str(average_ratio) + " ------------ \n")

            os.chdir(self.case_dir)

            '''total_average_target_radial_normal_stress = 0.0
            total_forces_modulus = 0.0

            for node in self.dem_solution.rigid_face_model_part.Nodes:

                total_radial_normal_stress_value = node.GetValue(DEM.RADIAL_NORMAL_STRESS_COMPONENT)
                nodal_area = node.GetSolutionStepValue(DEM.DEM_NODAL_AREA)
                total_average_target_radial_normal_stress += nodal_area * total_radial_normal_stress_value

            for element in self.dem_solution.spheres_model_part.Elements:

                total_force_X = element.GetNode(0).GetSolutionStepValue(Kratos.EXTERNAL_APPLIED_FORCE_X)
                total_force_Y = element.GetNode(0).GetSolutionStepValue(Kratos.EXTERNAL_APPLIED_FORCE_Y)
                total_force_Z = element.GetNode(0).GetSolutionStepValue(Kratos.EXTERNAL_APPLIED_FORCE_Z)
                forces_modulus = math.sqrt(total_force_X * total_force_X + total_force_Y * total_force_Y + total_force_Z * total_force_Z)
                total_forces_modulus += forces_modulus'''

    def DEMSolutionIsSteady(self):

        self.dem_time_steps_per_fem_time_step += 1

        if self.DEM_steps_counter == self.number_of_DEM_steps_between_steadiness_checks:

            self.DEM_steps_counter = 0

            print("\n*** DEM stationarity will now be checked...", flush=True)

            if not DEM.StationarityChecker().CheckIfVariableIsNullInModelPart(self.dem_solution.spheres_model_part, Kratos.TOTAL_FORCES_X, self.stationarity_measuring_tolerance, self.ignore_isolated_particles):
                print("  F_X is larger than the desired maximum\n", flush=True)
                self.stationarity_checking_is_activated = True
                return False
            if not DEM.StationarityChecker().CheckIfVariableIsNullInModelPart(self.dem_solution.spheres_model_part, Kratos.TOTAL_FORCES_Y, self.stationarity_measuring_tolerance, self.ignore_isolated_particles):
                print("  F_Y is larger than the desired maximum\n", flush=True)
                self.stationarity_checking_is_activated = True
                return False
            if not DEM.StationarityChecker().CheckIfVariableIsNullInModelPart(self.dem_solution.spheres_model_part, Kratos.TOTAL_FORCES_Z, self.stationarity_measuring_tolerance, self.ignore_isolated_particles):
                print("  F_Z is larger than the desired maximum\n", flush=True)
                self.stationarity_checking_is_activated = True
                return False

            if self.stationarity_checking_is_activated:
                print("\n  ********** DEM solution is steady! ********** \n", flush=True)
                with open('poro_solution_time_vs_dem_time_vs_dems_steps.txt', 'a') as time_and_steps_file:
                    time_and_steps_file.write(str(self.poromechanics_solution.time) + " " + str(self.dem_solution.time) + " " + str(self.dem_time_steps_per_fem_time_step) + '\n')
                return True
            else:
                return False
        else:
            self.DEM_steps_counter += 1
            return False

    def Finalize(self):
        self.poromechanics_solution.Finalize()
        self.dem_solution.Finalize()
        super().Finalize()
        if self.create_stress_ratio_animation:
            self.CreateReactionStressToTargetRatioAnimation()

    def CreateReactionStressToTargetRatioAnimation(self):
        # Use 'pip install opencv-python' at the command line in Windows to obtain cv2
        import cv2
        import matplotlib.pyplot as plt
        print_one_out_of_n_results = 1
        frames_per_second = 5
        n = 0
        list_of_pngs = []
        os.chdir(self.images_folder_path)

        for filename in os.listdir(self.images_folder_path):
            n = n + 1

            if n % print_one_out_of_n_results:
                continue

            if filename.endswith(".txt"): 
                X, Y = [], []
                for line in open(filename, 'r'):
                    values = [float(s) for s in line.split()]
                    if values[4]:
                        X.append(values[1])
                        Y.append(-values[5]/values[4])

                plt.axis((-1.5,1.5,0,2))
                plt.xlabel("Xcoord (m)")
                plt.ylabel("radial effective stress to target ratio")
                plt.plot(X, Y, '.', color='blue')
                plt.grid()
                printfile = filename[:-4]
                printfilename = os.path.join(self.images_folder_path, printfile + '.png')
                list_of_pngs.append(printfilename)
                plt.savefig(printfilename)
                plt.close()

        video_name = 'video.avi'
        images = [img for img in os.listdir(self.images_folder_path) if img.endswith(".png")]
        frame = cv2.imread(os.path.join(self.images_folder_path, images[0]))
        height, width, layers = frame.shape
        video = cv2.VideoWriter(video_name, 0, frames_per_second, (width,height))

        for image in images:
            video.write(cv2.imread(os.path.join(self.images_folder_path, image)))

        cv2.destroyAllWindows()
        video.release()
        os.chdir(self.case_dir)

    def _CreateSolver(self):
        return DemPoroMechanicsCouplingSolver(self.poromechanics_solution._GetSolver(), self.dem_solution._GetSolver(), self.parameters)

    def _CheckCoherentInputs(self):
        variables_exported_from_gauss_points_to_nodes = self.parameters["poromechanics_parameters"]["solver_settings"]["gp_to_nodal_variable_list"].GetStringArray()
        if not "LIQUID_PRESSURE_GRADIENT" in variables_exported_from_gauss_points_to_nodes or not "EFFECTIVE_STRESS_TENSOR" in variables_exported_from_gauss_points_to_nodes:
            raise Exception("Coupling DEM with Poromechanics Error: [\"poromechanics_parameters\"][\"solver_settings\"][\"gp_to_nodal_variable_list\"] must contain \"LIQUID_PRESSURE_GRADIENT\" and \"EFFECTIVE_STRESS_TENSOR\"\n")

    def _YieldDEMTime(self, current_time, current_time_plus_increment, delta_time):
        current_time += delta_time
        tolerance = 0.0001
        while current_time < (current_time_plus_increment - tolerance * delta_time):
            yield current_time
            current_time += delta_time

        current_time = current_time_plus_increment
        yield current_time

class DemPoroMechanicsCouplingSolver(PythonSolver):
    def __init__(self, poromechanics_solver, dem_solver, parameters):
        self.poromechanics_solver = poromechanics_solver
        self.dem_solver = dem_solver

    def ImportModelPart(self):
        pass

    def GetComputingModelPart(self):
        return self.poromechanics_solver.GetComputingModelPart()

if __name__ == "__main__":
    with open("ProjectParameters.json", 'r') as parameter_file:
        project_parameters = Kratos.Parameters(parameter_file.read())

    model = Kratos.Model()
    PoroMechanicsCouplingWithDemRadialMultiDofsControlModuleAnalysisStage(model, project_parameters).Run()
