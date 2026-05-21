import KratosMultiphysics
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
import KratosMultiphysics.DEMApplication as DEM
from KratosMultiphysics.DEMApplication.DEM_analysis_stage import DEMAnalysisStage
from KratosMultiphysics.DEMApplication import DEM_procedures as DEM_procedures
import math
import datetime

class DecompressedMaterialTriaxialTest(DEMAnalysisStage):

    def __init__(self, model, parameters):
        super().__init__(model, parameters)

        self.parameters = parameters
        self.compression_stage_completed = False
        self.decompression_stage_completed = False
        # Units in Pa and m/s respectively
        self.SigmaHorizontal = self.parameters["material_test_settings"]["SigmaHorizontal"].GetDouble()
        self.SigmaVertical = self.parameters["material_test_settings"]["SigmaVertical"].GetDouble()
        self.SigmaVerticalAlmostZero = self.parameters["material_test_settings"]["SigmaVerticalAlmostZero"].GetDouble()
        self.PlatesDecompressionVelocity = self.parameters["material_test_settings"]["PlatesDecompressionVelocity"].GetDouble()

    def Initialize(self):
        super().Initialize()
        self.GetMainProblemParameters()
        self.InitializeMaterialTest()
        self.PrepareDataForGraph()
        self._GetSolver().cplusplus_strategy.BreakAllBonds()
        self.ApplyPrecompression()
        self._GetSolver().cplusplus_strategy.HealAllBonds()
        ParallelBondUtilities().SetCurrentIndentationAsAReferenceInParallelBonds(self.spheres_model_part)
        PreUtilities().ResetSkinParticles(self.spheres_model_part)
        self._GetSolver().cplusplus_strategy.ComputeSkin(self.spheres_model_part, 1.5)
        self.ApplyDecompression()
        self.RestoreLoadingVelocity()

    def GetMainProblemParameters(self):
        list_of_material_relations = self.DEM_material_parameters["material_relations"]        
        for material_relation in list_of_material_relations:
            contact_properties = material_relation["Variables"]
        self.loose_young = contact_properties["LOOSE_MATERIAL_YOUNG_MODULUS"].GetDouble()
        self.bonded_young = contact_properties["BONDED_MATERIAL_YOUNG_MODULUS"].GetDouble()
        self.friction = contact_properties["STATIC_FRICTION"].GetDouble()
        self.contact_tau = contact_properties["CONTACT_TAU_ZERO"].GetDouble()

    def ApplyPrecompression(self):

        print("\n************************************ Applying Precompression...\n", flush=True)
        while not self.compression_stage_completed:
            self.time = self._GetSolver().AdvanceInTime(self.time)
            self.InitializeSolutionStep()
            self._GetSolver().Predict()
            self._GetSolver().SolveSolutionStep()
            self.FinalizeSolutionStepPrecompression()
            self.OutputSolutionStep()
        print("\n*************************** Finished Applying Precompression!!!\n", flush=True)

    def ResetLoadingVelocity(self):
        for smp in self.rigid_face_model_part.SubModelParts:
            if smp[IDENTIFIER] == 'TOP':
                smp[LINEAR_VELOCITY_Y] = 0.5 * self.PlatesDecompressionVelocity
            if smp[IDENTIFIER] == 'BOTTOM':
                smp[LINEAR_VELOCITY_Y] = -0.5 * self.PlatesDecompressionVelocity

    def RestoreLoadingVelocity(self):
        for smp in self.rigid_face_model_part.SubModelParts:
            if smp[IDENTIFIER] == 'TOP':
                smp[LINEAR_VELOCITY_Y] =  0.5 * self.LoadingVelocity
            if smp[IDENTIFIER] == 'BOTTOM':
                smp[LINEAR_VELOCITY_Y] = -0.5 * self.LoadingVelocity

    def ApplyDecompression(self):

        print("\n************************************ Applying Decompression...\n", flush=True)
        while not self.decompression_stage_completed:
            self.time = self._GetSolver().AdvanceInTime(self.time)
            self.InitializeSolutionStep()
            self._GetSolver().Predict()
            self.ResetLoadingVelocity()
            self._GetSolver().SolveSolutionStep()
            self.FinalizeSolutionStepDecompression()
            self.OutputSolutionStep()
        print("\n*************************** Finished Applying Decompression!!!\n", flush=True)

    def RunSolutionLoop(self):

        print("\n************************************ Applying standard triaxial...\n", flush=True)
        while self.KeepAdvancingSolutionLoop():
            self.time = self._GetSolver().AdvanceInTime(self.time)
            self.InitializeSolutionStep()
            self._GetSolver().Predict()
            self._GetSolver().SolveSolutionStep()
            self.FinalizeSolutionStep()
            self.OutputSolutionStep()
        print("\n*************************** Finished Applying standard triaxial...\n", flush=True)

    def OutputSolutionStep(self):
        super().OutputSolutionStep()
        self.PrintGraph(self.time)

    def FinalizeSolutionStepPrecompression(self):
        super().FinalizeSolutionStep()
        self.MeasureForcesAndPressurePrecompression()

    def FinalizeSolutionStepDecompression(self):
        super().FinalizeSolutionStep()
        self.MeasureForcesAndPressureDecompression()

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        self.MeasureForcesAndPressure()

    def Finalize(self):
        super().Finalize()

        # TODO: After self.CleanUpOperations() in base class!!
        self.FinalizeGraphs()
        self.PrintMachineLearningData()
    
    def PrintMachineLearningData(self):
        import numpy as np
        from scipy import interpolate

        filename1 = self.absolute_path_to_file1
        filename2 = os.path.join(os.getcwd(), 'triaxial_experiment_' + str(simulation_number) + '.grf')
        filename3 = os.path.join(os.getcwd(), 'triaxial_DEM_' + str(simulation_number) + '.grf')
        X, Y = [], []
        f1 = open(filename1, 'r')
        for line in f1:
            values = [float(s) for s in line.split()]
            X.append(values[0])
            q = 145.0 * values[1]
            Y.append(q)
        f1.close()
        minimum = min(Y)
        min_index = Y.index(minimum)
        X1 = X[:min_index]
        X2 = X[min_index:]
        Y1 = Y[:min_index]
        Y2 = Y[min_index:]
        Y2 = [x - 0.000145 * self.ConfinementPressure * 1e6 for x in Y2]
        index_zero = min(range(len(Y2)), key=lambda i: abs(Y2[i] - 0.0))
        close_zero = X2[index_zero]
        X2 = [x - close_zero for x in X2]
        X2 = X2[index_zero:]
        Y2 = Y2[index_zero:]
        f3 = open(filename3, 'w')
        for i in range(len(X2)):
            f3.write(str(X2[i]) + "  " + str(Y2[i]) + "\n")
        f3.close()
        X3, Y3 = [], []
        f2 = open(filename2, 'r')
        for line in f2:
            values = [float(s) for s in line.split()]
            X3.append(values[0])
            q = values[1]
            Y3.append(q)
        f2.close()
        X_0 = max(min(X2), min(X3))
        X_N = min(max(X2), max(X3))
        X4 = np.linspace(X_0, X_N, num_of_discretization_points)
        F2 = interpolate.interp1d(X2, Y2, kind = 'linear')
        F3 = interpolate.interp1d(X3, Y3, kind = 'linear')
        minimumX = 1.2 * min(min(X2), min(X3))
        maximumX = 1.2 * max(max(X2), max(X3))
        minimumY = 1.2 * min(min(Y2), min(Y3))
        maximumY = 1.2 * max(max(Y2), max(Y3))
        
        if print_images:
            import matplotlib.pyplot as plt
            plt.axis([minimumX, maximumX, minimumY, maximumY])
            plt.xlabel("vertical strain (%)")
            plt.ylabel("q' (psi)")
            plt.plot(X2, Y2, color='blue', linewidth=1, linestyle='solid', marker='None', label="DEM")
            plt.plot(X4, F2(X4), color='blue', markersize=2, linewidth=1, linestyle='dashed', marker='o', label="DEM interpolation")
            plt.plot(X3, Y3, color='red',  linewidth=1, linestyle='solid', marker='None', label="experiment")
            plt.plot(X4, F3(X4), color='red',  markersize=2, linewidth=1, linestyle='dashed', marker='o', label="experiment interpolation")
            plt.title('Triaxial deviatoric effective stress. Precompressed results vs experiments', fontdict = {'fontsize':8})
            plt.legend(loc='lower right')
            plt.grid()
            printfilename = 'dem_triaxial_vs_experiments_' + str(simulation_number) + '.png'
            plt.savefig(printfilename, dpi=300)
            plt.close()

        error = []
        length = len(F2(X4))
        for i in range(length):
            error.append(100 * (F2(X4[i]) - F3(X4[i]))/F3(X4[i]))
        if simulation_number == simulation_number_list[0]:
            machine_learning_data.write(str("%.8g"%self.loose_young).rjust(13) + '\n')
            machine_learning_data.write(str("%.8g"%self.bonded_young).rjust(13) + '\n')
            machine_learning_data.write(str("%.8g"%self.friction).rjust(13) + '\n')
            machine_learning_data.write(str("%.8g"%self.contact_tau).rjust(13) + '\n')
        for i in error:
            machine_learning_data.write(str("%.8g"%i).rjust(13) + '\n')
        machine_learning_data.flush()
        if simulation_number == simulation_number_list[-1]:
            machine_learning_data.close()

    def InitializeMaterialTest(self):

        self.top_mesh_nodes = []; self.bot_mesh_nodes = []; self.top_mesh_fem_nodes = []; self.bot_mesh_fem_nodes = []
        self.xtop_area = 0.0
        self.xbot_area = 0.0
        self.xlat_area = 0.0
        self.xtopcorner_area = 0.0
        self.xbotcorner_area = 0.0
        self.SKIN = list()
        self.LAT = list()
        self.BOT = list()
        self.TOP = list()
        self.XLAT = list()  # only lat, not the corner ones
        self.XTOP = list()  # only top, not corner ones...
        self.XBOT = list()
        self.XTOPCORNER = list()
        self.XBOTCORNER = list()
        self.bond_00_05 = list(); self.bond_05_10 = list(); self.bond_10_15 = list(); self.bond_15_20 = list(); self.bond_20_25 = list(); self.bond_25_30 = list(); self.bond_30_35 = list()
        self.bond_35_40 = list(); self.bond_40_45 = list(); self.bond_45_50 = list(); self.bond_50_55 = list(); self.bond_55_60 = list(); self.bond_60_65 = list(); self.bond_65_70 = list()
        self.bond_70_75 = list(); self.bond_75_80 = list(); self.bond_80_85 = list(); self.bond_85_90 = list()
        self.sizes = []
        self.sigma_mean_table = []; self.tau_mean_table = []; self.sigma_rel_std_dev_table = []; self.tau_rel_std_dev_table = []; self.sigma_ratio_table = [];

        for i in range(0,18):
            self.sizes.append(0.0)
            self.sigma_mean_table.append(0.0)
            self.tau_mean_table.append(0.0)
            self.sigma_rel_std_dev_table.append(0.0)
            self.tau_rel_std_dev_table.append(0.0)
            self.sigma_ratio_table.append(0.0)

        self.graph_counter = 0; self.renew_pressure = 0; self.Pressure = 0.0; self.pressure_to_apply = 0.0; self.CN_graph_counter = 0
        self.length_correction_factor = 1.0
        self.graph_frequency        = int(self.parameters["GraphExportFreq"].GetDouble()/self.spheres_model_part.ProcessInfo.GetValue(DELTA_TIME))
        self.strain = 0.0; self.strain_bts = 0.0; self.volumetric_strain = 0.0; self.radial_strain = 0.0; self.first_time_entry = 1; self.first_time_entry_2 = 1
        self.total_stress_top = 0.0; self.total_stress_bot = 0.0; self.total_stress_mean = 0.0
        self.LoadingVelocity = 0.0
        self.MeasuringSurface = 0.0

        if "material_test_settings" in self.parameters.keys():
            self.height = self.parameters["material_test_settings"]["SpecimenLength"].GetDouble()
            self.diameter = self.parameters["material_test_settings"]["SpecimenDiameter"].GetDouble()
            self.ConfinementPressure = self.parameters["material_test_settings"]["ConfinementPressure"].GetDouble()
            self.y_coordinate_of_cylinder_bottom_base = self.parameters["material_test_settings"]["YCoordinateOfCylinderBottomBase"].GetDouble()
            self.z_coordinate_of_cylinder_bottom_base = self.parameters["material_test_settings"]["ZCoordinateOfCylinderBottomBase"].GetDouble()
        else:
            self.height = self.parameters["SpecimenLength"].GetDouble()
            self.diameter = self.parameters["SpecimenDiameter"].GetDouble()
            self.ConfinementPressure = self.parameters["ConfinementPressure"].GetDouble()
            self.y_coordinate_of_cylinder_bottom_base = self.parameters["YCoordinateOfCylinderBottomBase"].GetDouble()
            self.z_coordinate_of_cylinder_bottom_base = self.parameters["ZCoordinateOfCylinderBottomBase"].GetDouble()
        
        self.ConfinementPressure = ConfinementPressure

        self.ComputeLoadingVelocity()
        self.ComputeMeasuringSurface()
        self.problem_name = self.parameters["problem_name"].GetString()
        self.initial_time = datetime.datetime.now()
        absolute_path_to_file = os.path.join(self.graphs_path, self.problem_name + "_Parameter_chart.grf")
        self.chart = open(absolute_path_to_file, 'w')
        self.aux = AuxiliaryUtilities()
        self.PreUtilities = PreUtilities()
        self.PrepareTests()
        self.PrepareTestTriaxialHydro()
        domain_volume = math.pi * 0.5 * 0.5 * self.diameter * self.diameter * self.height
        DEM_procedures.GranulometryUtils(domain_volume, self.spheres_model_part)

    def MeasureForcesAndPressurePrecompression(self):

        dt = self.spheres_model_part.ProcessInfo.GetValue(DELTA_TIME)
        self.strain += -100 * self.length_correction_factor * self.LoadingVelocity * dt / self.height

        total_force_top = 0.0
        for node in self.top_mesh_nodes:
            force_node_y = node.GetSolutionStepValue(ELASTIC_FORCES)[1]
            total_force_top += force_node_y
        self.total_stress_top = total_force_top / self.MeasuringSurface

        total_force_bot = 0.0
        for node in self.bot_mesh_nodes:
            force_node_y = -node.GetSolutionStepValue(ELASTIC_FORCES)[1]
            total_force_bot += force_node_y
        self.total_stress_bot = total_force_bot / self.MeasuringSurface

        self.total_stress_mean = 0.5 * (self.total_stress_bot + self.total_stress_top)

        if self.SigmaHorizontal:
            self.Pressure = min(self.total_stress_mean, self.SigmaHorizontal)
            self.ApplyLateralPressure(self.Pressure, self.XLAT, self.XBOT, self.XTOP, self.XBOTCORNER, self.XTOPCORNER,self.alpha_top,self.alpha_bot,self.alpha_lat)

        if self.total_stress_mean > self.SigmaVertical:
            self.compression_stage_completed = True

    def MeasureForcesAndPressureDecompression(self):

        dt = self.spheres_model_part.ProcessInfo.GetValue(DELTA_TIME)
        self.strain += -100 * self.length_correction_factor * self.PlatesDecompressionVelocity * dt / self.height

        total_force_top = 0.0
        for node in self.top_mesh_nodes:
            force_node_y = node.GetSolutionStepValue(ELASTIC_FORCES)[1]
            total_force_top += force_node_y
        self.total_stress_top = total_force_top / self.MeasuringSurface

        total_force_bot = 0.0
        for node in self.bot_mesh_nodes:
            force_node_y = -node.GetSolutionStepValue(ELASTIC_FORCES)[1]
            total_force_bot += force_node_y
        self.total_stress_bot = total_force_bot / self.MeasuringSurface

        self.total_stress_mean = 0.5 * (self.total_stress_bot + self.total_stress_top)

        if self.SigmaHorizontal:
            self.Pressure = min(self.total_stress_mean, self.SigmaHorizontal)
            self.ApplyLateralPressure(self.Pressure, self.XLAT, self.XBOT, self.XTOP, self.XBOTCORNER, self.XTOPCORNER,self.alpha_top,self.alpha_bot,self.alpha_lat)

        if self.total_stress_mean < self.SigmaVerticalAlmostZero:
            self.decompression_stage_completed = True

    def MeasureForcesAndPressure(self):

        dt = self.spheres_model_part.ProcessInfo.GetValue(DELTA_TIME)
        self.strain += -100 * self.length_correction_factor * self.LoadingVelocity * dt / self.height

        total_force_top = 0.0
        for node in self.top_mesh_nodes:
            force_node_y = node.GetSolutionStepValue(ELASTIC_FORCES)[1]
            total_force_top += force_node_y
        self.total_stress_top = total_force_top / self.MeasuringSurface

        total_force_bot = 0.0
        for node in self.bot_mesh_nodes:
            force_node_y = -node.GetSolutionStepValue(ELASTIC_FORCES)[1]
            total_force_bot += force_node_y
        self.total_stress_bot = total_force_bot / self.MeasuringSurface

        self.total_stress_mean = 0.5 * (self.total_stress_bot + self.total_stress_top)

        if self.ConfinementPressure:
            self.Pressure = min(self.total_stress_mean, self.ConfinementPressure * 1e6)
            self.ApplyLateralPressure(self.Pressure, self.XLAT, self.XBOT, self.XTOP, self.XBOTCORNER, self.XTOPCORNER,self.alpha_top,self.alpha_bot,self.alpha_lat)

    def ComputeLoadingVelocity(self):
        top_vel = bot_vel = 0.0
        for smp in self.rigid_face_model_part.SubModelParts:
            if smp[IDENTIFIER] == 'TOP':
                top_vel = smp[LINEAR_VELOCITY_Y]
            if smp[IDENTIFIER] == 'BOTTOM':
                bot_vel = smp[LINEAR_VELOCITY_Y]
        self.LoadingVelocity = top_vel - bot_vel

    def ComputeMeasuringSurface(self):
        self.MeasuringSurface = 0.25 * math.pi * self.diameter * self.diameter

    def PrepareTestTriaxialHydro(self):

        ####### Correction Coefs  TODO 0.25* for cylinder section EXXON
        self.alpha_top = math.pi*self.diameter*self.diameter*0.25/(self.xtop_area + 0.70710678*self.xtopcorner_area)
        self.alpha_bot = math.pi*self.diameter*self.diameter*0.25/(self.xbot_area + 0.70710678*self.xbotcorner_area)
        self.alpha_lat = math.pi*self.diameter*self.height/(self.xlat_area + 0.70710678*self.xtopcorner_area + 0.70710678*self.xbotcorner_area)

    def PrepareTests(self):

        absolute_path_to_file1 = os.path.join(self.graphs_path, self.problem_name + "_graph.grf")
        absolute_path_to_file2 = os.path.join(self.graphs_path, self.problem_name + "_graph_top.grf")
        absolute_path_to_file3 = os.path.join(self.graphs_path, self.problem_name + "_graph_bot.grf")
        absolute_path_to_file4 = os.path.join(self.graphs_path, self.problem_name + "_graph_strain_vs_q_in_psi.grf")
        self.graph_export_1 = open(absolute_path_to_file1, 'w')
        self.graph_export_2 = open(absolute_path_to_file2, 'w')
        self.graph_export_3 = open(absolute_path_to_file3, 'w')
        self.graph_export_4 = open(absolute_path_to_file4, 'w')

        (self.xtop_area,self.xbot_area,self.xlat_area,self.xtopcorner_area,self.xbotcorner_area,y_top_total,weight_top, y_bot_total, weight_bot) = self.CylinderSkinDetermination()

        initial_height_top = y_top_total/weight_top
        initial_height_bot = y_bot_total/weight_bot

        inner_initial_height = initial_height_top - initial_height_bot
        extended_length = self.height + (self.height - inner_initial_height)

        self.length_correction_factor = self.height/extended_length

        absolute_path_to_file = os.path.join(self.graphs_path, self.problem_name + "_CN.grf")
        self.CN_export = open(absolute_path_to_file, 'w')

    def CylinderSkinDetermination(self):

        # Cylinder dimensions
        h = self.height
        d = self.diameter
        y_min = self.y_coordinate_of_cylinder_bottom_base

        eps = 3.0 #2.0
        xlat_area = 0.0
        xbot_area = 0.0
        xtop_area = 0.0
        xbotcorner_area = 0.0
        xtopcorner_area = 0.0
        y_top_total = 0.0
        y_bot_total = 0.0
        weight_top = 0.0
        weight_bot = 0.0

        for element in self.spheres_model_part.Elements:

            element.GetNode(0).SetSolutionStepValue(SKIN_SPHERE, 0)

            node = element.GetNode(0)
            r = node.GetSolutionStepValue(RADIUS)
            x = node.X
            y = node.Y
            z = node.Z

            cross_section = math.pi * r * r

            if (x * x + z * z) >= ((0.5 * d - eps * r) * (0.5 * d - eps * r)):

                element.GetNode(0).SetSolutionStepValue(SKIN_SPHERE, 1)
                self.LAT.append(node)

                if (y > y_min + eps * r) and (y < y_min + (h - eps * r)):

                    self.SKIN.append(element)
                    self.XLAT.append(node)

                    xlat_area = xlat_area + cross_section

            if (y <= y_min + eps * r) or (y >= y_min + (h - eps * r)):

                element.GetNode(0).SetSolutionStepValue(SKIN_SPHERE, 1)
                self.SKIN.append(element)

                if y <= y_min + eps * r:

                    self.BOT.append(node)
                    y_bot_total += y*r
                    weight_bot += r

                elif y >= y_min + (h - eps * r):

                    self.TOP.append(node)

                    y_top_total += y*r
                    weight_top += r

                if (x * x + z * z) >= ((0.5 * d - eps * r) * (0.5 * d - eps * r)):

                    if y > y_min + h / 2:

                        self.XTOPCORNER.append(node)
                        xtopcorner_area = xtopcorner_area + cross_section

                    else:

                        self.XBOTCORNER.append(node)
                        xbotcorner_area = xbotcorner_area + cross_section
                else:

                    if y <= y_min + eps * r:

                        self.XBOT.append(node)
                        xbot_area = xbot_area + cross_section

                    elif y >= y_min + (h - eps * r):

                        self.XTOP.append(node)
                        xtop_area = xtop_area + cross_section
        #checks:
        if len(self.XLAT)==0:
            self.procedures.KratosPrintWarning("ERROR! in Cylinder Skin Determination - NO LATERAL PARTICLES" + "\n")
        else:
            self.procedures.KratosPrintInfo(str(h) + " * " + str(d) + " cylinder skin determination" + "\n")

        return (xtop_area, xbot_area, xlat_area, xtopcorner_area, xbotcorner_area, y_top_total, weight_top, y_bot_total, weight_bot)

    def ApplyLateralPressure(self, Pressure, XLAT, XBOT, XTOP, XBOTCORNER, XTOPCORNER, alpha_top, alpha_bot, alpha_lat):

        for node in XLAT:
            r = node.GetSolutionStepValue(RADIUS)
            x = node.X
            y = node.Y
            z = node.Z

            values = Array3()
            vect = Array3()

            cross_section = math.pi * r * r

            # normal vector to the center:
            vect_moduli = math.sqrt(x * x + z * z)

            if vect_moduli > 0.0:
                vect[0] = -x / vect_moduli
                vect[1] = 0
                vect[2] = -z / vect_moduli

            values[0] = cross_section * alpha_lat * Pressure * vect[0]
            values[1] = 0.0
            values[2] = cross_section * alpha_lat * Pressure * vect[2]

            node.SetSolutionStepValue(EXTERNAL_APPLIED_FORCE, values)

        for node in XTOPCORNER:

            r = node.GetSolutionStepValue(RADIUS)
            x = node.X
            y = node.Y
            z = node.Z

            values = Array3()
            vect = Array3()

            cross_section = math.pi * r * r

            # normal vector to the center:
            vect_moduli = math.sqrt(x * x + z * z)

            if vect_moduli > 0.0:
                vect[0] = -x / vect_moduli
                vect[1] = 0
                vect[2] = -z / vect_moduli

            values[0] = cross_section * alpha_lat * Pressure * vect[0] * 0.70710678
            values[1] = 0.0
            values[2] = cross_section * alpha_lat * Pressure * vect[2] * 0.70710678

            node.SetSolutionStepValue(EXTERNAL_APPLIED_FORCE, values)

        for node in XBOTCORNER:

            r = node.GetSolutionStepValue(RADIUS)
            x = node.X
            y = node.Y
            z = node.Z

            values = Array3()
            vect = Array3()

            cross_section = math.pi * r * r

            # vector normal al centre:
            vect_moduli = math.sqrt(x * x + z * z)

            if vect_moduli > 0.0:
                vect[0] = -x / vect_moduli
                vect[1] = 0
                vect[2] = -z / vect_moduli

            values[0] = cross_section * alpha_lat * Pressure * vect[0] * 0.70710678
            values[1] = 0.0
            values[2] = cross_section * alpha_lat * Pressure * vect[2] * 0.70710678

            node.SetSolutionStepValue(EXTERNAL_APPLIED_FORCE, values)

    def PrepareDataForGraph(self):

        prepare_check = [0,0,0,0]
        self.total_check = 0

        for smp in self.rigid_face_model_part.SubModelParts:
            if smp[IDENTIFIER] == 'TOP':
                self.top_mesh_nodes = smp.Nodes
                prepare_check[0] = 1
            if smp[IDENTIFIER] == 'BOTTOM':
                self.bot_mesh_nodes = smp.Nodes
                prepare_check[1] = 1

        for smp in self.spheres_model_part.SubModelParts:
            if smp[IDENTIFIER] == 'TOP':
                self.top_mesh_nodes = smp.Nodes
                prepare_check[2] = -1

            if smp[IDENTIFIER] == 'BOTTOM':
                self.bot_mesh_nodes = smp.Nodes
                prepare_check[3] = -1

        for it in range(len(prepare_check)):
            self.total_check += prepare_check[it]

        if math.fabs(self.total_check) != 2:
            self.Procedures.KratosPrintWarning(" ERROR in the definition of TOP BOT groups. Both groups are required to be defined, they have to be either on FEM groups or in DEM groups")

    def PrintGraph(self, time):

        if self.graph_counter == self.graph_frequency:
            self.graph_counter = 0
            self.graph_export_1.write(str("%.6g"%self.strain).rjust(13) + "  " + str("%.6g"%(self.total_stress_mean * 1e-6)).rjust(13) + "  " + str("%.8g"%time).rjust(12) + '\n')
            self.graph_export_2.write(str("%.8g"%self.strain).rjust(13) + "  " + str("%.6g"%(self.total_stress_top  * 1e-6)).rjust(13) + "  " + str("%.8g"%time).rjust(12) + '\n')
            self.graph_export_3.write(str("%.8g"%self.strain).rjust(13) + "  " + str("%.6g"%(self.total_stress_bot  * 1e-6)).rjust(13) + "  " + str("%.8g"%time).rjust(12) + '\n')
            self.graph_export_1.flush()
            self.graph_export_2.flush()
            self.graph_export_3.flush()
            self.graph_export_4.write(str("%.8g"%self.strain).rjust(15) + "  " + str("%.6g"%(self.total_stress_mean * 145 * 1e-6 - 145 * self.ConfinementPressure)).rjust(13) +  "  " + str("%.8g"%time).rjust(12) + '\n')
            self.graph_export_4.flush()
        self.graph_counter += 1

    def FinalizeGraphs(self):
        # Create a copy and renaming
        self.absolute_path_to_file1 = os.path.join(self.graphs_path, self.problem_name + "_graph.grf")
        absolute_path_to_file2 = os.path.join(self.graphs_path, self.problem_name + "_bts.grf")
        absolute_path_to_file3 = os.path.join(self.graphs_path, self.problem_name + "_graph_VOL.grf")
        for filename in os.listdir("."):
            if filename.startswith(self.absolute_path_to_file1):
                shutil.copy(filename, filename + "COPY")
                os.rename(filename+"COPY", self.absolute_path_to_file1 + str(self.initial_time).replace(":", "") + ".grf")
            if filename.startswith(absolute_path_to_file2):
                shutil.copy(filename, filename + "COPY")
                os.rename(filename+"COPY", absolute_path_to_file2 + str(self.initial_time).replace(":", "") + ".grf")
            if filename.startswith(absolute_path_to_file3):
                shutil.copy(filename, filename + "COPY")
                os.rename(filename+"COPY", absolute_path_to_file3 + str(self.initial_time).replace(":", "") + ".grf")
        self.graph_export_1.close()
        self.graph_export_2.close()
        self.graph_export_3.close()
        self.graph_export_4.close()

if __name__ == "__main__":

    with open("ProjectParametersDEM.json", 'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    machine_learning_file = "machine_learning_data.grf"
    machine_learning_data = open(machine_learning_file, 'w')

    confinement_pressure_list = [0.05, 0.1] # For the triaxials, in MPa
    num_of_discretization_points = 10 # To compute relative errors
    print_images = False
    
    number_of_triaxials = len(confinement_pressure_list)
    simulation_number_list = [i for i in range(1, number_of_triaxials + 1)]

    for simulation_number, ConfinementPressure in zip(simulation_number_list, confinement_pressure_list):
        model = KratosMultiphysics.Model()
        DecompressedMaterialTriaxialTest(model, parameters).Run()
