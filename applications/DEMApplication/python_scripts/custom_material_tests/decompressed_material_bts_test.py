import KratosMultiphysics
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.DEMApplication.DEM_analysis_stage import DEMAnalysisStage
from KratosMultiphysics.DEMApplication import DEM_procedures as DEM_procedures
import math
import datetime

class DecompressedMaterialBTSTest(DEMAnalysisStage):

    def __init__(self, model, parameters):
        super().__init__(model, parameters)

        self.parameters = parameters
        self.compression_stage_completed = False
        self.decompression_stage_completed = False
        self.time_to_print_bts_graph = False

        # Units in Pa and m/s respectively
        self.SigmaHorizontal = self.parameters["material_test_settings"]["SigmaHorizontal"].GetDouble()
        self.SigmaVertical = self.parameters["material_test_settings"]["SigmaVertical"].GetDouble()
        self.SigmaVerticalAlmostZero = self.parameters["material_test_settings"]["SigmaVerticalAlmostZero"].GetDouble()
        self.PlatesDecompressionVelocity = self.parameters["material_test_settings"]["PlatesDecompressionVelocity"].GetDouble()

    def Initialize(self):
        super().Initialize()
        self.InitializeMaterialTest()
        self.PrepareDataForGraph()
        self.ApplyPreCompression()
        self.ApplyDeCompression()
        self.PrepareBTSTest()
        self.ApplyLoadingVelocityToBTSPlates()

    def ApplyLoadingVelocityToBTSPlates(self):
        for smp in self.rigid_face_model_part.SubModelParts:
            if smp[IDENTIFIER] == 'TOP_BTS':
                smp[LINEAR_VELOCITY_Y] = 0.5 * self.LoadingVelocity
            if smp[IDENTIFIER] == 'BOTTOM_BTS':
                smp[LINEAR_VELOCITY_Y] = -0.5 * self.LoadingVelocity

    def ApplyPreCompression(self):

        self._GetSolver().cplusplus_strategy.BreakAllBonds()

        print("\n************************************ Applying PreCompression...\n", flush=True)
        while not self.compression_stage_completed:
            self.time = self._GetSolver().AdvanceInTime(self.time)
            self.InitializeSolutionStep()
            self._GetSolver().Predict()
            self._GetSolver().SolveSolutionStep()
            self.FinalizeSolutionStepPreCompression()
            self.OutputSolutionStep()
        print("\n*************************** Finished Applying PreCompression!!!\n", flush=True)

        self._GetSolver().cplusplus_strategy.HealAllBonds()
        ParallelBondUtilities().SetCurrentIndentationAsAReferenceInParallelBonds(self.spheres_model_part)
        PreUtilities().ResetSkinParticles(self.spheres_model_part)
        self._GetSolver().cplusplus_strategy.ComputeSkin(self.spheres_model_part, 1.5)

    def ResetLoadingVelocity(self):
        for smp in self.rigid_face_model_part.SubModelParts:
            if smp[IDENTIFIER] == 'TOP':
                smp[LINEAR_VELOCITY_Z] = 0.5 * self.PlatesDecompressionVelocity
            if smp[IDENTIFIER] == 'BOTTOM':
                smp[LINEAR_VELOCITY_Z] = -0.5 * self.PlatesDecompressionVelocity

    def ApplyDeCompression(self):

        print("\n************************************ Applying DeCompression...\n", flush=True)
        while not self.decompression_stage_completed:
            self.time = self._GetSolver().AdvanceInTime(self.time)
            self.InitializeSolutionStep()
            self._GetSolver().Predict()
            self.ResetLoadingVelocity()
            self._GetSolver().SolveSolutionStep()
            self.FinalizeSolutionStepDeCompression()
            self.OutputSolutionStep()
        print("\n*************************** Finished Applying DeCompression!!!\n", flush=True)

    def RunSolutionLoop(self):

        print("\n************************************ Applying standard BTS...\n", flush=True)
        while self.KeepAdvancingSolutionLoop():
            self.time = self._GetSolver().AdvanceInTime(self.time)
            self.InitializeSolutionStep()
            self._GetSolver().Predict()
            self._GetSolver().SolveSolutionStep()
            self.FinalizeSolutionStep()
            self.OutputSolutionStep()
        print("\n*************************** Finished Applying standard BTS...\n", flush=True)

    def OutputSolutionStep(self):
        super().OutputSolutionStep()
        self.PrintGraph(self.time)

    def FinalizeSolutionStepPreCompression(self):
        super().FinalizeSolutionStep()
        self.MeasureForcesAndPressurePreCompression()

    def FinalizeSolutionStepDeCompression(self):
        super().FinalizeSolutionStep()
        self.MeasureForcesAndPressureDeCompression()

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        self.MeasureForcesAndPressure()

    def Finalize(self):
        super().Finalize()

        # TODO: After self.CleanUpOperations() in base class!!
        self.FinalizeGraphs()

    def InitializeMaterialTest(self):

        self.top_mesh_nodes = []; self.bot_mesh_nodes = []; self.top_mesh_nodes_bts = []; self.bot_mesh_nodes_bts = []
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

    def MeasureForcesAndPressurePreCompression(self):

        dt = self.spheres_model_part.ProcessInfo.GetValue(DELTA_TIME)
        self.strain += -100 * self.length_correction_factor * self.LoadingVelocity * dt / self.height

        total_force_top = 0.0
        for node in self.top_mesh_nodes:
            force_node_z = node.GetSolutionStepValue(ELASTIC_FORCES)[2]
            total_force_top += force_node_z
        self.total_stress_top = total_force_top / self.MeasuringSurface

        total_force_bot = 0.0
        for node in self.bot_mesh_nodes:
            force_node_z = -node.GetSolutionStepValue(ELASTIC_FORCES)[2]
            total_force_bot += force_node_z
        self.total_stress_bot = total_force_bot / self.MeasuringSurface

        self.total_stress_mean = 0.5 * (self.total_stress_bot + self.total_stress_top)

        if self.SigmaHorizontal:
            self.Pressure = min(self.total_stress_mean, self.SigmaHorizontal)
            self.ApplyLateralPressure(self.Pressure, self.XLAT, self.XBOT, self.XTOP, self.XBOTCORNER, self.XTOPCORNER,self.alpha_top,self.alpha_bot,self.alpha_lat)

        if self.total_stress_mean > self.SigmaVertical:
            self.compression_stage_completed = True

    def MeasureForcesAndPressureDeCompression(self):

        dt = self.spheres_model_part.ProcessInfo.GetValue(DELTA_TIME)
        self.strain += -100 * self.length_correction_factor * self.PlatesDecompressionVelocity * dt / self.height

        total_force_top = 0.0
        for node in self.top_mesh_nodes:
            force_node_z = node.GetSolutionStepValue(ELASTIC_FORCES)[2]
            total_force_top += force_node_z
        self.total_stress_top = total_force_top / self.MeasuringSurface

        total_force_bot = 0.0
        for node in self.bot_mesh_nodes:
            force_node_z = -node.GetSolutionStepValue(ELASTIC_FORCES)[2]
            total_force_bot += force_node_z
        self.total_stress_bot = total_force_bot / self.MeasuringSurface

        self.total_stress_mean = 0.5 * (self.total_stress_bot + self.total_stress_top)

        if self.SigmaHorizontal:
            self.Pressure = min(self.total_stress_mean, self.SigmaHorizontal)
            self.ApplyLateralPressure(self.Pressure, self.XLAT, self.XBOT, self.XTOP, self.XBOTCORNER, self.XTOPCORNER,self.alpha_top,self.alpha_bot,self.alpha_lat)

        if self.total_stress_mean < self.SigmaVerticalAlmostZero:
            self.decompression_stage_completed = True

    def MeasureForcesAndPressure(self):

        dt = self.spheres_model_part.ProcessInfo.GetValue(DELTA_TIME)

        total_force_top = 0.0
        for node in self.top_mesh_nodes_bts:
            force_node_y = node.GetSolutionStepValue(ELASTIC_FORCES)[1]
            total_force_top += force_node_y

        total_force_bot = 0.0
        for node in self.bot_mesh_nodes_bts:
            force_node_y = -node.GetSolutionStepValue(ELASTIC_FORCES)[1]
            total_force_bot += force_node_y

        total_force_bts = 0.5 * (total_force_bot + total_force_top)

        self.total_stress_bts = 2.0 * total_force_bts / (math.pi * self.height * self.diameter)
        self.strain_bts += -100 * self.LoadingVelocity * dt / self.diameter

    def ComputeLoadingVelocity(self):
        top_vel = bot_vel = 0.0
        for smp in self.rigid_face_model_part.SubModelParts:
            if smp[IDENTIFIER] == 'TOP':
                top_vel = smp[LINEAR_VELOCITY_Z]
            if smp[IDENTIFIER] == 'BOTTOM':
                bot_vel = smp[LINEAR_VELOCITY_Z]
        self.LoadingVelocity = top_vel - bot_vel

    def ComputeMeasuringSurface(self):
        self.MeasuringSurface = 0.25 * math.pi * self.diameter * self.diameter

    def PrepareTestTriaxialHydro(self):

        ####### Correction Coefs  TODO 0.25* for cylinder section EXXON
        self.alpha_top = math.pi*self.diameter*self.diameter*0.25/(self.xtop_area + 0.70710678*self.xtopcorner_area)
        self.alpha_bot = math.pi*self.diameter*self.diameter*0.25/(self.xbot_area + 0.70710678*self.xbotcorner_area)
        self.alpha_lat = math.pi*self.diameter*self.height/(self.xlat_area + 0.70710678*self.xtopcorner_area + 0.70710678*self.xbotcorner_area)

    def PrepareBTSTest(self):
        absolute_path_to_file5 = os.path.join(self.graphs_path, self.problem_name + "_graph_bts.grf")
        self.graph_export_5 = open(absolute_path_to_file5, 'w')
        self.time_to_print_bts_graph = True

    def PrepareTests(self):

        absolute_path_to_file1 = os.path.join(self.graphs_path, self.problem_name + "_graph.grf")
        absolute_path_to_file2 = os.path.join(self.graphs_path, self.problem_name + "_graph_top.grf")
        absolute_path_to_file3 = os.path.join(self.graphs_path, self.problem_name + "_graph_bot.grf")
        absolute_path_to_file4 = os.path.join(self.graphs_path, self.problem_name + "_graph_strain_vs_q.grf")
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
        z_min = self.z_coordinate_of_cylinder_bottom_base

        eps = 3.0 #2.0
        xlat_area = 0.0
        xbot_area = 0.0
        xtop_area = 0.0
        xbotcorner_area = 0.0
        xtopcorner_area = 0.0
        z_top_total = 0.0
        z_bot_total = 0.0
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

            if (x * x + y * y) >= ((0.5 * d - eps * r) * (0.5 * d - eps * r)):

                element.GetNode(0).SetSolutionStepValue(SKIN_SPHERE, 1)
                self.LAT.append(node)

                if (z > z_min + eps * r) and (z < z_min + (h - eps * r)):

                    self.SKIN.append(element)
                    self.XLAT.append(node)

                    xlat_area = xlat_area + cross_section

            if (z <= z_min + eps * r) or (z >= z_min + (h - eps * r)):

                element.GetNode(0).SetSolutionStepValue(SKIN_SPHERE, 1)
                self.SKIN.append(element)

                if z <= z_min + eps * r:

                    self.BOT.append(node)
                    z_bot_total += z*r
                    weight_bot += r

                elif z >= z_min + (h - eps * r):

                    self.TOP.append(node)

                    z_top_total += z*r
                    weight_top += r

                if (x * x + y * y) >= ((0.5 * d - eps * r) * (0.5 * d - eps * r)):

                    if z > z_min + h / 2:

                        self.XTOPCORNER.append(node)
                        xtopcorner_area = xtopcorner_area + cross_section

                    else:

                        self.XBOTCORNER.append(node)
                        xbotcorner_area = xbotcorner_area + cross_section
                else:

                    if z <= z_min + eps * r:

                        self.XBOT.append(node)
                        xbot_area = xbot_area + cross_section

                    elif z >= z_min + (h - eps * r):

                        self.XTOP.append(node)
                        xtop_area = xtop_area + cross_section
        #checks:
        if len(self.XLAT)==0:
            self.procedures.KratosPrintWarning("ERROR! in Cylinder Skin Determination - NO LATERAL PARTICLES" + "\n")
        else:
            self.procedures.KratosPrintInfo(str(h) + " * " + str(d) + " cylinder skin determination" + "\n")

        return (xtop_area, xbot_area, xlat_area, xtopcorner_area, xbotcorner_area, z_top_total, weight_top, z_bot_total, weight_bot)

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
            vect_moduli = math.sqrt(x * x + y * y)

            if vect_moduli > 0.0:
                vect[0] = -x / vect_moduli
                vect[1] = -y / vect_moduli
                vect[2] = 0.0

            values[0] = cross_section * alpha_lat * Pressure * vect[0]
            values[1] = cross_section * alpha_lat * Pressure * vect[1]
            values[2] = 0.0

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
            vect_moduli = math.sqrt(x * x + y * y)

            if vect_moduli > 0.0:
                vect[0] = -x / vect_moduli
                vect[1] = -y / vect_moduli
                vect[2] = 0.0

            values[0] = cross_section * alpha_lat * Pressure * vect[0] * 0.70710678
            values[1] = cross_section * alpha_lat * Pressure * vect[1] * 0.70710678
            values[2] = 0.0

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
                vect[1] = -y / vect_moduli
                vect[2] = 0.0

            values[0] = cross_section * alpha_lat * Pressure * vect[0] * 0.70710678
            values[1] = cross_section * alpha_lat * Pressure * vect[1] * 0.70710678
            values[2] = 0.0

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

        for smp in self.rigid_face_model_part.SubModelParts:
            if smp[IDENTIFIER] == 'TOP_BTS':
                self.top_mesh_nodes_bts = smp.Nodes
            if smp[IDENTIFIER] == 'BOTTOM_BTS':
                self.bot_mesh_nodes_bts = smp.Nodes

        for it in range(len(prepare_check)):
            self.total_check += prepare_check[it]

        if math.fabs(self.total_check) != 2:
            self.Procedures.KratosPrintWarning(" ERROR in the definition of TOP BOT groups. Both groups are required to be defined, they have to be either on FEM groups or in DEM groups")

    def PrintGraph(self, time):

        if self.graph_counter == self.graph_frequency:
            self.graph_counter = 0
            total_stress_q = self.total_stress_mean * 1e-6 - self.ConfinementPressure
            self.graph_export_1.write(str("%.6g"%self.strain).rjust(13)     + "  " + str("%.6g"%(self.total_stress_mean * 1e-6)).rjust(13) + "  " + str("%.8g"%time).rjust(12) + '\n')
            self.graph_export_2.write(str("%.8g"%self.strain).rjust(13)     + "  " + str("%.6g"%(self.total_stress_top  * 1e-6)).rjust(13) + "  " + str("%.8g"%time).rjust(12) + '\n')
            self.graph_export_3.write(str("%.8g"%self.strain).rjust(13)     + "  " + str("%.6g"%(self.total_stress_bot  * 1e-6)).rjust(13) + "  " + str("%.8g"%time).rjust(12) + '\n')
            self.graph_export_4.write(str("%.8g"%self.strain).rjust(13)     + "  " + str("%.6g"%(total_stress_q)).rjust(13)                + "  " + str("%.8g"%time).rjust(12) + '\n')
            self.graph_export_1.flush()
            self.graph_export_2.flush()
            self.graph_export_3.flush()
            self.graph_export_4.flush()
            if self.time_to_print_bts_graph:
                self.graph_export_5.write(str("%.8g"%self.strain_bts).rjust(13) + "  " + str("%.6g"%(self.total_stress_bts  * 1e-6)).rjust(13) + "  " + str("%.8g"%time).rjust(12) + '\n')
                self.graph_export_5.flush()
        self.graph_counter += 1

    def FinalizeGraphs(self):
        # Create a copy and renaming
        absolute_path_to_file1 = os.path.join(self.graphs_path, self.problem_name + "_graph.grf")
        absolute_path_to_file2 = os.path.join(self.graphs_path, self.problem_name + "_bts.grf")
        absolute_path_to_file3 = os.path.join(self.graphs_path, self.problem_name + "_graph_VOL.grf")
        for filename in os.listdir("."):
            if filename.startswith(absolute_path_to_file1):
                shutil.copy(filename, filename + "COPY")
                os.rename(filename+"COPY", absolute_path_to_file1 + str(self.initial_time).replace(":", "") + ".grf")
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
        self.graph_export_5.close()

if __name__ == "__main__":

    with open("ProjectParametersDEM.json", 'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    DecompressedMaterialBTSTest(model, parameters).Run()
