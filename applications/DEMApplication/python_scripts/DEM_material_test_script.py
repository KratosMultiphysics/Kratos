import math
import datetime
import shutil
import weakref
import os

from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.DEMApplication import DEM_procedures as DEM_procedures

class MaterialTest():

    def __init__(self, DEM_parameters, procedures, solver, graphs_path, post_path, spheres_model_part, rigid_face_model_part):

        self.parameters = DEM_parameters
        self.graphs_path = graphs_path
        self.post_path = post_path
        self.spheres_model_part = spheres_model_part
        self.rigid_face_model_part = rigid_face_model_part
        self.Procedures = weakref.proxy(procedures)
        self.solver = weakref.proxy(solver)

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

        self.graph_frequency        = int(self.parameters["GraphExportFreq"].GetDouble()/spheres_model_part.ProcessInfo.GetValue(DELTA_TIME))
        self.strain = 0.0; self.strain_bts = 0.0; self.volumetric_strain = 0.0; self.radial_strain = 0.0; self.first_time_entry = 1; self.first_time_entry_2 = 1
        self.total_stress_top = 0.0; self.total_stress_bot = 0.0; self.total_stress_mean = 0.0

        self.new_strain = 0.0
        self.LoadingVelocity = 0.0
        self.MeasuringSurface = 1.0

        # for the graph plotting
        if "material_test_settings" in DEM_parameters.keys():
            self.height = self.parameters["material_test_settings"]["SpecimenLength"].GetDouble()
            self.diameter = self.parameters["material_test_settings"]["SpecimenDiameter"].GetDouble()
            self.ConfinementPressure = self.parameters["material_test_settings"]["ConfinementPressure"].GetDouble()
            self.test_type = self.parameters["material_test_settings"]["TestType"].GetString()
            self.y_coordinate_of_cylinder_bottom_base = self.parameters["material_test_settings"]["YCoordinateOfCylinderBottomBase"].GetDouble()
            self.z_coordinate_of_cylinder_bottom_base = self.parameters["material_test_settings"]["ZCoordinateOfCylinderBottomBase"].GetDouble()
        else:
            self.height = self.parameters["SpecimenLength"].GetDouble()
            self.diameter = self.parameters["SpecimenDiameter"].GetDouble()
            self.ConfinementPressure = self.parameters["ConfinementPressure"].GetDouble()
            self.test_type = self.parameters["TestType"].GetString()
            self.y_coordinate_of_cylinder_bottom_base = self.parameters["YCoordinateOfCylinderBottomBase"].GetDouble()
            self.z_coordinate_of_cylinder_bottom_base = self.parameters["ZCoordinateOfCylinderBottomBase"].GetDouble()

        self.ComputeLoadingVelocity()
        self.ComputeMeasuringSurface()
        self.problem_name = self.parameters["problem_name"].GetString()
        self.initial_time = datetime.datetime.now()

        # self.energy_plot = open(energy_plot, 'w')
        absolute_path_to_file = os.path.join(graphs_path, self.problem_name + "_Parameter_chart.grf")
        self.chart = open(absolute_path_to_file, 'w')
        self.aux = AuxiliaryUtilities()
        self.PreUtilities = PreUtilities()

    def Initialize(self):
        self.PrepareTests()
        self.PrepareTestTriaxialHydro()
        self.PrepareTestOedometric()

        domain_volume = math.pi * 0.5 * 0.5 * self.diameter * self.diameter * self.height
        DEM_procedures.GranulometryUtils(domain_volume, self.spheres_model_part)

    def BreakBondUtility(self, spheres_model_part):
        self.PreUtilities.BreakBondUtility(self.spheres_model_part)

    def Flush(self,a):
        a.flush()

    def PrepareTestOedometric(self):

        if self.test_type == "Oedometric":

            for node in self.LAT:

                node.SetSolutionStepValue(VELOCITY_X, 0.0)
                node.SetSolutionStepValue(VELOCITY_Z, 0.0)
                node.Fix(VELOCITY_X)
                node.Fix(VELOCITY_Z)

    def PrepareTestTriaxialHydro(self):

        if self.test_type == "Triaxial" or self.test_type == "Hydrostatic":
            ####### Correction Coefs  TODO 0.25* for cylinder section EXXON
            self.alpha_top = math.pi*self.diameter*self.diameter*0.25/(self.xtop_area + 0.70710678*self.xtopcorner_area)
            self.alpha_bot = math.pi*self.diameter*self.diameter*0.25/(self.xbot_area + 0.70710678*self.xbotcorner_area)
            self.alpha_lat = math.pi*self.diameter*self.height/(self.xlat_area + 0.70710678*self.xtopcorner_area + 0.70710678*self.xbotcorner_area)

    def PrepareTests(self):

        ##Fixing horizontally top and bot
        if self.test_type != "BTS":
            for node in self.TOP:
                node.SetSolutionStepValue(VELOCITY_X, 0.0)
                node.SetSolutionStepValue(VELOCITY_Z, 0.0)
                node.Fix(VELOCITY_X)
                node.Fix(VELOCITY_Z)

            for node in self.BOT:
                node.SetSolutionStepValue(VELOCITY_X, 0.0)
                node.SetSolutionStepValue(VELOCITY_Z, 0.0)
                node.Fix(VELOCITY_X)
                node.Fix(VELOCITY_Z)

        if self.test_type == "BTS":
            absolute_path_to_file = os.path.join(self.graphs_path, self.problem_name + ".grf")
            self.bts_export = open(absolute_path_to_file, 'w')
            self.BtsSkinDetermination()

        elif self.test_type == "Shear":
            self.BreakBondUtility(self.spheres_model_part)
            absolute_path_to_file1 = os.path.join(self.graphs_path, self.problem_name + "_graph.grf")
            absolute_path_to_file2 = os.path.join(self.graphs_path, self.problem_name + "_graph_top.grf")
            absolute_path_to_file3 = os.path.join(self.graphs_path, self.problem_name + "_graph_bot.grf")
            self.graph_export_1 = open(absolute_path_to_file1, 'w')
            self.graph_export_2 = open(absolute_path_to_file2, 'w')
            self.graph_export_3 = open(absolute_path_to_file3, 'w')

        else:
            absolute_path_to_file1 = os.path.join(self.graphs_path, self.problem_name + "_graph.grf")
            absolute_path_to_file2 = os.path.join(self.graphs_path, self.problem_name + "_graph_top.grf")
            absolute_path_to_file3 = os.path.join(self.graphs_path, self.problem_name + "_graph_bot.grf")
            absolute_path_to_file4 = os.path.join(self.graphs_path, self.problem_name + "_graph_strain_vs_q_in_psi.grf")
            self.graph_export_1 = open(absolute_path_to_file1, 'w')
            self.graph_export_2 = open(absolute_path_to_file2, 'w')
            self.graph_export_3 = open(absolute_path_to_file3, 'w')
            self.graph_export_4 = open(absolute_path_to_file4, 'w')

            if self.test_type == "Hydrostatic":
                absolute_path_to_file = os.path.join(self.graphs_path, self.problem_name + "_graph_VOL.grf")
                self.graph_export_volumetric   = open(absolute_path_to_file, 'w')

            self.Procedures.KratosPrintInfo('Initial Height of the Model: ' + str(self.height)+'\n')

            (self.xtop_area,self.xbot_area,self.xlat_area,self.xtopcorner_area,self.xbotcorner_area,y_top_total,weight_top, y_bot_total, weight_bot) = self.CylinderSkinDetermination()

            initial_height_top = y_top_total/weight_top
            initial_height_bot = y_bot_total/weight_bot

            inner_initial_height = initial_height_top - initial_height_bot
            extended_length = self.height + (self.height - inner_initial_height)

            self.length_correction_factor = self.height/extended_length

        absolute_path_to_file = os.path.join(self.graphs_path, self.problem_name + "_CN.grf")
        self.CN_export = open(absolute_path_to_file, 'w')

    def ComputeLoadingVelocity(self):
        top_vel = bot_vel = 0.0
        for smp in self.rigid_face_model_part.SubModelParts:
            if smp[TOP]:
                top_vel = smp[LINEAR_VELOCITY_Y]
            if smp[BOTTOM]:
                bot_vel = smp[LINEAR_VELOCITY_Y]
        self.LoadingVelocity = top_vel - bot_vel

    def ComputeMeasuringSurface(self):
        self.MeasuringSurface = 0.25 * math.pi * self.diameter * self.diameter

    def CylinderSkinDetermination(self): #model_part, solver, DEM_parameters):

        # SKIN DETERMINATION
        total_cross_section = 0.0

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
            self.Procedures.KratosPrintWarning("ERROR! in Cylinder Skin Determination - NO LATERAL PARTICLES" + "\n")
        else:
            self.Procedures.KratosPrintInfo(str(h) + " * " + str(d) + " cylinder skin determination" + "\n")

        return (xtop_area, xbot_area, xlat_area, xtopcorner_area, xbotcorner_area, y_top_total, weight_top, y_bot_total, weight_bot)

    def BtsSkinDetermination(self):

        # BTS SKIN DETERMINATION
        # Cylinder dimensions
        h = self.height
        d = self.diameter
        eps = 3.0 #2.0
        z_min = self.z_coordinate_of_cylinder_bottom_base

        for element in self.spheres_model_part.Elements:

            element.GetNode(0).SetSolutionStepValue(SKIN_SPHERE, 0)
            node = element.GetNode(0)
            r = node.GetSolutionStepValue(RADIUS)
            x = node.X
            y = node.Y
            z = node.Z

            if (x * x + y * y) >= ((0.5 * d - eps * r) * (0.5 * d - eps * r)):
                element.GetNode(0).SetSolutionStepValue(SKIN_SPHERE, 1)

            if (z <= z_min + eps * r) or (z >= z_min + (h - eps * r)):
                element.GetNode(0).SetSolutionStepValue(SKIN_SPHERE, 1)

        self.Procedures.KratosPrintInfo("Finished computing the skin of the BTS specimen..." + "\n")

    def PrepareDataForGraph(self):

        prepare_check = [0,0,0,0]
        self.total_check = 0

        for smp in self.rigid_face_model_part.SubModelParts:
            if smp[TOP]:
                self.top_mesh_nodes = smp.Nodes
                prepare_check[0] = 1
            if smp[BOTTOM]:
                self.bot_mesh_nodes = smp.Nodes
                prepare_check[1] = 1

        for smp in self.spheres_model_part.SubModelParts:
            if smp[TOP]:
                self.top_mesh_nodes = smp.Nodes
                prepare_check[2] = -1

            if smp[BOTTOM]:
                self.bot_mesh_nodes = smp.Nodes
                prepare_check[3] = -1

        for it in range(len(prepare_check)):

            self.total_check += prepare_check[it]

        if math.fabs(self.total_check) != 2:

            self.Procedures.KratosPrintWarning(" ERROR in the definition of TOP BOT groups. Both groups are required to be defined, they have to be either on FEM groups or in DEM groups")

    def MeasureForcesAndPressure(self):

        dt = self.spheres_model_part.ProcessInfo.GetValue(DELTA_TIME)

        self.strain += -100 * self.length_correction_factor * self.LoadingVelocity * dt / self.height

        if self.test_type =="BTS":

            total_force_bts = 0.0

            for node in self.top_mesh_nodes:

                force_node_y = node.GetSolutionStepValue(ELASTIC_FORCES)[1]
                total_force_bts += force_node_y

            self.total_stress_bts = 2.0 * total_force_bts / (math.pi * self.height * self.diameter)
            self.strain_bts += -100 * self.LoadingVelocity * dt / self.diameter

        else:

            if self.test_type =="Hydrostatic":
                radial_strain = -100*self.MeasureRadialStrain()
                self.volumetric_strain = self.strain + 2.0*radial_strain

            total_force_top = 0.0
            total_force_bot = 0.0

            for node in self.top_mesh_nodes:
                force_node_y = node.GetSolutionStepValue(ELASTIC_FORCES)[1]
                total_force_top += force_node_y

            self.total_stress_top = total_force_top / self.MeasuringSurface

            for node in self.bot_mesh_nodes:
                force_node_y = -node.GetSolutionStepValue(ELASTIC_FORCES)[1]
                total_force_bot += force_node_y

            self.total_stress_bot = total_force_bot / self.MeasuringSurface

            self.total_stress_mean = 0.5 * (self.total_stress_bot + self.total_stress_top)

            if self.test_type =="Shear":
                self.strain += dt
                self.total_stress_top = total_force_top/1.0 # applied force divided by efective shear cylinder area 2*pi*0.0225*0.08
                self.total_stress_mean = self.total_stress_top

            if (self.test_type == "Triaxial" or self.test_type == "Hydrostatic") and self.ConfinementPressure:

                self.Pressure = min(self.total_stress_mean, self.ConfinementPressure * 1e6)

                if self.test_type == "Hydrostatic":
                    self.Pressure = self.total_stress_mean

                self.ApplyLateralPressure(self.Pressure, self.XLAT, self.XBOT, self.XTOP, self.XBOTCORNER, self.XTOPCORNER,self.alpha_top,self.alpha_bot,self.alpha_lat)

    def PrintGraph(self, time):

        if self.graph_counter == self.graph_frequency:
            self.graph_counter = 0

            if self.test_type == "BTS":
                self.bts_export.write(str("%.8g"%time).rjust(12) + "  " + str("%.6g"%(self.total_stress_bts * 1e-6)).rjust(13) + '\n')
                self.Flush(self.bts_export)
            else:
                self.graph_export_1.write(str("%.6g"%self.strain).rjust(13) + "  " + str("%.6g"%(self.total_stress_mean * 1e-6)).rjust(13) + "  " + str("%.8g"%time).rjust(12) + '\n')
                self.graph_export_2.write(str("%.8g"%self.strain).rjust(15) + "  " + str("%.6g"%(self.total_stress_top * 1e-6)).rjust(13)+'\n')
                self.graph_export_3.write(str("%.8g"%self.strain).rjust(15) + "  " + str("%.6g"%(self.total_stress_bot * 1e-6)).rjust(13)+'\n')
                self.Flush(self.graph_export_1)
                self.Flush(self.graph_export_2)
                self.Flush(self.graph_export_3)

                if self.test_type != "Shear":
                    self.graph_export_4.write(str("%.8g"%self.strain).rjust(15) + "  " + str("%.6g"%(self.total_stress_mean * 1e-6 - self.ConfinementPressure)).rjust(13) + '\n')
                    self.Flush(self.graph_export_4)

                if self.test_type == "Hydrostatic":
                    self.graph_export_volumetric.write(str("%.8g"%self.volumetric_strain).rjust(12) + "    " + str("%.6g"%(self.total_stress_mean * 1e-6)).rjust(13) + '\n')
                    self.Flush(self.graph_export_volumetric)

        self.graph_counter += 1

    def PrintCoordinationNumberGraph(self, time, solver):

        if self.CN_graph_counter == self.graph_frequency:
            self.CN_graph_counter = 0
            dummy = 0
            CN = self.solver.cplusplus_strategy.ComputeCoordinationNumber(dummy)
            self.CN_export.write(str("%.8g"%time).rjust(12) + "  " + str(CN) + '\n')
            self.Flush(self.CN_export)

        self.CN_graph_counter += 1

    def PrintChart(self):

        loading_velocity = self.LoadingVelocity

        print ('************DEM VIRTUAL LAB******************'+'\n')
        print ('Loading velocity: ' + str(loading_velocity) + '\n')
        print ('Expected maximum deformation: ' + str(-loading_velocity*self.parameters["FinalTime"].GetDouble() /self.height*100) +'%'+'\n'+'\n'  )

        self.chart.write(("***********PARAMETERS*****************")+'\n')
        self.chart.write( "                                    " +'\n')
        self.chart.write( "    DENSI  = " + (str(self.spheres_model_part.GetProperties()[1][PARTICLE_DENSITY]).rjust(3))+" Kg/m3     "+'\n')
        self.chart.write( "    STAFRC = " + (str(self.spheres_model_part.GetProperties()[1][STATIC_FRICTION]).rjust(3))+"           "+'\n')
        self.chart.write( "    DYNFRC = " + (str(self.spheres_model_part.GetProperties()[1][DYNAMIC_FRICTION]).rjust(3))+"          " +'\n')
        self.chart.write( "    FRCDEC = " + (str(self.spheres_model_part.GetProperties()[1][FRICTION_DECAY]).rjust(3))+"          " +'\n')
        self.chart.write( "    YOUNG  = " + (str(self.spheres_model_part.GetProperties()[1][YOUNG_MODULUS]/1e9).rjust(3))+" GPa"+"     " +'\n')
        self.chart.write( "    POISS  = " + (str(self.spheres_model_part.GetProperties()[1][POISSON_RATIO]).rjust(3))+"           " +'\n')
        self.chart.write( "    FTS    = " + (str(self.spheres_model_part.GetProperties()[1][CONTACT_SIGMA_MIN]).rjust(3))+" Pa        " +'\n')
        self.chart.write( "    LCS1   = " + (str(self.spheres_model_part.GetProperties()[1][SLOPE_LIMIT_COEFF_C1]).rjust(3))+" Pa       " +'\n')
        self.chart.write( "    LCS2   = " + (str(self.spheres_model_part.GetProperties()[1][SLOPE_LIMIT_COEFF_C2]).rjust(3))+" Pa       " +'\n')
        self.chart.write( "    LCS3   = " + (str(self.spheres_model_part.GetProperties()[1][SLOPE_LIMIT_COEFF_C3]).rjust(3))+" Pa       " +'\n')
        self.chart.write( "    YRC1   = " + (str(self.spheres_model_part.GetProperties()[1][SLOPE_FRACTION_N1]).rjust(3))+"           " +'\n')
        self.chart.write( "    YRC2   = " + (str(self.spheres_model_part.GetProperties()[1][SLOPE_FRACTION_N2]).rjust(3))+"           " +'\n')
        self.chart.write( "    YRC3   = " + (str(self.spheres_model_part.GetProperties()[1][SLOPE_FRACTION_N3]).rjust(3))+"           " +'\n')
        self.chart.write( "    FSS    = " + (str(self.spheres_model_part.GetProperties()[1][CONTACT_TAU_ZERO]).rjust(3))+" Pa       " +'\n')
        self.chart.write( "    YEP    = " + (str(self.spheres_model_part.GetProperties()[1][YOUNG_MODULUS_PLASTIC]/1e9).rjust(3))+" GPa"+"     " +'\n')
        self.chart.write( "    YIELD  = " + (str(self.spheres_model_part.GetProperties()[1][PLASTIC_YIELD_STRESS]).rjust(3))+" Pa       " +'\n')
        self.chart.write( "    EDR    = " + (str(self.spheres_model_part.GetProperties()[1][DAMAGE_FACTOR]).rjust(3))+"           " +'\n')
        self.chart.write( "    SEC    = " + (str(self.spheres_model_part.GetProperties()[1][SHEAR_ENERGY_COEF]).rjust(3))+"           " +'\n')
        self.chart.write( "                                    " +'\n')
        self.chart.write( "**************************************" +'\n')
        self.chart.close()

        absolute_path_to_file = os.path.join(self.graphs_path, self.problem_name + "_Parameter_chart.grf")
        data_extract_for_print = open(absolute_path_to_file,"r")

        for line in data_extract_for_print.readlines():
            self.Procedures.KratosPrintInfo(line)
        data_extract_for_print.close()

    def FinalizeGraphs(self):

        #Create a copy and renaming
        absolute_path_to_file1 = os.path.join(self.graphs_path, self.problem_name + "_graph.grf")
        absolute_path_to_file2 = os.path.join(self.graphs_path, self.problem_name + "_bts.grf")
        absolute_path_to_file3 = os.path.join(self.graphs_path, self.problem_name + "_graph_VOL.grf")
        for filename in os.listdir("."):
            if filename.startswith(absolute_path_to_file1):
                shutil.copy(filename, filename+"COPY")
                os.rename(filename+"COPY", absolute_path_to_file1 + str(self.initial_time).replace(":", "") + ".grf")
            if filename.startswith(absolute_path_to_file2):
                shutil.copy(filename, filename+"COPY")
                os.rename(filename+"COPY", absolute_path_to_file2 + str(self.initial_time).replace(":", "") + ".grf")
            if filename.startswith(absolute_path_to_file3):
                shutil.copy(filename, filename+"COPY")
                os.rename(filename+"COPY", absolute_path_to_file3 + str(self.initial_time).replace(":", "") + ".grf")

        if self.test_type == "BTS":
            self.bts_export.close()
        else:
            self.graph_export_1.close()
            self.graph_export_2.close()
            self.graph_export_3.close()

            if self.test_type != "Shear":
                self.graph_export_4.close()

            if self.test_type == "Hydrostatic":
                self.graph_export_volumetric.close()

    def OrientationStudy(self,contact_model_part,step):

        absolute_path_to_file = os.path.join(self.graphs_path, "OrientationChart_"+str(step))
        OrientationChart = open(absolute_path_to_file, 'w')
        counter = 1

        for element in contact_model_part.Elements:
            u1 = element.GetNode(1).X - element.GetNode(0).X
            u2 = element.GetNode(1).Y - element.GetNode(0).Y
            u3 = element.GetNode(1).Z - element.GetNode(0).Z

            alpha = abs(math.asin(abs(u2)/math.sqrt((u1*u1)+(u2*u2)+(u3*u3))))

            alpha_deg = alpha/math.pi*180

            element.SetValue(CONTACT_ORIENTATION,alpha_deg)

            sigma = element.GetValue(CONTACT_SIGMA)

            OrientationChart.write(str(counter)+"    "+str(sigma/(self.total_stress_mean))+'\n')
            counter += 1

            if alpha_deg >= 0.0 and alpha_deg < 5.0:
                self.bond_00_05.append(element)

            if alpha_deg >= 5.0 and alpha_deg < 10.0:
                self.bond_05_10.append(element)

            if alpha_deg >= 10.0 and alpha_deg < 15.0:
                self.bond_10_15.append(element)

            if alpha_deg >= 15.0 and alpha_deg < 20.0:
                self.bond_15_20.append(element)

            if alpha_deg >= 20.0 and alpha_deg < 25.0:
                self.bond_20_25.append(element)

            if alpha_deg >= 25.0 and alpha_deg < 30.0:
                self.bond_25_30.append(element)

            if alpha_deg >= 30.0 and alpha_deg < 35.0:
                self.bond_30_35.append(element)

            if alpha_deg >= 35.0 and alpha_deg < 40.0:
                self.bond_35_40.append(element)

            if alpha_deg >= 40.0 and alpha_deg < 45.0:
                self.bond_40_45.append(element)

            if alpha_deg >= 45.0 and alpha_deg < 50.0:
                self.bond_45_50.append(element)

            if alpha_deg >= 50.0 and alpha_deg < 55.0:
                self.bond_50_55.append(element)

            if alpha_deg >= 55.0 and alpha_deg < 60.0:
                self.bond_55_60.append(element)

            if alpha_deg >= 60.0 and alpha_deg < 65.0:
                self.bond_60_65.append(element)

            if alpha_deg >= 65.0 and alpha_deg < 70.0:
                self.bond_65_70.append(element)

            if alpha_deg >= 70.0 and alpha_deg < 75.0:
                self.bond_70_75.append(element)

            if alpha_deg >= 75.0 and alpha_deg < 80.0:
                self.bond_75_80.append(element)

            if alpha_deg >= 80.0 and alpha_deg < 85.0:
                self.bond_80_85.append(element)

            if alpha_deg >= 85.0 and alpha_deg < 90.0:
                self.bond_85_90.append(element)

        ii=0
        for item in [self.bond_00_05, self.bond_05_10, self.bond_10_15, self.bond_15_20, self.bond_20_25, self.bond_25_30, self.bond_30_35, self.bond_35_40, self.bond_40_45,  self.bond_45_50, self.bond_50_55, self.bond_55_60, self.bond_60_65, self.bond_65_70, self.bond_70_75, self.bond_75_80, self.bond_80_85, self.bond_85_90]:

            self.sizes[ii] = len(item)

            i = 0.0
            sigma_sum =0.0
            tau_sum = 0.0

            sigma_total_sum_squared = 0
            tau_total_sum_squared = 0.0

            volume = 0.0
            area = 0.0

            for element in item:
                sigma_normal = element.GetValue(CONTACT_SIGMA)
                sigma_tau = element.GetValue(CONTACT_TAU)

                sigma_sum += sigma_normal
                tau_sum += sigma_tau

                sigma_partial_sum_squared = sigma_normal ** 2.0
                sigma_total_sum_squared += sigma_partial_sum_squared

                tau_partial_sum_squared = sigma_tau ** 2.0
                tau_total_sum_squared += tau_partial_sum_squared

                i += 1.0

            sigma_mean = sigma_sum / len(item)
            sigma_var = sigma_total_sum_squared / len(item) - sigma_mean ** 2.0

            sigma_std_dev = 0.0

            if abs(sigma_var) > 1e-9:
                std_dev = sigma_var ** 0.5

            sigma_rel_std_dev = sigma_std_dev / sigma_mean

            tau_mean = tau_sum/ len(item)
            tau_var = tau_total_sum_squared / len(item) - tau_mean ** 2.0

            tau_std_dev = 0.0

            if abs(tau_var) > 1e-9:
                tau_std_dev = tau_var ** 0.5

            tau_rel_std_dev = tau_std_dev / tau_mean

            self.sigma_mean_table[ii] = sigma_mean
            self.sigma_rel_std_dev_table[ii] = sigma_rel_std_dev
            self.tau_mean_table[ii] = tau_mean
            self.tau_rel_std_dev_table[ii] = tau_rel_std_dev
            self.sigma_ratio_table[ii]=sigma_mean/(self.total_stress_mean)
            ii+=1

        self.Procedures.KratosPrintInfo(self.sigma_ratio_table)
        OrientationChart.close()

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

    def MeasureRadialStrain(self):

        mean_radial_strain = 0.0
        radial_strain = 0.0
        weight = 0.0

        for node in self.XLAT:

            r = node.GetSolutionStepValue(RADIUS)
            x = node.X
            z = node.Z

            x0 = node.X0
            z0 = node.Z0

            dist_initial = math.sqrt(x0 * x0 + z0 * z0)
            dist_now = math.sqrt(x * x + z * z)
            node_radial_strain = (dist_now - dist_initial) / dist_initial
            mean_radial_strain += node_radial_strain

            weight += 1.0

        radial_strain = mean_radial_strain/weight

        return radial_strain

    def PoissonMeasure(self):

        self.Procedures.KratosPrintWarning("Not Working now")

        #left_nodes = list()
        #right_nodes = list()

        #xleft_weight  = 0.0
        #xright_weight  = 0.0

        #left_counter = 0.0
        #right_counter = 0.0

        #if(self.parameters.PoissonMeasure == "ON"):

            #for node in spheres_model_part.Nodes:

                #if (node.GetSolutionStepValue(GROUP_ID)==4):

                #left_nodes.append(node)
                #xleft_weight = +(node.X0 - node.GetSolutionStepValue(RADIUS))*node.GetSolutionStepValue(RADIUS)
                #left_counter = +node.GetSolutionStepValue(RADIUS)

                #elif(node.GetSolutionStepValue(GROUP_ID)==8):

                #right_nodes.append(node)
                #xright_weight = +(node.X + node.GetSolutionStepValue(RADIUS))*node.GetSolutionStepValue(RADIUS)
                #right_counter = +node.GetSolutionStepValue(RADIUS)

            #width_ini = xright_weight/right_counter - xleft_weight/left_counter

    ##################################POISSON##################################

        #if(self.parameters.PoissonMeasure == "ON"):

        #xleft_weight  = 0.0
        #xright_weight  = 0.0

        #left_counter = 0.0
        #right_counter = 0.0

        #for node in left_nodes:

            #xleft_weight = +(node.X - node.GetSolutionStepValue(RADIUS))*node.GetSolutionStepValue(RADIUS)
            #left_counter = +node.GetSolutionStepValue(RADIUS)

        #for node in right_nodes:

            #xright_weight = +(node.X + node.GetSolutionStepValue(RADIUS))*node.GetSolutionStepValue(RADIUS)
            #right_counter = +node.GetSolutionStepValue(RADIUS)

        #width_now = xright_weight/right_counter - xleft_weight/left_counter
        #measured_poisson =  ((width_now-width_ini)/width_ini)/strain

        #graph_export_poisson.write(str(strain)+"  "+str(measured_poisson)+'\n')

    #-------------------------------------------------------------------------------------#

    def GenerateGraphics(self):

        ## PROBLEM DATA
        area = 0.000001 ### 1mm2
        grad_p = 1 ## Pa/m

        ## Read Data
        data_file_name0 = "test.grf"
        data0 = loadtxt(data_file_name0)
        strain = array(data0[:,0])
        stress = array(data0[:,1])

        data_file_name1 = "test.grf"
        data1 = loadtxt(data_file_name1)
        strain1 = array(data1[:,0])
        stress1 = array(data1[:,1])

        data_file_name2 = "test.grf"
        data2 = loadtxt(data_file_name2)
        strain2 = array(data2[:,0])
        stress2 = array(data2[:,1])

        # setting to be changed#############################3
        set_mode = 'extralarge'  # large; publishable; medium
        legend_position = 'lower left'

        ##graph_name = ""
        x_name = 'Axial Strain (%)'
        y_name = 'Stress (MPa) - Load-axis'
        ####################################################################
        ####################################################################

        clf()
        plot_settings.set_mode(set_mode)
        #plt.semilogx()
        plot(strain, stress, 'k:s', strain1, stress1, 'r--v', strain2, stress2, 'b-.o',linewidth=1 )
        legend(('test', 'test'), legend_position, numpoints=1,)
        ##       bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2, mode="expand", borderaxespad=0.)
        grid(True)
        #insert name ######################################################
        savedname = "stress_graph"
        ####################################################################
        ##graphtitle = graph_name
        ##title(graphtitle)
        xlabel(x_name)
        ylabel(y_name)
        ##xlim(0.0, 1.0)
        ##ylim(0.0, 1.0)
        ##savefig(savedname + '.eps')
        savefig(savedname + '.png')

        ####################################################################
        ####################################################################

        clf()
        plot_settings.set_mode(set_mode)
        #plt.semilogx()
        plot(strain, stress, 'k:s', strain1, stress1, 'r--v',linewidth=2 )
        legend(( 'IFT variation', 'Viscosity variation'), legend_position, numpoints=1,)
        ##       bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2, mode="expand", borderaxespad=0.)
        grid(True)
        #insert name ######################################################
        savedname = "stress_graph2"
        ####################################################################
        ##graphtitle = graph_name
        ##title(graphtitle)
        xlabel(x_name)
        ylabel(y_name)
        ##xlim(0.0, 1.0)
        ##ylim(0.0, 1.0)
        ##savefig(savedname + '.eps')
        savefig(savedname + '.png')

class PreUtils():

    def __init__(self, spheres_model_part):

        self.spheres_model_part = spheres_model_part
        self.PreUtilities = PreUtilities()

    def BreakBondUtility(self, spheres_model_part):
        self.PreUtilities.BreakBondUtility(self.spheres_model_part)
