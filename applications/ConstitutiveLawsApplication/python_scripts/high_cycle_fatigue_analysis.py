# Importing Kratos
import KratosMultiphysics
import KratosMultiphysics.ConstitutiveLawsApplication as CLA

# Importing the base class
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

class HighCycleFatigueAnalysis(StructuralMechanicsAnalysis):
    """This class is used to complement the structurea_mechanics_analysis
    when using the HCF constitutive law
    """

    def RunSolutionLoop(self):
        """This function executes the solution loop of the AnalysisStage
        It can be overridden by derived classes
        """
        while self.KeepAdvancingSolutionLoop():
            self.time = self._GetSolver().AdvanceInTime(self.time)
            process = CLA.AdvanceInTimeHighCycleFatigueProcess(self._GetSolver().GetComputingModelPart(), self.project_parameters)
            if self.project_parameters["fatigue"]["advancing_strategy"].GetBool():
                process.Execute()
                time_incr = self._GetSolver().GetComputingModelPart().ProcessInfo[CLA.TIME_INCREMENT]
                self.time += time_incr
                self._GetSolver().GetComputingModelPart().ProcessInfo[CLA.TIME_INCREMENT] = 0.0
            self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME] = self.time
            self.InitializeSolutionStep()
            self._GetSolver().Predict()
            self._GetSolver().SolveSolutionStep()
            self.FinalizeSolutionStep()
            self.OutputSolutionStep()


    def OutputSolutionStep(self):

        super(HighCycleFatigueAnalysis, self).OutputSolutionStep()
        interval = 0
        plot_output = self.project_parameters["fatigue"]["plot_output"].GetBool()
        if plot_output:
            if self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.STEP] == 1:
                first_code_line = True
                restart_code_line = False
                open("PlotElement.txt","w")
            elif self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
                first_code_line = False
                restart_code_line = True
            else:
                first_code_line = False
                restart_code_line = False

            id_for_print = self.project_parameters["fatigue"]["element_gausspoint_print"][0].GetInt()
            id_gauss_point = self.project_parameters["fatigue"]["element_gausspoint_print"][1].GetInt()
            self.main_model_part = self.model.GetModelPart(self.project_parameters["solver_settings"]["model_part_name"].GetString())
            
            output_interval_by_cycle = self.project_parameters["fatigue"]["output_interval_by_cycle"].GetBool()

            if output_interval_by_cycle:
                plot_file = open("PlotElement.txt","a")
                for elem in self.main_model_part.Elements:
                    if elem.Id == id_for_print:
                        if self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.STEP] == 1:
                            previous_number_of_cycles = 0
                        elif self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
                            previous_number_of_cycles = elem.CalculateOnIntegrationPoints(KratosMultiphysics.NUMBER_OF_CYCLES, self.main_model_part.ProcessInfo)[id_gauss_point] - 1
                        else:                        
                            previous_number_of_cycles = self.project_parameters["fatigue"]["previous_cycle"].GetInt()
                        
                        number_of_cycles = elem.CalculateOnIntegrationPoints(KratosMultiphysics.NUMBER_OF_CYCLES, self.main_model_part.ProcessInfo)
                        if number_of_cycles[id_gauss_point] > previous_number_of_cycles:
                        
                            local_number_of_cycles = elem.CalculateOnIntegrationPoints(CLA.LOCAL_NUMBER_OF_CYCLES, self.main_model_part.ProcessInfo)
                            damage = elem.CalculateOnIntegrationPoints(CLA.DAMAGE,self.main_model_part.ProcessInfo)
                            f_red = elem.CalculateOnIntegrationPoints(CLA.FATIGUE_REDUCTION_FACTOR, self.main_model_part.ProcessInfo)
                            reversion_factor = elem.CalculateOnIntegrationPoints(CLA.INFINITY_YIELD_STRESS,self.main_model_part.ProcessInfo)

                            if first_code_line == True:
                                plot_file.write("    " + "Element Id".rjust(20) + "    " + "{0:6d}".format(id_for_print).rjust(20) + "    " + "Gauss Point Id".rjust(20) + "    " + "{0:6d}".format(id_gauss_point).rjust(20) +  "\n")
                                plot_file.write("    " + "Time".rjust(20) + "    " + "Nc Global".rjust(20) + "    " + "Nc Local".rjust(20) + "    " + "Damage".rjust(20) + "    " + "Fred".rjust(20) + "    " + "Reversion Factor".rjust(20) + "\n")
                            elif restart_code_line == True: 
                                self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.IS_RESTARTED] = False
                                plot_file.write("\n" + "       RestartUtility Activated" + "\n")
                                plot_file.write("    " + "Element Id".rjust(20) + "    " + "{0:6d}".format(id_for_print).rjust(20) + "    " + "Gauss Point Id".rjust(20) + "    " + "{0:6d}".format(id_gauss_point).rjust(20) +  "\n")
                                plot_file.write("    " + "Time".rjust(20) + "    " + "Nc Global".rjust(20) + "    " + "Nc Local".rjust(20) + "    " + "Damage".rjust(20) + "    " + "Fred".rjust(20) + "    " + "Reversion Factor".rjust(20) + "\n")
                            plot_file.write("    " + "{0:.11e}".format(self.time).rjust(20) + "    " + "{0:.11e}".format(number_of_cycles[id_gauss_point]).rjust(20) + "    " + "{0:.11e}".format(local_number_of_cycles[id_gauss_point]).rjust(20) + "    " + "{0:.4e}".format(damage[id_gauss_point]).rjust(20) + "    " + "{0:.8e}".format(f_red[id_gauss_point]).rjust(20) + "    " + "{0:.4e}".format(reversion_factor[id_gauss_point]).rjust(20) + "\n")
                            plot_file.close()
                            self.project_parameters["fatigue"].AddEmptyValue("previous_cycle").SetInt(number_of_cycles[id_gauss_point])
                        break

            else:
                plot_file = open("PlotElement.txt","a")
                for elem in self.main_model_part.Elements:
                    if elem.Id == id_for_print:  
                        number_of_cycles = elem.CalculateOnIntegrationPoints(KratosMultiphysics.NUMBER_OF_CYCLES, self.main_model_part.ProcessInfo)
                        local_number_of_cycles = elem.CalculateOnIntegrationPoints(CLA.LOCAL_NUMBER_OF_CYCLES, self.main_model_part.ProcessInfo)
                        uniaxial_stresses = elem.CalculateOnIntegrationPoints(CLA.UNIAXIAL_STRESS,self.main_model_part.ProcessInfo)
                        damage = elem.CalculateOnIntegrationPoints(CLA.DAMAGE,self.main_model_part.ProcessInfo)
                        f_red = elem.CalculateOnIntegrationPoints(CLA.FATIGUE_REDUCTION_FACTOR, self.main_model_part.ProcessInfo)
                        reversion_factor = elem.CalculateOnIntegrationPoints(CLA.INFINITY_YIELD_STRESS,self.main_model_part.ProcessInfo)

                        if first_code_line == True:
                            plot_file.write("    " + "Element Id".rjust(20) + "    " + "{0:6d}".format(id_for_print).rjust(20) + "    " + "Gauss Point Id".rjust(20) + "    " + "{0:6d}".format(id_gauss_point).rjust(20) +  "\n")
                            plot_file.write("    " + "Time".rjust(20) + "    " + "Nc Global".rjust(20) + "    " + "Nc Local".rjust(20) + "    " + "Uniaxial Stresses".rjust(20) + "    " + "Damage".rjust(20) + "    " + "Fred".rjust(20) + "    " + "Reversion Factor".rjust(20) + "\n")
                        elif restart_code_line == True:
                            self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.IS_RESTARTED] = False
                            plot_file.write("\n" + "       RestartUtility Activated" + "\n")
                            plot_file.write("    " + "Element Id".rjust(20) + "    " + "{0:6d}".format(id_for_print).rjust(20) + "    " + "Gauss Point Id".rjust(20) + "    " + "{0:6d}".format(id_gauss_point).rjust(20) +  "\n")
                            plot_file.write("    " + "Time".rjust(20) + "    " + "Nc Global".rjust(20) + "    " + "Nc Local".rjust(20) + "    " + "Uniaxial Stresses".rjust(20) + "    " + "Damage".rjust(20) + "    " + "Fred".rjust(20) + "    " + "Reversion Factor".rjust(20) + "\n")
                        plot_file.write("    " + "{0:.11e}".format(self.time).rjust(20) + "    " + "{0:.11e}".format(number_of_cycles[id_gauss_point]).rjust(20) + "    " + "{0:.11e}".format(local_number_of_cycles[id_gauss_point]).rjust(20) + "    " + "{0:.14e}".format(uniaxial_stresses[id_gauss_point]).rjust(20) + "    " + "{0:.4e}".format(damage[id_gauss_point]).rjust(20) + "    " + "{0:.8e}".format(f_red[id_gauss_point]).rjust(20) + "    " + "{0:.4e}".format(reversion_factor[id_gauss_point]).rjust(20) + "\n")
                        plot_file.close()
                        break

                id_for_print = self.project_parameters["fatigue"]["node_print"].GetInt()
                plot_file = open("PlotNode.txt","a")   
                for node in self.main_model_part.Nodes:
                    if node.Id == id_for_print:
                        if first_code_line == True:
                            plot_file.write("    " + "Node Id".rjust(20) + "    " + "{0:6d}".format(id_for_print).rjust(20) +  "\n")
                            plot_file.write("    " + "Time".rjust(20) + "    " + "Displacement".rjust(20) + "\n")                        
                        if restart_code_line == True:
                            plot_file.write("\n" + "       RestartUtility Activated" + "\n")
                            plot_file.write("    " + "Node Id".rjust(20) + "    " + "{0:6d}".format(id_for_print).rjust(20) +  "\n")
                            plot_file.write("    " + "Time".rjust(20) + "    " + "Displacement".rjust(20) + "\n")   
                        displacement = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)
                        plot_file.write("    " + "{0:.11e}".format(self.time).rjust(20) + "    " + "{0:.5e}".format(displacement[0]).rjust(20) + "\n")
                        plot_file.close()
                        break
