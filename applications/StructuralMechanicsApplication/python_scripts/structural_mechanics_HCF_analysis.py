from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing Kratos
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as SMA

# Importing the base class
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

class StructuralMechanicsHCFAnalysis(StructuralMechanicsAnalysis):
    """
    This class is used to complement the structurea_mechanics_analysis 
    
    when using the HCF constitutive law
    """
    def __init__(self, model, project_parameters):
        super(StructuralMechanicsHCFAnalysis, self).__init__(model, project_parameters)

    def RunSolutionLoop(self):
        """This function executes the solution loop of the AnalysisStage
        It can be overridden by derived classes
        """
        while self.KeepAdvancingSolutionLoop():
            self.time = self._GetSolver().AdvanceInTime(self.time)
            process = SMA.AdvanceInTimeStrategyHighCycleFatigueProcess(self._GetSolver().GetComputingModelPart(), self.project_parameters)
            if self.project_parameters["fatigue"]["advancing_strategy"].GetBool():
                process.Execute()
                time_incr = self._GetSolver().GetComputingModelPart().ProcessInfo[SMA.TIME_INCREMENT]
                self.time += time_incr
                self._GetSolver().GetComputingModelPart().ProcessInfo[SMA.TIME_INCREMENT] = 0.0
            self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME] = self.time
            self.InitializeSolutionStep()
            self._GetSolver().Predict()
            self._GetSolver().SolveSolutionStep()
            self.FinalizeSolutionStep()
            self.OutputSolutionStep()

    def FinalizeSolutionStep(self):
        """This function performs all the required operations that should be executed
        (for each step) AFTER solving the solution step.
        """
        self._GetSolver().FinalizeSolutionStep()

        for process in self._GetListOfProcesses():
            process.ExecuteFinalizeSolutionStep()

        interval = 0
        if self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.STEP] == 1:
            self.TimePreviousPlotting=0

        if self.time - self.TimePreviousPlotting >= interval:
            if self.TimePreviousPlotting == 0:
                first_code_line = True
            else:
                first_code_line = False

            self.TimePreviousPlotting = self.time
            id_for_print = self.project_parameters["fatigue"]["element_gausspoint_print"][0].GetInt()
            id_gauss_point = self.project_parameters["fatigue"]["element_gausspoint_print"][1].GetInt()
            self.main_model_part = self.model.GetModelPart(self.project_parameters["solver_settings"]["model_part_name"].GetString())
            for elem in self.main_model_part.Elements:
                if elem.Id == id_for_print:
                    plot_file = open("PlotElement.txt","a")
                    number_of_cycles = elem.CalculateOnIntegrationPoints(SMA.NUMBER_OF_CYCLES, self.main_model_part.ProcessInfo)
                    local_number_of_cycles = elem.CalculateOnIntegrationPoints(SMA.LOCAL_NUMBER_OF_CYCLES, self.main_model_part.ProcessInfo)
                    uniaxial_stresses = elem.CalculateOnIntegrationPoints(SMA.UNIAXIAL_STRESS,self.main_model_part.ProcessInfo)
                    damage = elem.CalculateOnIntegrationPoints(SMA.DAMAGE,self.main_model_part.ProcessInfo)
                    f_red = elem.CalculateOnIntegrationPoints(SMA.FATIGUE_REDUCTION_FACTOR, self.main_model_part.ProcessInfo)
                    stress_tensor = elem.CalculateOnIntegrationPoints(KratosMultiphysics.CAUCHY_STRESS_TENSOR, self.main_model_part.ProcessInfo)
                    strain_vector = elem.CalculateOnIntegrationPoints(KratosMultiphysics.GREEN_LAGRANGE_STRAIN_TENSOR, self.main_model_part.ProcessInfo)
                    wohler_stress = elem.CalculateOnIntegrationPoints(SMA.WOHLER_STRESS, self.main_model_part.ProcessInfo)
                    yield_stress = self.main_model_part.GetProperties()[1][KratosMultiphysics.YIELD_STRESS]
                    
                    if first_code_line == True:
                        plot_file.write("    " + "id_for_print".rjust(20) + "    " + "{0:.4e}".format(id_for_print).rjust(20) +  "\n")    
                        plot_file.write("    " + "self.time".rjust(20) + "    " + "N_c GLOBAL".rjust(20) + "    " + "N_c LOCAL".rjust(20) + "    " + "uniaxial_stresses".rjust(20) + "    " + "damage".rjust(20) + "    " + "f_red".rjust(20) + "    " + "wohler_stress".rjust(20) + "    " + "strain_vector".rjust(20) + "    " + "stress_tensor".rjust(20) + "    " + "stress_tensor_norm".rjust(20) + "    " + "uniaxial_stress_norm".rjust(20) + "\n")    
                    
                    plot_file.write("    " + "{0:.11e}".format(self.time).rjust(20) + "    " + "{0:.11e}".format(number_of_cycles[id_gauss_point]).rjust(20) + "    " + "{0:.11e}".format(local_number_of_cycles[id_gauss_point]).rjust(20) + "    " + "{0:.14e}".format(uniaxial_stresses[id_gauss_point]).rjust(20) + "    " + "{0:.4e}".format(damage[id_gauss_point]).rjust(20) + "    " + "{0:.8e}".format(f_red[id_gauss_point]).rjust(20) + "    " + "{0:.4e}".format(wohler_stress[id_gauss_point]).rjust(20) + "    " + "{0:.4e}".format(strain_vector[id_gauss_point][0,0]).rjust(20) + "   " + "{0:.4e}".format(stress_tensor[id_gauss_point][0,0]).rjust(20) + "   " + "{0:.4e}".format(stress_tensor[id_gauss_point][0,0] / yield_stress).rjust(20) + "   " + "{0:.4e}".format(uniaxial_stresses[id_gauss_point] / yield_stress).rjust(20) + "\n")
                    plot_file.close()

            id_for_print = self.project_parameters["fatigue"]["node_print"].GetInt()
            for node in self.main_model_part.Nodes:
                if node.Id == id_for_print:
                    plot_file = open("PlotNode.txt","a")
                    displacement = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)

                    plot_file.write("    " + "{0:.4e}".format(self.time).rjust(11) + "    " + "{0:.4e}".format(displacement[0]).rjust(11) + "\n")
                    plot_file.close()