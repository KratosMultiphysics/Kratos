from KratosMultiphysics.CoSimulationApplication.solver_wrappers.kratos import structural_mechanics_wrapper
import KratosMultiphysics

class VolumeCouplingStructuralWrapper(structural_mechanics_wrapper.StructuralMechanicsWrapper):

    def Initialize(self):
        # Set ACTIVATION_LEVEL to 0 at the start
        self._analysis_stage._solver.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.ACTIVATION_LEVEL,0)
        #print("ACTIVATION_LEVEL INITIALIZE ################################################ = ", self._analysis_stage._solver.main_model_part.ProcessInfo[KratosMultiphysics.ACTIVATION_LEVEL])
        super().Initialize()

    def OutputSolutionStep(self):
        # Manipulate ACTIVATION_LEVEL as required
        self._analysis_stage._solver.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.ACTIVATION_LEVEL,1)
        #print("before output_ACTIVATION_LEVEL = ", self._analysis_stage._solver.main_model_part.ProcessInfo[KratosMultiphysics.ACTIVATION_LEVEL])
        super().OutputSolutionStep()
        self._analysis_stage._solver.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.ACTIVATION_LEVEL,0)
        # print activation level
        #print("afteroutput_ACTIVATION_LEVEL = ", self._analysis_stage._solver.main_model_part.ProcessInfo[KratosMultiphysics.ACTIVATION_LEVEL])

def Create(settings, model, solver_name):
    return VolumeCouplingStructuralWrapper(settings, model, solver_name)

