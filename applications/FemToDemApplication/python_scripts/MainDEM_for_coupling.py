
import KratosMultiphysics.DEMApplication as DEM
# import KratosMultiphysics.DEMApplication.main_script as MainDEM
import KratosMultiphysics.FemToDemApplication as FEMDEM
import KratosMultiphysics.DEMApplication.DEM_analysis_stage as MainDEM

class DEM_for_coupling_Solution(MainDEM.DEMAnalysisStage):

    def SetAnalyticParticleWatcher(self):
        pass

    def AddVariables(self):
        super(DEM_for_coupling_Solution, self).AddVariables()
        # For averaging forces when substepping
        self.spheres_model_part.AddNodalSolutionStepVariable(DEM.CONTACT_IMPULSE)
