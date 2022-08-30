
import KratosMultiphysics.DEMApplication as DEM
import KratosMultiphysics.DEMApplication.DEM_analysis_stage as MainDEM

class DEM_for_coupling_Solution(MainDEM.DEMAnalysisStage):

    def SetAnalyticParticleWatcher(self):
        pass

    def AddVariables(self):
        super(DEM_for_coupling_Solution, self).AddVariables()
        # For averaging forces when substepping
        self.spheres_model_part.AddNodalSolutionStepVariable(DEM.CONTACT_IMPULSE)

    def GraphicalOutputInitialize(self):
        pass

    def PrintResultsForGid(self, time):
        pass

    def GraphicalOutputFinalize(self):
        pass

    def PrintResults(self):
        pass

    def RunAnalytics(self, time, is_time_to_print=True):
        pass