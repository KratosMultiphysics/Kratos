import swimming_DEM_procedures as SDP
import swimming_DEM_algorithm
BaseAlgorithm = swimming_DEM_algorithm.Algorithm

class Algorithm(BaseAlgorithm):
    def __init__(self, pp):
        BaseAlgorithm.__init__(self, pp)

    def GetFluidSolveCounter(self, pp):
        return SDP.Counter(is_dead = True)

    def GetEmbeddedCounter(self, pp):
        return SDP.Counter(1, 3, pp.CFD_DEM.embedded_option)  # MA: because I think DISTANCE,1 (from previous time step) is not calculated correctly for step=1

    def GetBackwardCouplingCounter(self, pp):
        return SDP.Counter(1, 4, 0)

    def GetDebugInfo(self, pp):
        return SDP.Counter(pp.CFD_DEM.debug_tool_cycle, 1, is_dead = 1)
