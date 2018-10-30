from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

#Activate it to import in the gdb path:
#import sys
#sys.path.append('/home/cpuigbo/kratos')
#x = input("stopped to allow debug: set breakpoints and press enter to continue");
import time as timer
# Import system python
import os

# Import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.ExternalSolversApplication    as KratosSolvers
import KratosMultiphysics.DelaunayMeshingApplication    as KratosDelaunay
import KratosMultiphysics.PfemFluidDynamicsApplication  as KratosPfemFluid
import KratosMultiphysics.FluidDynamicsApplication as KratosFluidDynamics
import KratosMultiphysics.DEMApplication as Dem_application
import MainFluidPFEM

class Solution(MainFluidPFEM.Solution):

    def __init__(self, model):
        self.pp = self.ProblemParameters()
        super(Solution,self).__init__(model)

        self.fluid_model_part = self.main_model_part

        self.pp.nodal_results = []
        output_settings = self.ProjectParameters["output_configuration"]["result_file_configuration"]
        for i in range(output_settings["nodal_results"].size()):
            self.pp.nodal_results.append(output_settings["nodal_results"][i].GetString())

        self.pp.gauss_points_results = []
        for i in range(output_settings["gauss_point_results"].size()):
            self.pp.nodal_results.append(output_settings["gauss_point_results"][i].GetString())


        self.pp.problem_name = self.problem_name
        self.pp.Start_time = self.ProjectParameters["problem_data"]["start_time"].GetDouble()
        self.pp.VolumeOutput = True

        if output_settings["gidpost_flags"]["GiDPostMode"].GetString() == "GiD_PostBinary":
            self.pp.GiDPostMode  = "Binary"
        else:
            self.pp.GiDPostMode  = "Ascii"

        if output_settings["gidpost_flags"]["MultiFileFlag"].GetString() == "MultipleFiles":
            self.pp.GiDMultiFileFlag  = "Multiples"
        else:
            self.pp.GiDMultiFileFlag  = "Single"

        if output_settings["gidpost_flags"]["WriteDeformedMeshFlag"].GetString() == "WriteDeformed":
            self.pp.GiDWriteMeshFlag  = True
        else:
            self.pp.GiDWriteMeshFlag  = False

        if output_settings["gidpost_flags"]["WriteConditionsFlag"].GetString() == "WriteConditions":
            self.pp.GiDWriteConditionsFlag  = True
        else:
            self.pp.GiDWriteConditionsFlag  = False

    def AddNodalVariablesToModelPart(self):

        super(Solution,self).AddNodalVariablesToModelPart()

        # This is the Kinematic Viscosity
        self.fluid_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)

        self.AddFluidVariablesBySwimmingDEMAlgorithm()

    def CalculateNodalArea(self):

        KratosMultiphysics.CalculateNodalAreaProcess(self.main_model_part,self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]).Execute()

    def AddFluidVariablesBySwimmingDEMAlgorithm(self):
        self.vars_man.AddNodalVariables(self.main_model_part, self.pp.fluid_vars)

    def GetDeltaTimeFromParameters(self):
        return self.ProjectParameters["problem_data"]["time_step"].GetDouble()

    def GraphicalOutputExecuteInitialize(self):
        pass

    def GraphicalOutputExecuteBeforeSolutionLoop(self):
        pass

    def GraphicalOutputExecuteInitializeSolutionStep(self):
        pass

    def GraphicalOutputExecuteFinalizeSolutionStep(self):
        pass

    def GraphicalOutputPrintOutput(self):
        pass

    def GraphicalOutputExecuteFinalize(self):
        pass

    class ProblemParameters:
        def __init__(self):
            pass


if __name__ == "__main__":
    Solution().Run()

