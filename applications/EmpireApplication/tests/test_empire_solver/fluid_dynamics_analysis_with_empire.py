from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication
import KratosMultiphysics.EmpireApplication

from fluid_dynamics_analysis import FluidDynamicsAnalysis

import sys
import time

USE_EMPIRE = True # Flag for testing, if set to False then this is a pure CFD-Analysis

class FluidDynamicsAnalysisWithEmpire(FluidDynamicsAnalysis):
    '''This class is for testing only
    It enhances the standard FluidDynamicsAnalysis to work with /
    test the EmpireSolver
    '''

    def Initialize(self):
        super(FluidDynamicsAnalysisWithEmpire,self).Initialize()

        if USE_EMPIRE:
            print(self.__class__.__name__ + ":" + " Starting to initialize Empire")
            import empire_wrapper
            print(self.__class__.__name__ + ":" + " Wrapper-Import Successful")
            self.empire = empire_wrapper.EmpireWrapper()
            print(self.__class__.__name__ + ":" + " Wrapper Created")
            xml_file_name = self.project_parameters["xml_file_name"].GetString()
            self.empire.Connect(xml_file_name)

            self.empire.SendMesh("Structural_Mesh_FSI_Interface",
                                 self.model[self.project_parameters["solver_settings"]["model_part_name"].GetString()])

    def RunSolutionLoop(self):
        """This function executes the solution loop of the AnalysisStage
        It can be overridden by derived classes
        """
        while self.time < self.end_time:
            self.time = self._GetSolver().AdvanceInTime(self.time)
            self.InitializeSolutionStep()
            self._GetSolver().Predict()

            if USE_EMPIRE:
                print("BEFFOORREE receiving")
                self.empire.ReceiveDataField("Structural_Mesh_FSI_Interface",
                                            "displacement",
                                            KratosMultiphysics.DISPLACEMENT)
                print("Afterrrrrr receiving")

            self._GetSolver().SolveSolutionStep()

            if USE_EMPIRE:
                print("BEFFOORREE SENDING")
                self.empire.SendDataField("Structural_Mesh_FSI_Interface",
                                          "displacement",
                                          KratosMultiphysics.DISPLACEMENT)
                print("Afterrrrrr SENDING")

            self.FinalizeSolutionStep()
            self.OutputSolutionStep()

if __name__ == "__main__":
    with open("ProjectParametersCFD.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = FluidDynamicsAnalysisWithEmpire(model,parameters)
    simulation.Run()
