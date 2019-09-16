from __future__ import print_function, absolute_import, division  #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.PfemFluidDynamicsApplication.pfem_fluid_dynamics_analysis as PfemFluidDynamicsAnalysis

def Wait():
    input("Press Something")

class MainPFEM_for_coupling_solution(PfemFluidDynamicsAnalysis):
    """The derived class for the PfemFluidDynamicsAnalysis
    """

    def OutputSolutionStep(self):
        """This function printed / writes output files after the solution of a step
        """
        pass

    def KratosPrintInfo(self, message):
        """This function prints info on screen
        """
        KratosMultiphysics.Logger.Print(message, label="")
        KratosMultiphysics.Logger.Flush()