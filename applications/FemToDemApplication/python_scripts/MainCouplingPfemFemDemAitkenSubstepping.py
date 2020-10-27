
import KratosMultiphysics as KM
import KratosMultiphysics.FemToDemApplication as FEMDEM
import KratosMultiphysics.PfemFluidDynamicsApplication as PFEM
import KratosMultiphysics.FemToDemApplication.MainCouplingFemDemSubstepping_for_PFEM_coupling as MainCouplingFemDemSubstepping_for_PFEM_coupling
import KratosMultiphysics.FemToDemApplication.MainPFEM_for_coupling as MainPFEM_for_coupling
import KratosMultiphysics.FemToDemApplication.MainCouplingPfemFemDem as MainCouplingPfemFemDem
import math as math

def Wait():
    input("PFEM-FEMDEM Aitken -> Press Something")

def KratosPrintInfo(message):
    """This function prints info on screen
    """
    KM.Logger.Print("", message)
    KM.Logger.Flush()

#============================================================================================================================
class MainCouplingPfemFemDemAitkenSubstepping_Solution(MainCouplingPfemFemDemAitken.MainCouplingPfemFemDemAitken_Solution):
#============================================================================================================================

    def __init__(self, Model, PFEMparameters):
        # Initialize solutions of the FEMDEM and PFEM
        self.model = Model
        self.FEMDEM_Solution = MainCouplingFemDemSubstepping_for_PFEM_coupling.MainCoupledFemDemSubstepping_for_PFEM_coupling_Solution(Model)
        self.FEMDEM_Solution.is_slave = True
        self.FEMDEM_Solution.Initialize()

        self.PFEM_Solution = MainPFEM_for_coupling.MainPFEM_for_coupling_solution(Model, 
                                                                                  self.FEMDEM_Solution.FEM_Solution.main_model_part,
                                                                                  PFEMparameters)
        KratosPrintInfo("    ___                  _          _            _ _   _         ___  ___  __       "  + "\n" +
                       "    / __\___  _   _ _ __ | | ___  __| | __      _(_) |_| |__     / _ \/ __\/__\/\/\   " + "\n" +
                       "   / /  / _ \| | | | '_ \| |/ _ \/ _` | \ \ /\ / / | __| '_ \   / /_)/ _\ /_\ /    \  " + "\n" +
                       "  / /__| (_) | |_| | |_) | |  __/ (_| |  \ V  V /| | |_| | | | / ___/ /  //__/ /\/\ \ " + "\n" +
                       "  \____/\___/ \__,_| .__/|_|\___|\__,_|   \_/\_/ |_|\__|_| |_| \/   \/   \__/\/    \/ " + "\n")

        project_parameters = self.FEMDEM_Solution.FEM_Solution.ProjectParameters
        if (project_parameters.Has("Aitken_parameters")):
            if (project_parameters["Aitken_parameters"].Has("tolerance")):
                self.aitken_residual_dof_tolerance = project_parameters["Aitken_parameters"]["tolerance"].GetDouble()
            else:
                self.aitken_residual_dof_tolerance = 1.0e-7
            if (project_parameters["Aitken_parameters"].Has("max_iterations")):
                self.aitken_max_iterations = project_parameters["Aitken_parameters"]["max_iterations"].GetInt()
            else:
                self.aitken_max_iterations = 10
            if (project_parameters["Aitken_parameters"].Has("max_relaxation")):
                max_relaxation = project_parameters["Aitken_parameters"]["max_relaxation"].GetDouble()
            else:
                max_relaxation = 0.9
            if (project_parameters["Aitken_parameters"].Has("min_relaxation")):
                min_relaxation = project_parameters["Aitken_parameters"]["min_iterations"].GetDouble()
            else:
                min_relaxation = 0.2
            if (project_parameters["Aitken_parameters"].Has("initial_relaxation")):
                initial_relaxation = project_parameters["Aitken_parameters"]["initial_relaxation"].GetDouble()
            else:
                initial_relaxation = 0.825
        else:
            max_relaxation = 0.95
            min_relaxation = 0.1
            initial_relaxation = 0.825
            self.aitken_max_iterations = 20
            self.aitken_residual_dof_tolerance = 1e-7

        self.FSI_aitken_utility = FEMDEM.AitkenRelaxationFEMDEMUtility(initial_relaxation, max_relaxation, min_relaxation)

        self.developer_mode = False
        if self.developer_mode:
            self.pressure_plot = open("pressure_plot.txt", "w")
            self.pressure_plot.write("This File prints the pressures at the interface nodes!\n\n")
            self.pressure_plot.close()