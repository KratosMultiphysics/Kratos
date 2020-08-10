from __future__ import print_function, absolute_import, division  #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics as KM
import KratosMultiphysics.FemToDemApplication as FEMDEM
import KratosMultiphysics.PfemFluidDynamicsApplication as PFEM
import KratosMultiphysics.FemToDemApplication.MainCouplingFemDem_for_PFEM_coupling as MainCouplingFemDem_for_PFEM_coupling
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
class MainCouplingPfemFemDemAitken_Solution(MainCouplingPfemFemDem.MainCouplingPfemFemDem_Solution):
#============================================================================================================================

    def __init__(self, Model, PFEMparameters):
        super(MainCouplingPfemFemDemAitken_Solution, self).__init__(Model, PFEMparameters)

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
                min_relaxation = project_parameters["Aitken_parameters"]["min_relaxation"].GetDouble()
            else:
                min_relaxation = 0.2
            if (project_parameters["Aitken_parameters"].Has("initial_relaxation")):
                initial_relaxation = project_parameters["Aitken_parameters"]["initial_relaxation"].GetDouble()
            else:
                initial_relaxation = 0.825
        else:
            max_relaxation = 0.9
            min_relaxation = 0.1
            initial_relaxation = 0.825
            self.aitken_max_iterations = 10
            self.aitken_residual_dof_tolerance = 1e-5

        self.FSI_aitken_utility = FEMDEM.AitkenRelaxationFEMDEMUtility(initial_relaxation, max_relaxation, min_relaxation)

        self.developer_mode = False
        if self.developer_mode:
            self.pressure_plot = open("pressure_plot.txt", "w")
            self.pressure_plot.write("This File prints the pressures at the interface nodes!\n\n")
            self.pressure_plot.close()

#============================================================================================================================
    def SolveSolutionStep(self):

        aitken_iteration = 1
        residual_norm = 1.0
        residual_norm_old = 1.0
        is_converged = False

        self.FSI_aitken_utility.InitializeSolutionStep()
        solid_model_part = self.FEMDEM_Solution.FEM_Solution.main_model_part
        self.FSI_aitken_utility.InitializeInterfaceSubModelPart(solid_model_part)
        self.FSI_aitken_utility.ResetNodalValues(solid_model_part)

        # We perform the Aitken iterative scheme...
        while (not is_converged and aitken_iteration < self.aitken_max_iterations):

            KratosPrintInfo("")
            KratosPrintInfo(" ---> Aitken fsi iteration number: " + str(aitken_iteration))
            KratosPrintInfo("")

            residual_norm_old = residual_norm

            # It's necessary to Fix in order to maintain the FEMDEM Kinematics
            self.FixNodesModelPart(self.FEMDEM_Solution.FEM_Solution.main_model_part)

            KratosPrintInfo("==============================================" + "\n" +
                           " ==== SOLVING PFEM PART OF THE CALCULATION ====" + "\n" +
                           " ==============================================")
            self.SolveSolutionStepPFEM()

            if self.developer_mode:
                self.PrintPressures()
            
            # Now we Free the nodes to be calculated by the FEMDEM
            self.FreeNodesModelPart(self.FEMDEM_Solution.FEM_Solution.main_model_part)

            # Transfer pressure forces
            self.RegenerateAndUpdatePFEMPressureConditions()

            KratosPrintInfo("================================================" + "\n" +
                           " ==== SOLVING FEMDEM PART OF THE CALCULATION ====" + "\n" +
                           " ================================================")
            self.SolveSolutionStepFEMDEM()

            # If there are no interface nodes yet
            if (solid_model_part.GetSubModelPart("fsi_interface_model_part").NumberOfNodes() < 2):
                is_converged = True
                self.FSI_aitken_utility.FinalizeNonLinearIteration()
                break

            # We relax the obtained solid velocities "v" and give them to the fluid
            residual_norm = self.UpdateAndRelaxSolution(solid_model_part)

            self.FSI_aitken_utility.FinalizeNonLinearIteration()

            # Check convergence
            is_converged = self.CheckConvergence(residual_norm, residual_norm_old, aitken_iteration)
            if (aitken_iteration == 1): # We want at least 2 iterations
                is_converged = False

            aitken_iteration += 1

            if (not is_converged and aitken_iteration != self.aitken_max_iterations): # We reset the kinematics of fluid
                self.FSI_aitken_utility.ResetPFEMkinematicValues(self.PFEM_Solution.main_model_part)

        if (aitken_iteration == self.aitken_max_iterations):
            KratosPrintInfo(" Aitken reached max iterations!!!")

        self.FSI_aitken_utility.FinalizeSolutionStep()

        self.PFEM_Solution.main_model_part.RemoveNodes(KM.TO_ERASE)

        # We update the NO_MESH flag in the FEMDEM skin
        self.UpdateFEMDEMBoundary()

#============================================================================================================================
    def UpdateAndRelaxSolution(self, solid_model_part):

        self.FSI_aitken_utility.SavePreviousRelaxedValues(solid_model_part)

        vector_size = self.FSI_aitken_utility.GetVectorSize()
        iteration_value_vector = KM.Vector(vector_size)
        residual_value_vector  = KM.Vector(vector_size)
        residual_vector_norm   = 1.0

        self.FSI_aitken_utility.FillOldRelaxedValuesVector(solid_model_part, iteration_value_vector)
        residual_vector_norm = self.FSI_aitken_utility.ComputeInterfaceResidualVector(solid_model_part,
                                                                                      residual_value_vector)

        self.FSI_aitken_utility.UpdateSolution(residual_value_vector, iteration_value_vector)
        self.FSI_aitken_utility.UpdateInterfaceValues(solid_model_part, iteration_value_vector)

        return residual_vector_norm

#============================================================================================================================
    def CheckConvergence(self, residual_new, residual_old, iteration):
        is_conv = False
        error = residual_new / math.sqrt(self.FSI_aitken_utility.GetVectorSize())
        if (error < self.aitken_residual_dof_tolerance):
            is_conv = True
            if (iteration > 1):
                KratosPrintInfo(" Aitken converged with an error of " + str(error) + " and " + str(iteration) + " iterations.")
            else:
                KratosPrintInfo(" Aitken error of " + str(error))
            return is_conv
        else:
            KratosPrintInfo(" Aitken error of " + str(error))
        return is_conv

#============================================================================================================================
    # Just for developers
    def PrintPressures(self):
        self.pressure_plot = open("pressure_plot.txt","a")
        self.pressure_plot.write("STEP:" + str(self.FEMDEM_Solution.FEM_Solution.step) + "\n")

        for node in self.FEMDEM_Solution.FEM_Solution.main_model_part.GetSubModelPart("fsi_interface_model_part").Nodes:
            pressure = node.GetSolutionStepValue(KM.PRESSURE)
            text =  "   " + str(node.Id) + "   " + str(pressure) + "\n"
            self.pressure_plot.write(text)

        self.pressure_plot.close()



