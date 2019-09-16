from __future__ import print_function, absolute_import, division  #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import time as timer
import KratosMultiphysics.PfemFluidDynamicsApplication.pfem_fluid_dynamics_analysis as PfemFluidDynamicsAnalysis
import os

def Wait():
    input("Press Something")

class MainPFEM_for_coupling_solution(PfemFluidDynamicsAnalysis.PfemFluidDynamicsAnalysis):
    """The derived class for the PfemFluidDynamicsAnalysis
    """

    def __init__(self, model, parameters):
        """The constructor of the AnalysisStage-Object.
        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- The Model to be used
        parameters -- The ProjectParameters used
        """
        self.model = model
        #### TIME MONITORING START ####
        # Time control starts
        self.KratosPrintInfo(timer.ctime())
        # Measure process time
        self.t0p = timer.clock()
        # Measure wall time
        self.t0w = timer.time()
        #### TIME MONITORING END ####

        #### PARSING THE PARAMETERS ####

        #set echo level
        self.echo_level = parameters["problem_data"]["echo_level"].GetInt()

        # Print solving time
        self.report = False
        if( self.echo_level > 0 ):
            self.report = True

        self.KratosPrintInfo(" ")

        # defining the number of threads:
        num_threads = parameters["problem_data"]["threads"].GetInt()
        self.SetParallelSize(num_threads)
        self.KratosPrintInfo("::[KPFEM Simulation]:: [OMP USING" + str(num_threads) + "THREADS ]")
        #parallel.PrintOMPInfo()

        self.KratosPrintInfo(" ")
        self.KratosPrintInfo("::[KPFEM Simulation]:: [Time Step:" + str(parameters["solver_settings"]["time_stepping"]["time_step"].GetDouble()) + " echo:" +  str(self.echo_level) + "]")

        #### Model_part settings start ####
        super(PfemFluidDynamicsAnalysis.PfemFluidDynamicsAnalysis,self).__init__(model,parameters)
        # Defining the model_part
        self.main_model_part = self.model.GetModelPart(parameters["solver_settings"]["model_part_name"].GetString())

        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.SPACE_DIMENSION, parameters["solver_settings"]["domain_size"].GetInt())
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, parameters["solver_settings"]["domain_size"].GetInt())
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, parameters["solver_settings"]["time_stepping"]["time_step"].GetDouble())
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, parameters["problem_data"]["start_time"].GetDouble())
        if parameters["problem_data"].Has("gravity_vector"):
             self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.GRAVITY_X, parameters["problem_data"]["gravity_vector"][0].GetDouble())
             self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.GRAVITY_Y, parameters["problem_data"]["gravity_vector"][1].GetDouble())
             self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.GRAVITY_Z, parameters["problem_data"]["gravity_vector"][2].GetDouble())

        self.problem_path = os.getcwd()
        self.problem_name = parameters["problem_data"]["problem_name"].GetString()

    def OutputSolutionStep(self):
        """This function printed / writes output files after the solution of a step
        """
        pass

    def KratosPrintInfo(self, message):
        """This function prints info on screen
        """
        KratosMultiphysics.Logger.Print(message, label="")
        KratosMultiphysics.Logger.Flush()