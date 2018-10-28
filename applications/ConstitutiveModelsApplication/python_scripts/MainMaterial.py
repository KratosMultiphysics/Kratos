from __future__ import print_function, absolute_import, division # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Import system python modules
import time as timer
import sys
import os

# Import kratos core and applications
import KratosMultiphysics

class Solution(object):

    def __init__(self, input_file = "ProjectParameters.json", materials_file = "Materials.json"):

        # Time control starts
        print(timer.ctime())

        sys.stdout.flush()

        # Measure process time
        self.t0p = timer.clock()

        # Measure wall time
        self.t0w = timer.time()

        # Read input parameters
        settings = None
        if os.path.isfile(input_file):
            parameters_file = open(input_file,'r')
            settings = KratosMultiphysics.Parameters(parameters_file.read())
        else:
            raise Exception("Input file "+input_file+" does not exist")

        # Read material parameters
        parameters = None
        if os.path.isfile(materials_file):
            parameters_file = open(materials_file,'r')
            parameters = KratosMultiphysics.Parameters(parameters_file.read())
        else:
            raise Exception("Materials file "+materials_file+" does not exist")

        # Set materials
        if( parameters.Has("material_models_list") ):
            self.materials = parameters["material_models_list"]
        else:
            raise Exception("No material model supplied")

        # Set solver process
        if( settings.Has("solver_process") ):
            self.solver_settings = settings["solver_process"]
        else:
            raise Exception("Solver process not supplied")

        # Print solving time
        self.report = False
        self.echo_level = 0
        if( self.echo_level > 0 ):
            self.report = True

    def Run(self):

        self.Initialize()

        self.Solve()

        self.Finalize()


    def Initialize(self):

        # start model
        model_part_name = self.materials[0]["Parameters"]["model_part_name"].GetString()
        self.model = KratosMultiphysics.Model()
        self.model_part = self.model.CreateModelPart(model_part_name)

        # construct material process
        self._construct_material_processes(self.materials)

        # set solver
        self.solver = self._get_solver_process(self.solver_settings)

        # set time parameters
        self.process_info = self.model_part.ProcessInfo

        # Set time settings
        self.step       = self.process_info[KratosMultiphysics.STEP]
        self.time       = self.process_info[KratosMultiphysics.TIME]

        self.delta_time = self.process_info[KratosMultiphysics.DELTA_TIME]
        self.end_time   = self.solver.GetEndTime()

        # Initialize Solver
        self.solver.Initialize()

        print(" ")
        print("::[Material Modelling]:: -START- ")


    def Solve(self):

        # Solving the problem (integration)
        while(self.time < self.end_time):

            self.InitializeSolutionStep()
            self.SolveSolutionStep()
            self.FinalizeSolutionStep()

            sys.stdout.flush()


    def InitializeSolutionStep(self):

        self.clock_time = self._start_time_measuring();

        self.time = self.time + self.delta_time

        self.process_info[KratosMultiphysics.TIME] = self.time

        print(" [ MODEL STEP - TIME :","{0:1.{1}f}".format(self.time,6),"]")

        self.solver.InitializeSolutionStep()

        self._stop_time_measuring(self.clock_time,"Initialize Step", self.report);


    def SolveSolutionStep(self):

        self.clock_time = self._start_time_measuring();

        self.solver.Solve()

        self._stop_time_measuring(self.clock_time,"Solve Step", self.report);


    def FinalizeSolutionStep(self):

        self.clock_time = self._start_time_measuring();

        self.solver.FinalizeSolutionStep()

        self._stop_time_measuring(self.clock_time,"Finalize Step", self.report);


    def Finalize(self):

        # Ending the problem
        self.solver.Finalize()

        print("::[Material Modelling]:: -END- ")
        print(" ")

        # Measure process time
        tfp = timer.clock()

        # Measure wall time
        tfw = timer.time()

        print("::[Material Modeling]:: [Elapsed Time = %.2f" % (tfw - self.t0w),"seconds] (%.2f" % (tfp - self.t0p),"seconds of cpu/s time)")
        print(timer.ctime())

    #### Main internal methods ####


    def _get_solver_process(self, solver_settings):
        solver_module = __import__(solver_settings["solver_type"].GetString())
        return (solver_module.CreateSolver(self.model_part, solver_settings["Parameters"]))


    def _construct_material_processes(self, material_parameters):
        import process_factory
        return (process_factory.KratosProcessFactory(self.model).ConstructListOfProcesses(material_parameters))

    @classmethod
    def _start_time_measuring(self):
        # Measure process time
        time_ip = timer.clock()
        return time_ip

    @classmethod
    def _stop_time_measuring(self, time_ip, process, report):
        # Measure process time
        time_fp = timer.clock()
        if( report ):
            used_time = time_fp - time_ip
            print("::[Material Modelling]:: [  %.2f" % round(used_time,2),"s", process," ] ")

if __name__ == "__main__":
    Solution().Run()
