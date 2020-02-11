from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Import system python modules
import time as timer
import os

# Import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication     as KratosSolid
import KratosMultiphysics.ExternalSolversApplication as KratosSolvers
import KratosMultiphysics.FemToDemApplication as KratosFemDem
import KratosMultiphysics.FemToDemApplication.MainSolidFEM as MainSolidFEM
import KratosMultiphysics.process_factory as process_factory
import KratosMultiphysics.gid_output_process as gid_output_process


def Wait():
    input("Press Something")

class FEM_Solution(MainSolidFEM.Solution):

    def Info(self):
        KratosMultiphysics.Logger.PrintInfo("FEM part of the FEMDEM application")

    def KratosPrintInfo(self, message):
        KratosMultiphysics.Logger.Print(message, label="")
        KratosMultiphysics.Logger.Flush()
#============================================================================================================================                    
    def __init__(self, Model):

        #### TIME MONITORING START ####
        # Time control starts        
        self.KratosPrintInfo(timer.ctime())
        # Measure process time
        self.t0p = timer.clock()
        # Measure wall time
        self.t0w = timer.time()
        #### TIME MONITORING END ####

        #### PARSING THE PARAMETERS ####

        # Import input
        parameter_file = open("ProjectParameters.json",'r')
        self.ProjectParameters = KratosMultiphysics.Parameters(parameter_file.read())

        # set echo level
        self.echo_level = self.ProjectParameters["problem_data"]["echo_level"].GetInt()

        self.KratosPrintInfo(" ")

        # defining the number of threads:
        num_threads =  self.GetParallelSize()
        self.KratosPrintInfo("::[KSM Simulation]:: [OMP USING " + str(num_threads) +" THREADS ]")
        #parallel.PrintOMPInfo()

        #### Model_part settings start ####

        # Defining the model_part
        self.model = Model
        self.model.CreateModelPart(self.ProjectParameters["problem_data"]["model_part_name"].GetString())
        self.main_model_part = self.model.GetModelPart(self.ProjectParameters["problem_data"]["model_part_name"].GetString())

        if (self.ProjectParameters["solver_settings"]["solution_type"].GetString() == "Dynamic"):
            self.main_model_part.ProcessInfo.SetValue(KratosFemDem.IS_DYNAMIC, 1)
        else:
            self.main_model_part.ProcessInfo.SetValue(KratosFemDem.IS_DYNAMIC, 0)

        self.time_step  = self.ComputeDeltaTime()
        self.start_time = self.ProjectParameters["problem_data"]["start_time"].GetDouble()
        self.end_time   = self.ProjectParameters["problem_data"]["end_time"  ].GetDouble()

        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, self.ProjectParameters["problem_data"]["domain_size"].GetInt())
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, self.time_step)
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, self.start_time)


        ### replace this "model" for real one once available in kratos core
        self.Model = {self.ProjectParameters["problem_data"]["model_part_name"].GetString() : self.main_model_part}

        # Construct the solver (main setting methods are located in the solver_module)
        if self.ProjectParameters["solver_settings"]["solver_type"].GetString() == "FemDemDynamicSolver":
            import KratosMultiphysics.FemToDemApplication.FemDemDynamicSolver as FemDemDynamicSolver
            self.solver = FemDemDynamicSolver.CreateSolver(self.main_model_part, self.ProjectParameters["solver_settings"])
        elif self.ProjectParameters["solver_settings"]["solver_type"].GetString() == "FemDemStaticSolver":
            import KratosMultiphysics.FemToDemApplication.FemDemStaticSolver as FemDemStaticSolver
            self.solver = FemDemStaticSolver.CreateSolver(self.main_model_part, self.ProjectParameters["solver_settings"])

        #### Output settings start ####
        self.problem_path = os.getcwd()
        self.problem_name = self.ProjectParameters["problem_data"]["problem_name"].GetString()

        
#============================================================================================================================
    def AddMaterials(self):

        # Assign material to model_parts (if Materials.json exists)
        if os.path.isfile("Materials.json"):
            materials_file = open("Materials.json",'r')
            MaterialParameters = KratosMultiphysics.Parameters(materials_file.read())

            if(MaterialParameters.Has("material_models_list")):

                ## Get the list of the model_part's in the object Model
                for i in range(self.ProjectParameters["solver_settings"]["problem_domain_sub_model_part_list"].size()):
                    part_name = self.ProjectParameters["solver_settings"]["problem_domain_sub_model_part_list"][i].GetString()
                    if( self.main_model_part.HasSubModelPart(part_name) ):
                        self.Model.update({part_name: self.main_model_part.GetSubModelPart(part_name)})

                assign_materials_processes = process_factory.KratosProcessFactory(self.Model).ConstructListOfProcesses( MaterialParameters["material_models_list"] )
            for process in assign_materials_processes:
                process.Execute()
        else:
            self.KratosPrintInfo(" No Materials.json found ")
                
#============================================================================================================================          
    def AddProcesses(self):

        # Build sub_model_parts or submeshes (rearrange parts for the application of custom processes)
        ## Get the list of the submodel part in the object Model
        for i in range(self.ProjectParameters["solver_settings"]["processes_sub_model_part_list"].size()):
            part_name = self.ProjectParameters["solver_settings"]["processes_sub_model_part_list"][i].GetString()
            if( self.main_model_part.HasSubModelPart(part_name) ):
                self.Model.update({part_name: self.main_model_part.GetSubModelPart(part_name)})
        
        # Obtain the list of the processes to be applied
        import KratosMultiphysics.SolidMechanicsApplication.process_handler

        process_parameters = KratosMultiphysics.Parameters("{}") 
        process_parameters.AddValue("echo_level", self.ProjectParameters["problem_data"]["echo_level"])
        process_parameters.AddValue("constraints_process_list", self.ProjectParameters["constraints_process_list"])
        process_parameters.AddValue("loads_process_list", self.ProjectParameters["loads_process_list"])
        if( self.ProjectParameters.Has("problem_process_list") ):
            process_parameters.AddValue("problem_process_list", self.ProjectParameters["problem_process_list"])
        if( self.ProjectParameters.Has("output_process_list") ):
            process_parameters.AddValue("output_process_list", self.ProjectParameters["output_process_list"])

        return (KratosMultiphysics.SolidMechanicsApplication.process_handler.ProcessHandler(self.Model, process_parameters))

#============================================================================================================================    
    def Run(self):

        self.Initialize()
        self.RunMainTemporalLoop()
        self.Finalize()    
#============================================================================================================================        
    def Initialize(self):

        # Add variables (always before importing the model part)
        self.solver.AddVariables()
        
        # Read model_part (note: the buffer_size is set here) (restart is read here)
        self.solver.ImportModelPart()

        # Add dofs (always after importing the model part)
        if((self.main_model_part.ProcessInfo).Has(KratosMultiphysics.IS_RESTARTED)):
            if(self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] == False):
                self.solver.AddDofs()
        else:
            self.solver.AddDofs()

        # Add materials (assign material to model_parts if Materials.json exists)
        self.AddMaterials()
        
        # Add processes
        self.model_processes = self.AddProcesses()
        self.model_processes.ExecuteInitialize()

        # Print model_part and properties
        if(self.echo_level > 1):
            self.KratosPrintInfo("")
            self.KratosPrintInfo(self.main_model_part)
            for properties in self.main_model_part.Properties:
                self.KratosPrintInfo(properties)

        #### START SOLUTION ####
        self.computing_model_part = self.solver.GetComputingModelPart()

        ## Sets strategies, builders, linear solvers, schemes and solving info, and fills the buffer
        self.solver.Initialize()
        self.solver.SetEchoLevel(self.echo_level)

        # Initialize GiD  I/O (gid outputs, file_lists)
        self.SetGraphicalOutput()
        
        self.GraphicalOutputExecuteInitialize()

        self.model_processes.ExecuteBeforeSolutionLoop()

        self.GraphicalOutputExecuteBeforeSolutionLoop()        

        # Set time settings
        self.step       = self.main_model_part.ProcessInfo[KratosMultiphysics.STEP]
        self.time       = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]

        self.end_time   = self.ProjectParameters["problem_data"]["end_time"].GetDouble()
        self.delta_time = self.ProjectParameters["problem_data"]["time_step"].GetDouble()

#============================================================================================================================    
    def Run(self):

        self.Initialize()
        self.RunMainTemporalLoop()
        self.Finalize()    
#============================================================================================================================        
    def Initialize(self):

        # Add variables (always before importing the model part)
        self.solver.AddVariables()
        
        # Read model_part (note: the buffer_size is set here) (restart is read here)
        self.solver.ImportModelPart()

        # Add dofs (always after importing the model part)
        if((self.main_model_part.ProcessInfo).Has(KratosMultiphysics.IS_RESTARTED)):
            if(self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] == False):
                self.solver.AddDofs()
        else:
            self.solver.AddDofs()

        # Add materials (assign material to model_parts if Materials.json exists)
        self.AddMaterials()
        
        # Add processes
        self.model_processes = self.AddProcesses()
        self.model_processes.ExecuteInitialize()

        # Print model_part and properties
        if(self.echo_level > 1):
            self.KratosPrintInfo("")
            self.KratosPrintInfo(self.main_model_part)
            for properties in self.main_model_part.Properties:
                self.KratosPrintInfo(properties)

        #### START SOLUTION ####
        self.computing_model_part = self.solver.GetComputingModelPart()

        ## Sets strategies, builders, linear solvers, schemes and solving info, and fills the buffer
        self.solver.Initialize()
        self.solver.SetEchoLevel(self.echo_level)

        # Initialize GiD  I/O (gid outputs, file_lists)
        self.SetGraphicalOutput()
        
        self.GraphicalOutputExecuteInitialize()

        self.model_processes.ExecuteBeforeSolutionLoop()

        self.GraphicalOutputExecuteBeforeSolutionLoop()        

        # Set time settings
        self.step       = self.main_model_part.ProcessInfo[KratosMultiphysics.STEP]
        self.time       = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]

        self.end_time   = self.ProjectParameters["problem_data"]["end_time"].GetDouble()
        self.delta_time = self.ProjectParameters["problem_data"]["time_step"].GetDouble()


#============================================================================================================================
    def RunMainTemporalLoop(self):

        # Solving the problem (time integration)
        while(self.time < self.end_time):
            
            self.InitializeSolutionStep()
            self.SolveSolutionStep()
            self.FinalizeSolutionStep()

#============================================================================================================================
    def InitializeSolutionStep(self):

        self.KratosPrintInfo("[STEP: " + str(self.step) + "  --  TIME: " + str(self.time) +  "  --  TIME_STEP: " + str(self.delta_time) + "]")

        # processes to be executed at the begining of the solution step
        self.model_processes.ExecuteInitializeSolutionStep()
        self.GraphicalOutputExecuteInitializeSolutionStep()
        self.solver.InitializeSolutionStep()

#============================================================================================================================
    def SolveSolutionStep(self):
        self.clock_time = self.StartTimeMeasuring()
        self.solver.Solve()
        self.StopTimeMeasuring(self.clock_time,"Solving", False)

#============================================================================================================================
    def FinalizeSolutionStep(self):

        self.GraphicalOutputExecuteFinalizeSolutionStep()            

        # processes to be executed at the end of the solution step
        self.model_processes.ExecuteFinalizeSolutionStep()

        # processes to be executed before witting the output      
        self.model_processes.ExecuteBeforeOutputStep()

        # write output results GiD: (frequency writing is controlled internally)
        self.GraphicalOutputPrintOutput()            

        # processes to be executed after witting the output
        self.model_processes.ExecuteAfterOutputStep()

#============================================================================================================================
    def Finalize(self):
        
        # Ending the problem (time integration finished)
        self.GraphicalOutputExecuteFinalize()        
        self.model_processes.ExecuteFinalize()
        self.KratosPrintInfo(" ")
        self.KratosPrintInfo("=================================================")
        self.KratosPrintInfo(" - Kratos FemDem Application Calculation End   - ")
        self.KratosPrintInfo("=================================================")
        self.KratosPrintInfo(" ")
        #### END SOLUTION ####
        # Measure process time
        tfp = timer.clock()
        # Measure wall time
        tfw = timer.time()
        KratosMultiphysics.Logger.PrintInfo("::[KSM Simulation]:: [Elapsed Time = %.2f" % (tfw - self.t0w),"seconds] (%.2f" % (tfp - self.t0p),"seconds of cpu/s time)")
        KratosMultiphysics.Logger.PrintInfo(timer.ctime())

 #============================================================================================================================       
    def SetGraphicalOutput(self):
        self.output_settings = self.ProjectParameters["output_configuration"]
        self.graphical_output = gid_output_process.GiDOutputProcess(self.computing_model_part,
                                      self.problem_name,
                                      self.output_settings)        
    #============================================================================================================================
    def GraphicalOutputExecuteInitialize(self):
        self.graphical_output.ExecuteInitialize() 
    #============================================================================================================================        
    def GraphicalOutputExecuteBeforeSolutionLoop(self):
        # writing a initial state results file or single file (if no restart)
        if((self.main_model_part.ProcessInfo).Has(KratosMultiphysics.IS_RESTARTED)):
            if(self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] == False):
                self.graphical_output.ExecuteBeforeSolutionLoop()
    #============================================================================================================================               
    def GraphicalOutputExecuteInitializeSolutionStep(self):
        self.graphical_output.ExecuteInitializeSolutionStep()
    #============================================================================================================================        
    def GraphicalOutputExecuteFinalizeSolutionStep(self):
        self.graphical_output.ExecuteFinalizeSolutionStep() 
    #============================================================================================================================        
    def GraphicalOutputPrintOutput(self):
        if(self.graphical_output.IsOutputStep()):
                self.graphical_output.PrintOutput()
    #============================================================================================================================
    def GraphicalOutputExecuteFinalize(self):
        self.graphical_output.ExecuteFinalize()
    #============================================================================================================================   
    def SetParallelSize(self, num_threads):
        parallel = KratosMultiphysics.OpenMPUtils()
        parallel.SetNumThreads(int(num_threads))
    #============================================================================================================================
    def GetParallelSize(self):
        parallel = KratosMultiphysics.OpenMPUtils()
        return parallel.GetNumThreads()    
    #============================================================================================================================    
    def StartTimeMeasuring(self):
        # Measure process time
        time_ip = timer.clock()
        return time_ip
    #============================================================================================================================
    def StopTimeMeasuring(self, time_ip, process, report):
        # Measure process time
        time_fp = timer.clock()
        if( report ):
            used_time = time_fp - time_ip
            print("::[KSM Simulation]:: [ %.2f" % round(used_time,2),"s", process," ] ")

    #============================================================================================================================

if __name__ == "__main__": 
    Solution().Run()

#============================================================================================================================

    def ComputeDeltaTime(self):

        if self.ProjectParameters["problem_data"].Has("time_step"):
            return self.ProjectParameters["problem_data"]["time_step"].GetDouble()

        elif self.ProjectParameters["problem_data"].Has("variable_time_steps"):

            current_time = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]
            for key in self.ProjectParameters["problem_data"]["variable_time_steps"].keys():
                interval_settings = self.ProjectParameters["problem_data"]["variable_time_steps"][key]
                interval = KratosMultiphysics.IntervalUtility(interval_settings)            
                # Getting the time step of the interval
                if interval.IsInInterval(current_time):
                    return interval_settings["time_step"].GetDouble()
                # If we arrive here we raise an error because the intervals are not well defined
                raise Exception("::[MechanicalSolver]:: Time stepping not well defined!")
        else:
            raise Exception("::[MechanicalSolver]:: Time stepping not defined!")
