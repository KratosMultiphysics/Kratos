from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import *
CheckForPreviousImport()

class SolvingInfoUtility(object):

    #
    def __init__(self, model_part, SolverSettings):

        self.model_part  = model_part
        self.time_integration_method = SolverSettings.time_integration_method
        self.convergence_criterion = "Residual_criteria"
        self.max_iters   = 30
        self.total_iters = 0

        if(hasattr(SolverSettings, "convergence_criterion")):
            self.convergence_criterion = SolverSettings.convergence_criterion
        if(hasattr(SolverSettings, "max_iteration")):
            self.max_iters = SolverSettings.max_iteration

        self.model_part.ProcessInfo[NL_ITERATION_NUMBER] = 0 
        self.step_iters=[]
        
        self.step_iters.append(0.5) #zero iters maximum at the begining
                    
        for i in range(0,2*self.max_iters): #solving can be segragated (2 iteration steps)
            self.step_iters.append(0)
            
        self.simulated_time  = 0
        self.max_iter_number = 0

        self.execute_write = True
        self.write_id      = 0

        self.execute_save  = False
        self.restart_id    = 0

        self.write_performed = False
        self.number_of_result_files = 0;

        self.restart_performed = False
        self.number_of_restart_files = 0;

        self.convergence_warning = False
        self.number_of_non_converged_steps = 0;
       
    #
    def set_convergence_criterion(self, convergence_criterion):
        self.convergence_criterion = convergence_criterion
 
    #
    def set_max_iters(self, max_iters):
        self.max_iters = max_iters
    
    #
    def print_step_info(self, current_time, current_step):
        self.model_part.ProcessInfo[TIME_STEPS] = current_step
        self.simulated_time = current_time
        print("  [ TIME=%.8e" % current_time,"]: STEP=",current_step," [-----------------] ")

    #
    def update_solving_info(self):
        if(self.time_integration_method != "Explicit"):
            num_iters = self.model_part.ProcessInfo[NL_ITERATION_NUMBER]
 
            if( num_iters >= self.max_iters ):
                self.number_of_non_converged_steps += 1
                self.convergence_warning = True
            else:
                self.step_iters[num_iters] += 1
                self.total_iters += num_iters
 
                if( self.max_iter_number < num_iters ):
                    self.max_iter_number = num_iters

    #
    def set_print_info(self, execute_write, write_id):
        self.execute_write = execute_write
        self.write_id      = write_id
        self.number_of_result_files += 1
        self.write_performed = True

    #
    def set_restart_info(self, execute_save, restart_id):
        self.execute_save = execute_save
        self.restart_id   = restart_id
        self.number_of_restart_files += 1
        self.restart_performed = True
    #
    def print_solving_info(self):

        num_iters         = self.model_part.ProcessInfo[NL_ITERATION_NUMBER]
        
        convergence_ratio = self.model_part.ProcessInfo[CONVERGENCE_RATIO]
        residual_norm     = self.model_part.ProcessInfo[RESIDUAL_NORM]

        #average_iters    = self.total_iters/current_step
        average_iters     = self.step_iters.index(max(self.step_iters))

        if(self.time_integration_method != "Explicit"):
            if( num_iters < self.max_iters ):
                total_info = "  [ ITERS="+str(num_iters)
            else:
                total_info = "  [ NOT CONVERGED|ITERS="+str(self.max_iters)
 
            total_info = total_info+"|AVG-ITERS="+str(average_iters)
            total_info = total_info+"|RATIO=%.4e" % convergence_ratio
            total_info = total_info+"|NORM=%.4e" % residual_norm
            total_info = total_info+" ]"
            if(self.execute_write):
                total_info = total_info + "[ WRITE_ID="+str(self.write_id)+" ]"
            if(self.execute_save):
                total_info = total_info + "[ RESTART_ID="+str(self.restart_id)+" ]" 
            print(total_info)   
        else:
            total_info = "  [ EXPLICIT_STABLE=yes"
            #total_info = total_info+"|NORM=%.4e" % residual_norm
            total_info = total_info+" ]"
            if(self.execute_write):
                total_info = total_info + "[ WRITE_ID="+str(self.write_id)+" ]"
            if(self.execute_save):
                total_info = total_info + "[ RESTART_ID="+str(self.restart_id)+" ]" 
            print(total_info)   
          
        self.execute_write = False
        self.execute_save  = False

    #
    def info_check(self):
            
        current_step = self.model_part.ProcessInfo[TIME_STEPS] 
        total_info = "::[Solving Info]:: [ SIMULATED-TIME=%.4e" % self.simulated_time
        total_info = total_info+" seconds | STEPS-NUMBER="+str(current_step)
        total_info = total_info+" ]"

        print(total_info)

        average_iters = self.step_iters.index(max(self.step_iters))
        
        total_info = "::[Solving Info]:: [ AVG-ITERS:"+str(average_iters)
        total_info = total_info+"|MAX-ITERS="+str(self.max_iter_number) 
        
        if( self.write_performed ):
            total_info = total_info+"|RESULTS-FILES="+str(self.number_of_result_files)
        else:
            total_info = total_info+"|WARNING: NO RESULTS WRITTEN"
            
        if( self.restart_performed ):
            total_info = total_info+"|RESTART-FILES="+str(self.number_of_restart_files)                

        total_info = total_info+" ]"

        print(total_info)

        if(self.convergence_warning):
            total_info = "::[Solving Info]:: [ NON-CONVERGED-STEPS="+str(self.number_of_non_converged_steps)
            total_info = total_info+" (WARNING::RESULTS WILL NOT BE ACCURATE)"
            total_info = total_info+" ]"
            print(total_info)

        if(self.write_performed):
            print( "::[Solving Info]:: -NOTE: run_GID_for_viewing_results-")
                    
        print("")
