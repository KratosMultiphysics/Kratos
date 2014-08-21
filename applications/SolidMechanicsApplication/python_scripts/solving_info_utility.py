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
        for i in range(0,self.max_iters+1):
            self.step_iters.append(0)


        self.execute_write = True
        self.write_id      = 0

        self.execute_save  = False
        self.restart_id    = 0

        self.write_perfomed = False
        self.convergence_warning = False
       
    #
    def set_convergence_criterion(self, convergence_criterion):
        self.convergence_criterion = convergence_criterion
 
    #
    def set_max_iters(self, max_iters):
        self.max_iters = max_iters
    
    #
    def print_step_info(self, current_time, current_step):
        self.model_part.ProcessInfo[TIME_STEPS] = current_step
        print("  [ TIME=%.8e" % current_time,"]: STEP=",current_step," [-----------------] ")

    #
    def update_solving_info(self):
        if(self.time_integration_method != "Explicit"):
            num_iters = self.model_part.ProcessInfo[NL_ITERATION_NUMBER]
            self.step_iters[num_iters] += 1
            self.total_iters += num_iters
 
            if( num_iters < self.max_iters ):
                convergence_warning = True

    #
    def set_print_info(self, execute_write, write_id):
        self.execute_write = execute_write
        self.write_id      = write_id
        self.write_performed = True

    #
    def set_restart_info(self, execute_save, restart_id):
        self.execute_save = execute_save
        self.restart_id   = restart_id

    #
    def print_solving_info(self):

        num_iters         = self.model_part.ProcessInfo[NL_ITERATION_NUMBER]
        
        convergence_ratio = self.model_part.ProcessInfo[CONVERGENCE_RATIO]
        residual_norm     = self.model_part.ProcessInfo[RESIDUAL_NORM]

        #average_iters     = self.total_iters/current_step
        average_iters     = self.step_iters.index(max(self.step_iters))

        if(self.time_integration_method != "Explicit"):
            if( num_iters < self.max_iters ):
                total_info = "  [ ITERS="+str(num_iters)
            else:
                total_info = "  [ NO CONVERGED|ITERS="+str(self.max_iters)
 
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
            print("  [ STABLE= yes|NORM=",residual_norm,"]")    
          
        self.execute_write = False
        self.execute_save  = False
    #
    def info_check(self):

        if(self.convergence_warning):
            print( "::[Solving Info]:: Some steps did not converge: RESULTS ARE NOT ACCURATE ")
        else:
            if(self.write_performed):
                print( "::[Solving Info]:: -NOTE: run_GID_for_viewing_results-")
            else:
                print( "::[Solving Info]:: -WARNING: no results written-")
