# ==============================================================================
#  TopologyOptimizationApplication
#
#  License:         BSD License
#                   license: TopologyOptimizationApplication/license.txt
#
#  Main authors:    Baumgärtner Daniel, https://github.com/dbaumgaertner
#                   Octaviano Malfavón Farías
#
# ==============================================================================

# ------------------------------------------------------------------------------
# Imports
# ------------------------------------------------------------------------------

# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.TopologyOptimizationApplication import *

# For GID output
from gid_output import GiDOutput

# Further necessary imports
import csv
import math
import time

# ==============================================================================
def ConstructOptimizer( opt_model_part, config, analyzer ):

    # Creat optimizer according to selected optimization method
    if( config.optimization_method == "simp_method" ):
        optimizer = SIMPMethod( opt_model_part, config, analyzer )
        return optimizer

    else:
        sys.exit( "Specified optimization_method not implemented" )

# ==============================================================================
class SIMPMethod:

    # --------------------------------------------------------------------------
    def __init__(self, opt_model_part, config, analyzer):

        # Set Topology Optimization configurations
        self.config = config

        # For GID output
        self.gid_io = GiDOutput(config.GiD_output_file_name,
                        config.VolumeOutput,
                        config.GiDPostMode,
                        config.GiDMultiFileFlag,
                        config.GiDWriteMeshFlag,
                        config.GiDWriteConditionsFlag)

        # Set analyzer
        self.analyzer = analyzer

        # Set response functions
        self.objectives = config.objectives
        self.constraints = config.constraints

        print("\n::[Initializing Topology Optimization Application]::")
        print("  The following objectives are defined:")
        for func_id in config.objectives:
            print("   ",func_id,":",config.objectives[func_id],"\n")

        print("  The following constraints are defined:")
        for func_id in config.constraints:
            print("   ",func_id,":",config.constraints[func_id])

        # Create controller object
        self.controller = Controller( config );

        # Model parameters
        self.opt_model_part = opt_model_part

        # Initialize element variables
        for element_i in opt_model_part.Elements:
            element_i.SetValue(E_MIN, config.E_min)
            element_i.SetValue(PENAL, config.penalty)
            element_i.SetValue(X_PHYS, config.initial_volume_fraction)
            element_i.SetValue(X_PHYS_OLD, config.initial_volume_fraction)
            element_i.SetValue(E_0, opt_model_part.Properties[config.simp_property].GetValue(YOUNG_MODULUS))
            element_i.SetValue(YOUNG_MODULUS, opt_model_part.Properties[config.simp_property].GetValue(YOUNG_MODULUS))

        # Only happens if continuation strategy is activated (Initialization of penalty factor)
        if(self.config.continuation_strategy == 1):
            for element_i in self.opt_model_part.Elements:
                element_i.SetValue(PENAL,1)

        # Add toolbox for topology filtering utilities
        self.filter_utils = TopologyFilteringUtilities( opt_model_part,
                                                        self.config.filter_radius,
                                                        self.config.max_elements_in_filter_radius )

        # Add toolbox for topology updating utilities
        self.design_update_utils = TopologyUpdatingUtilities( opt_model_part )

        # Add toolbox for I/O
        self.io_utils = IOUtilities()

    # --------------------------------------------------------------------------
    def optimize(self):

        print("\n> ==============================================================================================================")
        print("> Starting optimization using the following algorithm: ",self.config.optimization_algorithm)
        print("> ==============================================================================================================")

        # Start timer and assign to object such that total time of opt may be measured at each step
        self.opt_start_time = time.time()

        # Initialize the design output in GiD format and print initial 0 state
        self.gid_io.initialize_results(self.opt_model_part)
        self.gid_io.write_results(0, self.opt_model_part, self.config.nodal_results, self.config.gauss_points_results)

        # Call for the specified optimization algorithm
        if(self.config.optimization_algorithm == "oc_algorithm"):
           self.start_oc_algorithm()

        else:
            raise TypeError("Specified optimization_algorithm not implemented!")

        # Finalize the design output in GiD format
        self.gid_io.finalize_results()

        # Stop timer
        opt_end_time = time.time()

        print("\n> ==============================================================================================================")
        print("> Finished optimization in ",round(opt_end_time - self.opt_start_time,1)," s!")
        print("> ==============================================================================================================")

    # --------------------------------------------------------------------------
    def start_oc_algorithm(self):

        # Get Id of objective & constraint
        only_F_id = None
        only_C_id = None
        for F_id in self.objectives:
            only_F_id = F_id
            break
        for C_id in self.constraints:
            only_C_id = C_id
            break

        # Initialize variables for comparison purposes in Topology Optimization Tool
        pmax                          = self.config.penalty   # Maximum penalty value used for continuation strategy
        Obj_Function                  = None
        Obj_Function_old              = None
        Obj_Function_initial          = None
        Obj_Function_relative_change  = None
        Obj_Function_absolute_change  = None

        # Print the Topology Optimization Settings that will be used in the program
        print("\n::[Topology Optimization Settings]::")
        print("  E_min:          ", self.config.E_min)
        print("  Filter radius:  ", self.config.filter_radius)
        print("  Penalty factor: ", self.config.penalty)
        print("  Rel. Tolerance: ", self.config.relative_tolerance)
        print("  Volume Fraction:", self.config.initial_volume_fraction)
        print("  Max. number of iterations:", self.config.max_opt_iterations)

        if (self.config.restart_write_frequency < self.config.max_opt_iterations):
            if (self.config.restart_write_frequency == 1):
                print("  Make a restart file every iteration")
            elif (self.config.restart_write_frequency > 1):
                print("  Make a restart file every", self.config.restart_write_frequency, "iterations")
            else:
                print("  No restart file will be done during the simulation")
        else:
            print("  No restart file will be done during the simulation")

        # Start optimization loop
        for opt_itr in range(1,self.config.max_opt_iterations+1):

            # Some output
            print("\n> ==============================================================================================")
            print("> Starting optimization iteration ",opt_itr)
            print("> ==============================================================================================\n")

            # Start measuring time needed for current optimization step
            start_time = time.time()

            # Initialize response container
            response = self.controller.create_response_container()

            # Set controller to evaluate objective & constraint
            self.controller.initialize_controls()
            self.controller.get_controls()[only_F_id]["calc_func"] = 1
            self.controller.get_controls()[only_C_id]["calc_func"] = 1

            # Set to evaluate objective & constraint gradient if provided
            if(self.objectives[only_F_id]["grad"]=="provided"):
                self.controller.get_controls()[only_F_id]["calc_grad"] = 1
            if(self.constraints[only_C_id]["grad"]=="provided"):
                self.controller.get_controls()[only_C_id]["calc_grad"] = 1

            # RUN FEM: Call analyzer with current X to compute response (global_strain_energy, dcdx)
            self.analyzer(self.controller.get_controls(), response, opt_itr)

            # Filter sensitivities
            print("\n::[Filter Sensitivities]::")
            self.filter_utils.ApplyFilter(self.config.filter_type , self.config.filter_kernel )

            # Update design variables ( densities )  --> new X by:
            print("\n::[Update Densities]::")
            self.design_update_utils.UpdateDensitiesUsingOCMethod( self.config.optimization_algorithm,
                                                                   self.config.initial_volume_fraction,
                                                                   self.config.grey_scale_filter,
                                                                   opt_itr,
                                                                   self.config.q_max )

            # Print of results
            print("\n::[RESULTS]::")
            Obj_Function = response[only_F_id]["func"]
            C_Function = response[only_C_id]["func"]

            print("  Obj. function value           = ", math.ceil(Obj_Function*1000000)/1000000 )
            print("  Const. function value         = ", math.ceil(C_Function*1000000)/1000000 )

            if opt_itr == 1:
                Obj_Function_initial = Obj_Function

            if opt_itr > 1:
                Obj_Function_relative_change = (Obj_Function - Obj_Function_old) / Obj_Function_initial
                print("  Relative Obj. Function change =", math.ceil((Obj_Function_relative_change*100)*10000)/10000, "%" )

                Obj_Function_absolute_change = (Obj_Function - Obj_Function_initial) / Obj_Function_initial
                print("  Absolute Obj. Function change =", math.ceil((Obj_Function_absolute_change*100)*10000)/10000, "%" )

            Obj_Function_old = Obj_Function

            # Write design in GiD format
            self.gid_io.write_results(opt_itr, self.opt_model_part, self.config.nodal_results, self.config.gauss_points_results)

            # Continuation Strategy
            if(self.config.continuation_strategy == 1):
                print("  Continuation Strategy for current iteration was ACTIVE")
                if opt_itr < 20:
                    for element_i in self.opt_model_part.Elements:
                        element_i.SetValue(PENAL, 1)
                else:
                    for element_i in self.opt_model_part.Elements:
                        element_i.SetValue(PENAL, min(pmax,1.02*element_i.GetValue(PENAL)))
            else:
                print("  Continuation Strategy for current iteration was UNACTIVE")

            # Write restart file every selected number of iterations
            restart_filename = self.config.restart_output_file.replace(".mdpa","_"+str(opt_itr)+".mdpa")
            if (self.config.restart_write_frequency > 0):
                if (opt_itr % self.config.restart_write_frequency == False):
                    print("\n::[Restart File]::")
                    print("  Saving file at iteration", opt_itr)
                    self.io_utils.SaveOptimizationResults(self.config.restart_input_file, self.opt_model_part, restart_filename)

            # Check convergence
            if opt_itr > 1:
                # Check if maximum iterations were reached
                if(opt_itr==self.config.max_opt_iterations):
                    end_time = time.time()
                    print("\n  Time needed for current optimization step = ",round(end_time - start_time,1),"s")
                    print("  Time needed for total optimization so far = ",round(end_time - self.opt_start_time,1),"s")
                    print("\n  Maximal iterations of optimization problem reached!")
                    self.io_utils.SaveOptimizationResults(self.config.restart_input_file, self.opt_model_part, restart_filename)
                    break

                # Check for relative tolerance
                if(abs(Obj_Function_relative_change)<self.config.relative_tolerance):
                    end_time = time.time()
                    print("\n  Time needed for current optimization step = ",round(end_time - start_time,1),"s")
                    print("  Time needed for total optimization so far = ",round(end_time - self.opt_start_time,1),"s")
                    print("\n  Optimization problem converged within a relative objective tolerance of",self.config.relative_tolerance)
                    self.io_utils.SaveOptimizationResults(self.config.restart_input_file, self.opt_model_part, restart_filename)
                    break

            # Set X_PHYS_OLD to X_PHYS to update the value for the next simulation's "change percentage"
            for element_i in self.opt_model_part.Elements:
                element_i.SetValue(X_PHYS_OLD, element_i.GetValue(X_PHYS))

            # Take time needed for current optimization step
            end_time = time.time()
            print("\n  Time needed for current optimization step = ",round(end_time - start_time,1),"s")
            print("  Time needed for total optimization so far = ",round(end_time - self.opt_start_time,1),"s")

# ==============================================================================
class Controller:

    # --------------------------------------------------------------------------
    def __init__( self, config ):

        # Create and initialize controller
        self.controls = {}
        for func_id in config.objectives:
            self.controls[func_id] = {"calc_func": 0, "calc_grad": 0}
        for func_id in config.constraints:
            self.controls[func_id] = {"calc_func": 0, "calc_grad": 0}

        # Initialize response container to provide storage for any response
        self.response_container = {}
        for func_id in config.objectives:
            self.response_container[func_id] = {"func": None, "grad": None}
        for func_id in config.constraints:
            self.response_container[func_id] = {"func": None, "grad": None}

    # --------------------------------------------------------------------------
    def initialize_controls( self ):

        # Sets
        for func_id in self.controls:
            self.controls[func_id] = {"calc_func": 0, "calc_grad": 0}

    # --------------------------------------------------------------------------
    def get_controls( self ):

        return self.controls

    # --------------------------------------------------------------------------
    def create_response_container( self ):

        # Create and initialize container to store any response defined
        for func_id in self.response_container:
            self.response_container[func_id] = {"func": None, "grad": None}

        # Return container
        return self.response_container

# ==============================================================================