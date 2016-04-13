# ==============================================================================
'''
 KratosShapeOptimizationApplication
 A library based on:
 Kratos
 A General Purpose Software for Multi-Physics Finite Element Analysis
 (Released on march 05, 2007).

 Copyright (c) 2016: Daniel Baumgaertner
                     daniel.baumgaertner@tum.de
                     Chair of Structural Analysis
                     Technische Universitaet Muenchen
                     Arcisstrasse 21 80333 Munich, Germany

 Permission is hereby granted, free  of charge, to any person obtaining
 a  copy  of this  software  and  associated  documentation files  (the
 "Software"), to  deal in  the Software without  restriction, including
 without limitation  the rights to  use, copy, modify,  merge, publish,
 distribute,  sublicense and/or  sell copies  of the  Software,  and to
 permit persons to whom the Software  is furnished to do so, subject to
 the following condition:

 Distribution of this code for  any  commercial purpose  is permissible
 ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

 The  above  copyright  notice  and  this permission  notice  shall  be
 included in all copies or substantial portions of the Software.

 THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
 EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
 CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
 TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
 SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
'''
#==============================================================================
#
#   Project Name:        KratosShape                            $
#   Created by:          $Author:    daniel.baumgaertner@tum.de $
#   Last modified by:    $Co-Author: daniel.baumgaertner@tum.de $
#   Date:                $Date:                      March 2016 $
#   Revision:            $Revision:                         0.0 $
#
# ==============================================================================

# ------------------------------------------------------------------------------
# Imports
# ------------------------------------------------------------------------------
# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division 

# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.ShapeOptimizationApplication import *

# check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

# For GID output
from gid_output import GiDOutput

# Further necessary imports
import csv
import math
import time

# ==============================================================================
class VertexMorphingMethod:

    # --------------------------------------------------------------------------
    def __init__(self,opt_model_part,config,analyzer):

        # For GID output
        self.gid_io = GiDOutput(config.design_surface_name,
                                config.VolumeOutput,
                                config.GiDPostMode,
                                config.GiDMultiFileFlag,
                                config.GiDWriteMeshFlag,
                                config.GiDWriteConditionsFlag)


        # Add variables needed for shape optimization
        opt_model_part.AddNodalSolutionStepVariable(NORMAL)
        opt_model_part.AddNodalSolutionStepVariable(NORMALIZED_SURFACE_NORMAL)
        opt_model_part.AddNodalSolutionStepVariable(OBJECTIVE_SENSITIVITY)
        opt_model_part.AddNodalSolutionStepVariable(CONSTRAINT_SENSITIVITY) 
        opt_model_part.AddNodalSolutionStepVariable(SHAPE_UPDATE) 
        opt_model_part.AddNodalSolutionStepVariable(SHAPE_CHANGE_ABSOLUTE)
        opt_model_part.AddNodalSolutionStepVariable(IS_ON_BOUNDARY)
        opt_model_part.AddNodalSolutionStepVariable(BOUNDARY_PLANE) 
        opt_model_part.AddNodalSolutionStepVariable(SHAPE_UPDATES_DEACTIVATED) 
        opt_model_part.AddNodalSolutionStepVariable(SENSITIVITIES_DEACTIVATED) 

        
        # Read and print model part (design surface)
        buffer_size = 1
        model_part_io = ModelPartIO(config.design_surface_name)
        model_part_io.ReadModelPart(opt_model_part)
        opt_model_part.SetBufferSize(buffer_size)
        print("\nThe following design surface was defined:\n\n",opt_model_part)

        # Set configurations
        self.config = config

        # Set analyzer
        self.analyzer = analyzer

        # Set response functions 
        self.objectives = config.objectives
        self.constraints = config.constraints

        print("> The following objectives are defined:\n")
        for func_id in config.objectives:
            print(func_id,":",config.objectives[func_id],"\n")

        print("> The following constraints are defined:\n")
        for func_id in config.constraints:
            print(func_id,":",config.constraints[func_id],"\n")

        # Create controller object
        self.controller = Controller( config );  

        # Model parameters 
        self.opt_model_part = opt_model_part

        # Add toolbox for vertex morphing    
        max_nodes_affected = 10000 # Specification required by spatial (tree) search 
                                   # Defines maximum nodes that may be considered within the filter radius
        self.vm_utils = VertexMorphingUtilities( self.opt_model_part,
                                                 self.config.domain_size,
                                                 self.objectives,
                                                 self.constraints,
                                                 config.filter_size,
                                                 max_nodes_affected)

    # --------------------------------------------------------------------------
    def optimize(self):

        print("\n> ==============================================================================================================")
        print("> Starting optimization using the following algorithm: ",self.config.optimization_algorithm)
        print("> ==============================================================================================================\n")

        # Start timer and assign to object such that total time of opt may be measured at each step
        self.opt_start_time = time.time()

        # Initialize design output in GID Format
        self.gid_io.initialize_results(self.opt_model_part)

        # Call for for specified optimization algorithm
        if(self.config.optimization_algorithm == "steepest_descent"):
           self.start_steepest_descent()

        elif(self.config.optimization_algorithm == "augmented_lagrange"):
           self.start_augmented_lagrange()

        elif(self.config.optimization_algorithm == "penalized_projection"):
           self.start_penalized_projection()           

        else:
            sys.exit("Specified optimization_algorithm not implemented!")

        # Finalize design output in GID formad
        self.gid_io.finalize_results()

        # Stop timer
        opt_end_time = time.time()

        print("\n> ==============================================================================================================")
        print("> Finished optimization in ",round(opt_end_time - self.opt_start_time,1)," s!")
        print("> ==============================================================================================================\n")

    # --------------------------------------------------------------------------
    def start_steepest_descent(self):

        # Flags to trigger proper function calls
        constraints_given = False

        # Get Id of objective
        only_F_id = None
        for F_id in self.objectives:
            only_F_id = F_id
            break

        # Initialize file where design evolution is recorded
        with open(self.config.design_history_file, 'w') as csvfile:
            historyWriter = csv.writer(csvfile, delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL)
            row = []
            row.append("itr")
            row.append("\tf")
            row.append("\t\N{GREEK CAPITAL LETTER DELTA}f_absolute[%]")
            row.append("\t\N{GREEK CAPITAL LETTER DELTA}f_relative[%]")
            historyWriter.writerow(row)    
        
        # Define initial design (initial design corresponds to a zero shape update)
        # Note that we assume an incremental design variable in vertex morphing
        X = {}
        for node in self.opt_model_part.Nodes:
            X[node.Id] = [0.,0.,0.]

        # Miscellaneous working variables for data management
        initial_f = 0.0
        previous_f = 0.0

        # Start optimization loop
        for opt_itr in range(1,self.config.max_opt_iterations+1):

            # Some output
            print("\n>===============================================")
            print("> Starting optimization iteration ",opt_itr)
            print(">===============================================\n")

            # Start measuring time needed for current optimization step
            start_time = time.time()

            # Set controller to evaluate objective
            self.controller.initialize_controls()
            self.controller.get_controls()[only_F_id]["calc_func"] = 1

            # Set to evaluate objective gradient if provided
            if(self.objectives[only_F_id]["grad"]=="provided"):
                self.controller.get_controls()[only_F_id]["calc_grad"] = 1

            # Initialize response container
            response = self.controller.create_response_container()          

            # Start analyzer according to specified controls
            iterator = str(opt_itr) + ".0"
            self.analyzer( X, self.controller.get_controls(), iterator, response )

            # Store gradients on the nodes of the model_part
            self.store_grads_on_nodes( response[only_F_id]["grad"] )

            # Compute unit surface normals at each node of current design
            self.vm_utils.compute_unit_surface_normals()

            # Project received gradients on unit surface normal at each node to obtain normal gradients
            self.vm_utils.project_grad_on_unit_surface_normal( constraints_given )

            # Compute filtered gradients (do backward mapping)
            self.vm_utils.filter_gradients( constraints_given )

            # Compute search direction
            self.vm_utils.compute_search_direction_steepest_descent()

            # Compute design variable update (do one step in the optimization algorithm)
            self.vm_utils.update_design_variable( self.config.step_size_0 )

            # Compute shape update (do forward mapping)
            self.vm_utils.update_shape()

            # Compute and output some measures to track changes in the objective function
            delta_f_absolute = 0.0
            delta_f_relative = 0.0
            print("\n> Current value of objective function = ",response[only_F_id]["func"])
            if(opt_itr>1):
                delta_f_absolute = 100* ( response[only_F_id]["func"]/initial_f - 1 )
                delta_f_relative = 100*( response[only_F_id]["func"]/previous_f - 1 )
                print("\n> Absolut change of objective function = ",round(delta_f_absolute,6)," [%]")
                print("\n> Relative change of objective function = ",round(delta_f_relative,6)," [%]")            

            # Write design history to file
            with open(self.config.design_history_file, 'a') as csvfile:
                historyWriter = csv.writer(csvfile, delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL)
                row = []
                row.append(str(opt_itr)+"\t")
                row.append("\t"+str("%.12f"%(response[only_F_id]["func"]))+"\t")
                row.append("\t"+str("%.2f"%(delta_f_absolute))+"\t")
                row.append("\t"+str("%.6f"%(delta_f_relative))+"\t")
                historyWriter.writerow(row)

            # Write design in GID format
            self.gid_io.write_results(opt_itr, self.opt_model_part, self.config.nodal_results, [])                

            # Check convergence
            if(opt_itr>1):

                # Check if maximum iterations were reached
                if(opt_itr==self.config.max_opt_iterations):
                    print("\n> Maximal iterations of supoptimization problem reached!")
                    break

                # Check for relative tolerance
                if(abs(delta_f_relative)<self.config.relative_tolerance_objective):
                    print("\n> Optimization subproblem converged within a relative objective tolerance of ",self.config.relative_tolerance_objective,".")
                    break

                # # Check if value of objective increases
                # if(response[only_F_id]["func"]>previous_f):
                #     print("\n> Value of objective function increased!")
                #     break

            # Update design
            X = self.get_design()

            # Store values of for next iteration
            previous_f = response[only_F_id]["func"]

            # Store initial objective value
            if(opt_itr==1):
                initial_f = response[only_F_id]["func"]

            # Take time needed for current optimization step
            end_time = time.time()
            print("\n> Time needed for current optimization step = ",round(end_time - start_time,1),"s")
            print("\n> Time needed for total optimization so far = ",round(end_time - self.opt_start_time,1),"s")

    # --------------------------------------------------------------------------
    def start_augmented_lagrange(self):

        # README!!!
        # Note that the current implementation assumes that only one scalar objective & constraint is present

        # Flags to trigger proper function calls
        constraints_given = True

        # Get Id of objective
        only_F_id = None
        for F_id in self.objectives:
            only_F_id = F_id
            break

        # Get Id of constraint
        only_C_id = None
        for C_id in self.constraints:
            only_C_id = C_id
            break   

        # Initialize the optimization algorithm
        self.vm_utils.initialize_augmented_lagrange( self.constraints, 
                                                     self.config.penalty_fac_0, 
                                                     self.config.gamma, 
                                                     self.config.penalty_fac_max, 
                                                     self.config.lambda_0 )

        # Initialize file where design evolution is recorded
        with open(self.config.design_history_file, 'w') as csvfile:
            historyWriter = csv.writer(csvfile, delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL)
            row = []
            row.append("itr\t")
            row.append("\tsub_itr\t")
            row.append("\tl\t")
            row.append("\t\N{GREEK CAPITAL LETTER DELTA}l_relative[%]\t")
            row.append("\tf\t")
            row.append("\t\N{GREEK CAPITAL LETTER DELTA}f_absolute[%]\t")
            row.append("\tpenalty_fac\t")
            row.append("\tC["+str(only_C_id)+"]: "+str(self.constraints[only_C_id]["type"])+"\t") 
            row.append("\tlambda["+str(only_C_id)+"]\t")
            historyWriter.writerow(row)  

        # Define initial design (initial design corresponds to a zero shape update)
        # Note that we assume an incremental design variable in vertex morphing
        X = {}
        for node in self.opt_model_part.Nodes:
            X[node.Id] = [0.,0.,0.]   

        # Miscellaneous working variables for data management
        initial_f = 0.0
        previous_l = 0.0         

        # Start primary optimization loop
        for opt_itr in range(1,self.config.max_opt_iterations+1):

            # Some output
            print("\n>===============================================")
            print("> Starting optimization iteration ",opt_itr)
            print(">===============================================\n")

            # Start measuring time needed for current optimization step
            start_time = time.time()

            # Solve optimization subproblem
            for sub_opt_itr in range(1,self.config.max_sub_opt_iterations+1):

                # Some output
                print("\n>===============================================")
                print("> Starting suboptimization iteration ",sub_opt_itr)
                print(">===============================================\n")

                # Start measuring time needed for current suboptimization step
                subopt_start_time = time.time()

                # Set controller to evaluate objective and constraint
                self.controller.initialize_controls()
                self.controller.get_controls()[only_F_id]["calc_func"] = 1
                self.controller.get_controls()[only_C_id]["calc_func"] = 1

                # Set controller to evaluate objective and constraint gradient if provided
                if(self.objectives[only_F_id]["grad"]=="provided"):
                    self.controller.get_controls()[only_F_id]["calc_grad"] = 1
                if(self.constraints[only_C_id]["grad"]=="provided"):
                    self.controller.get_controls()[only_C_id]["calc_grad"] = 1                

                # Initialize response container
                response = self.controller.create_response_container()          

                # Start analyzer according to specified controls
                iterator = str(opt_itr) + "." + str(sub_opt_itr)
                self.analyzer( X, self.controller.get_controls(), iterator, response )

                # Evaluate Lagrange function
                l = self.vm_utils.get_value_of_augmented_lagrangian( only_F_id, self.constraints, response )

                # Store gradients on the nodes of the model_part
                self.store_grads_on_nodes( response[only_F_id]["grad"], response[only_C_id]["grad"] )

                # Compute unit surface normals at each node of current design
                self.vm_utils.compute_unit_surface_normals()

                # Project received gradients on unit surface normal at each node to obtain normal gradients
                dJdX_n = self.vm_utils.project_grad_on_unit_surface_normal( constraints_given )

                # Compute filtered gradients (do backward mapping)
                self.vm_utils.filter_gradients( constraints_given )

                # Compute search direction
                self.vm_utils.compute_search_direction_augmented_lagrange( self.constraints, response )

                # Compute design variable update (do one step in the optimization algorithm)
                self.vm_utils.update_design_variable( self.config.step_size_0 )

                # Compute shape update (do forward mapping)
                self.vm_utils.update_shape()

                # Compute and output some measures to track changes in the objective function
                delta_f_absolute = 0.0
                delta_l_relative = 0.0
                print("\n> Current value of Lagrange function = ",round(l,12))
                if(sub_opt_itr>1):
                    delta_f_absolute = 100* ( response[only_F_id]["func"]/initial_f - 1 )
                    delta_l_relative = 100*( l/previous_l - 1 )
                    print("\n> Relative change of Lagrange function = ",round(delta_l_relative,6)," [%]")

                # We write every major and every suboptimization iteration in design history
                with open(self.config.design_history_file, 'a') as csvfile:
                   historyWriter = csv.writer(csvfile, delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL)
                   row = []
                   row.append(str(opt_itr)+"\t")
                   row.append("\t"+str(sub_opt_itr)+"\t")
                   row.append("\t"+str("%.12f"%(l))+"\t")
                   row.append("\t"+str("%.6f"%(delta_l_relative))+"\t")
                   row.append("\t"+str("%.12f"%(response[only_F_id]["func"]))+"\t")
                   row.append("\t"+str("%.2f"%(delta_f_absolute))+"\t")
                   row.append("\t"+str("%.2f"%(self.vm_utils.get_penalty_fac()))+"\t")
                   row.append("\t"+str("%.12f"%(response[only_C_id]["func"]))+"\t")
                   row.append("\t"+str("%.12f"%(self.vm_utils.get_lambda(only_C_id)))+"\t")
                   historyWriter.writerow(row)

                # Write design in GID format
                write_itr = float(iterator)
                self.gid_io.write_results(write_itr, self.opt_model_part, self.config.nodal_results, [])

                # Check convergence (We ensure that at least 2 subiterations are done)
                if(sub_opt_itr>3):

                    # Check if maximum iterations were reached
                    if(sub_opt_itr==self.config.max_sub_opt_iterations): 
                        print("\n> Maximal iterations of supoptimization problem reached!")
                        break

                    # Check for relative tolerance
                    if(abs(delta_l_relative)<self.config.relative_tolerance_sub_opt): 
                        print("\n> Optimization subproblem converged within a relative objective tolerance of ",self.config.relative_tolerance_sub_opt,".")
                        break

                    # Check if value of lagrangian increases
                    if(l>previous_l):
                        print("\n> Value of Lagrange function increased!")
                        break

                # Update design
                X = self.get_design()   

                # Store value of Lagrange function for next iteration
                previous_l = l

                # Store initial objective value
                if(opt_itr==1 and sub_opt_itr==1):
                    initial_f = response[only_F_id]["func"]

                # Take time needed for current suboptimization step as well as for the overall opt so far
                subopt_end_time = time.time()
                print("\n> Time needed for current suboptimization step = ",round(subopt_end_time - subopt_start_time,1),"s")
                print("\n> Time needed for total optimization so far = ",round(subopt_end_time - self.opt_start_time,1),"s")

            # Check Convergence (More convergence criterion for major optimization iteration to be implemented!)

            # Check if maximum iterations were reached
            if(opt_itr==self.config.max_opt_iterations):
                print("\n> Maximal iterations of optimization problem reached!")
                break
            
            # Update lagrange multipliers and penalty factor
            self.vm_utils.udpate_augmented_lagrange_parameters( self.constraints, response )  

            # Take time needed for current optimization step
            end_time = time.time()
            print("\n> Time needed for current optimization step = ",round(end_time - start_time,1),"s")

    # --------------------------------------------------------------------------
    def start_penalized_projection(self):

        # README!!!
        # Note that the current implementation assumes that only one scalar objective & constraint is present

        # Flags to trigger proper function calls
        constraints_given = True

        # Get Id of objective
        only_F_id = None
        for F_id in self.objectives:
            only_F_id = F_id
            break

        # Get Id of constraint
        only_C_id = None
        for C_id in self.constraints:
            only_C_id = C_id
            break            

        # Initialize file where design evolution is recorded
        with open(self.config.design_history_file, 'w') as csvfile:
            historyWriter = csv.writer(csvfile, delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL)
            row = []
            row.append("itr")
            row.append("\tf")
            row.append("\t\N{GREEK CAPITAL LETTER DELTA}f_absolute[%]")
            row.append("\t\N{GREEK CAPITAL LETTER DELTA}f_relative[%]")
            row.append("\tC["+str(only_C_id)+"]: "+str(self.constraints[only_C_id]["type"])+"\t")           
            historyWriter.writerow(row)    
        
        # Define initial design (initial design corresponds to a zero shape update)
        # Note that we assume an incremental design variable in vertex morphing
        X = {}
        for node in self.opt_model_part.Nodes:
            X[node.Id] = [0.,0.,0.]

        # Miscellaneous working variables for data management
        initial_f = 0.0
        previous_f = 0.0            

        # Start optimization loop
        for opt_itr in range(1,self.config.max_opt_iterations+1):

            # Some output
            print("\n>===============================================")
            print("> Starting optimization iteration ",opt_itr)
            print(">===============================================\n")

            # Start measuring time needed for current optimization step
            start_time = time.time()

            # Set controller to evaluate objective and constraint
            self.controller.initialize_controls()
            self.controller.get_controls()[only_F_id]["calc_func"] = 1
            self.controller.get_controls()[only_C_id]["calc_func"] = 1

            # Set controller to evaluate objective and constraint gradient if provided
            if(self.objectives[only_F_id]["grad"]=="provided"):
                self.controller.get_controls()[only_F_id]["calc_grad"] = 1
            if(self.constraints[only_C_id]["grad"]=="provided"):
                self.controller.get_controls()[only_C_id]["calc_grad"] = 1                

            # Initialize response container
            response = self.controller.create_response_container()          

            # Start analyzer according to specified controls
            iterator = str(opt_itr) + ".0"
            self.analyzer( X, self.controller.get_controls(), iterator, response )

            # Check if constraint is active
            for func_id in self.config.constraints:
                if(self.config.constraints[func_id]["type"] == "eq"):
                    constraints_given = True
                elif(response[func_id]["func"]>0):
                    constraints_given = True
                else:
                    constraints_given = False
                break

            # Store gradients on the nodes of the model_part
            self.store_grads_on_nodes( response[only_F_id]["grad"], response[only_C_id]["grad"] )

            # Compute unit surface normals at each node of current design
            self.vm_utils.compute_unit_surface_normals()

            # Project received gradients on unit surface normal at each node to obtain normal gradients
            self.vm_utils.project_grad_on_unit_surface_normal( constraints_given )

            # Compute filtered gradients (do backward mapping)
            self.vm_utils.filter_gradients( constraints_given )

            # Compute search direction
            if(constraints_given):
                self.vm_utils.compute_search_direction_penalized_projection( response[only_C_id]["func"]*self.config.constraint_scaling)
            else:
                self.vm_utils.compute_search_direction_steepest_descent()

            # Compute design variable update (do one step in the optimization algorithm)
            self.vm_utils.update_design_variable( self.config.step_size_0 )

            # Compute shape update (do forward mapping)
            self.vm_utils.update_shape()

            # Compute and output some measures to track changes in the objective function
            delta_f_absolute = 0.0
            delta_f_relative = 0.0
            print("\n> Current value of objective function = ",response[only_F_id]["func"])
            if(opt_itr>1):
                delta_f_absolute = 100* ( response[only_F_id]["func"]/initial_f - 1 )
                delta_f_relative = 100*( response[only_F_id]["func"]/previous_f - 1 )
                print("\n> Absolut change of objective function = ",round(delta_f_absolute,6)," [%]")
                print("\n> Relative change of objective function = ",round(delta_f_relative,6)," [%]")           
                print("\n> Current value of constraint function = ",round(response[only_C_id]["func"],12))

            # Write design history to file
            with open(self.config.design_history_file, 'a') as csvfile:
                historyWriter = csv.writer(csvfile, delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL)
                row = []
                row.append(str(opt_itr)+"\t")
                row.append("\t"+str("%.12f"%(response[only_F_id]["func"]))+"\t")
                row.append("\t"+str("%.2f"%(delta_f_absolute))+"\t")
                row.append("\t"+str("%.6f"%(delta_f_relative))+"\t")
                row.append("\t"+str("%.12f"%(response[only_C_id]["func"]))+"\t")
                historyWriter.writerow(row)

            # Write design in GID format
            self.gid_io.write_results(opt_itr, self.opt_model_part, self.config.nodal_results, [])                

            # Check convergence (Further convergence criterions to be implemented )
            if(opt_itr>1):

                # Check if maximum iterations were reached
                if(opt_itr==self.config.max_opt_iterations):
                    print("\n> Maximal iterations of supoptimization problem reached!")
                    break

            # Update design
            X = self.get_design()

            # Store values of for next iteration
            previous_f = response[only_F_id]["func"]

            # Store initial objective value
            if(opt_itr==1):
                initial_f = response[only_F_id]["func"]            

            # Take time needed for current optimization step
            end_time = time.time()
            print("\n> Time needed for current optimization step = ",round(end_time - start_time,1),"s")
            print("\n> Time needed for total optimization so far = ",round(end_time - self.opt_start_time,1),"s")

    # --------------------------------------------------------------------------
    def store_grads_on_nodes(self,objective_grads,constraint_grads={}):

        # Read objective gradients
        eucledian_norm_obj_sens = 0.0
        for node_Id in objective_grads:

            # If deactivated, nodal sensitivities will not be assigned
            if(self.opt_model_part.Nodes[node_Id].GetSolutionStepValue(SENSITIVITIES_DEACTIVATED)):
                continue

            # If not deactivated, nodal sensitivities will be assigned
            sens_i = Vector(3)
            sens_i[0] = objective_grads[node_Id][0]
            sens_i[1] = objective_grads[node_Id][1]
            sens_i[2] = objective_grads[node_Id][2]           
            self.opt_model_part.Nodes[node_Id].SetSolutionStepValue(OBJECTIVE_SENSITIVITY,0,sens_i)
            eucledian_norm_obj_sens += sens_i[0] * sens_i[0] + sens_i[1] * sens_i[1] + sens_i[2] * sens_i[2]
        eucledian_norm_obj_sens = math.sqrt(eucledian_norm_obj_sens) 

        # When constraint_grads is defined also store constraint sensitivities (bool returns false if dictionary is empty)
        if(bool(constraint_grads)):
            eucledian_norm_cons_sens = 0.0
            for node_Id in constraint_grads:

                # If deactivated, nodal sensitivities will not be assigned
                if(self.opt_model_part.Nodes[node_Id].GetSolutionStepValue(SENSITIVITIES_DEACTIVATED)):
                    continue

                # If not deactivated, nodal sensitivities will be assigned
                sens_i = Vector(3)
                sens_i[0] = constraint_grads[node_Id][0]
                sens_i[1] = constraint_grads[node_Id][1]
                sens_i[2] = constraint_grads[node_Id][2]           
                self.opt_model_part.Nodes[node_Id].SetSolutionStepValue(CONSTRAINT_SENSITIVITY,0,sens_i)
                eucledian_norm_cons_sens += sens_i[0] * sens_i[0] + sens_i[1] * sens_i[1] + sens_i[2] * sens_i[2]
            eucledian_norm_cons_sens = math.sqrt(eucledian_norm_cons_sens)                  

    # --------------------------------------------------------------------------
    def get_design(self):

        # Read and return the current design in the corresponding mode
        X = {}

        if(self.config.design_output_mode=="relative"):
            for node in self.opt_model_part.Nodes:
                X[node.Id] = node.GetSolutionStepValue(SHAPE_UPDATE)

        elif(self.config.design_output_mode=="total"):
            for node in self.opt_model_part.Nodes:
                X[node.Id] = node.GetSolutionStepValue(SHAPE_CHANGE_ABSOLUTE)

        elif(self.config.design_output_mode=="absolute"):
            for node in self.opt_model_part.Nodes:
                X[node.Id] = [node.X,node.Y,node.Z]

        else:
            sys.exit("Wrong definition of design_output_mode!")

        return X

    # --------------------------------------------------------------------------

# ==============================================================================
class Controller:

    # --------------------------------------------------------------------------
    def __init__( self, config ):

        # Create and initialize controler
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

    # --------------------------------------------------------------------------

# ==============================================================================
def CreateOptimizer( opt_model_part, config, analyzer ):

    # Create folder where all functions includig the optimizer may store their design history in
    os.system( "rm -rf " + config.design_history_directory )
    os.system( "mkdir -p " + config.design_history_directory )

    # Creat optimizer according to selected optimization method
    if( config.design_control == "vertex_morphing" ):
        optimizer = VertexMorphingMethod( opt_model_part, config, analyzer )
        return optimizer

    else:
        sys.exit( "Specified design_control not implemented" )

# ==============================================================================
