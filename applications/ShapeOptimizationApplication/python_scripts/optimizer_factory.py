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
#   Date:                $Date:                   December 2016 $
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
    def __init__( self, inputModelPart, optimizationSettings ):

        self.gidIO = GiDOutput( optimizationSettings.design_surface_submodel_part_name,
                                optimizationSettings.VolumeOutput,
                                optimizationSettings.GiDPostMode,
                                optimizationSettings.GiDMultiFileFlag,
                                optimizationSettings.GiDWriteMeshFlag,
                                optimizationSettings.GiDWriteConditionsFlag)

        self.inputModelPart = inputModelPart
        self.optimizationSettings = optimizationSettings
        self.objectives = optimizationSettings.objectives
        self.constraints = optimizationSettings.constraints

        self.addVariablesNeededForOptimization()

        self.outputInformationAboutResponseFunctions()

    # --------------------------------------------------------------------------
    def importModelPart( self ):
        model_part_io = ModelPartIO( self.optimizationSettings.inputModelPart_name )
        model_part_io.ReadModelPart( self.inputModelPart )
        buffer_size = 1
        self.inputModelPart.SetBufferSize( buffer_size )
        self.inputModelPart.ProcessInfo.SetValue(DOMAIN_SIZE,self.optimizationSettings.domain_size)

    # --------------------------------------------------------------------------
    def importFunctionForDesignAnalysis( self, importedFunctionForDesignAnalysis ): 
        self.analyzer = Analyzer( importedFunctionForDesignAnalysis, self.optimizationSettings ) 

    # --------------------------------------------------------------------------
    def optimize( self ):
        
        print("\n> ==============================================================================================================")
        print("> Starting optimization using the following algorithm: ",self.optimizationSettings.optimization_algorithm)
        print("> ==============================================================================================================\n")

        print(time.ctime())
        
        self.timeAtStartOfOptimization = time.time()

        self.prepareOptimization()

        iteratorForInitialDesign = 0
        self.gidIO.initialize_results(self.optimizationModelPart)
        self.gidIO.write_results(iteratorForInitialDesign, self.optimizationModelPart, self.optimizationSettings.nodal_results, [])   

        self.runOptimizationAlgorithmSpecifiedInTheOptimizationSettings()

        self.gidIO.finalize_results()

        timeAtEndOfOptimization = time.time()

        print("\n> ==============================================================================================================")
        print("> Finished optimization in ",round(timeAtEndOfOptimization - self.timeAtStartOfOptimization,2)," s!")
        print("> ==============================================================================================================\n")

    # --------------------------------------------------------------------------
    def prepareOptimization( self ):

        # Check if part to be optimized is defined as sub-model in the input model part and extract it as separate model part
        if( self.inputModelPart.HasSubModelPart( self.optimizationSettings.design_surface_submodel_part_name) ):
            self.optimizationModelPart = self.inputModelPart.GetSubModelPart( self.optimizationSettings.design_surface_submodel_part_name )
            print("\nThe following design surface was defined:\n\n",self.optimizationModelPart)
        else:
            raise RuntimeError("Sub-model part specified for optimization does not exist!")

        # Identify possible damping regions and if necessary list corresponding sub-model parts and settings
        damping_regions = []
        if(self.optimizationSettings.perform_damping):
            print("\n> The following damping regions are defined: \n")
            for region in self.optimizationSettings.damping_regions:
                sub_mdpa_name = region[0]
                damp_in_X = region[1]
                damp_in_Y = region[2]
                damp_in_Z = region[3]
                damping_function = region[4]
                damping_radius = region[5]
                if( self.inputModelPart.HasSubModelPart(sub_mdpa_name) ):
                    print(region)
                    damping_regions.append( [ self.inputModelPart.GetSubModelPart(sub_mdpa_name), 
                                            damp_in_X, 
                                            damp_in_Y, 
                                            damp_in_Z, 
                                            damping_function, 
                                            damping_radius] )
                else:
                    raise ValueError("The following sub-model part specified for damping does not exist: ",sub_mdpa_name)

        # Create communicator object
        self.communicator = Communicator(self.optimizationSettings);                     

        # Create mapper to map between geometry and design space 
        self.mapper = VertexMorphingMapper(self.optimizationModelPart,
                                           self.optimizationSettings.filter_function,
                                           self.optimizationSettings.use_mesh_preserving_filter_matrix,
                                           self.optimizationSettings.filter_size,
                                           self.optimizationSettings.perform_damping,
                                           damping_regions)
        
        # Toolbox to perform optimization
        self.opt_utils = OptimizationUtilities(self.optimizationModelPart,
                                               self.objectives,
                                               self.constraints,
                                               self.optimizationSettings.step_size,
                                               self.optimizationSettings.normalize_search_direction)

        # Toolbox to pre & post process geometry data
        self.geom_utils = GeometryUtilities(self.optimizationModelPart)

    # --------------------------------------------------------------------------
    def runOptimizationAlgorithmSpecifiedInTheOptimizationSettings( self ):
        if(self.optimizationSettings.optimization_algorithm == "steepest_descent"):
           self.runSteepestDescentAlgorithm()
        elif(self.optimizationSettings.optimization_algorithm == "augmented_lagrange"):
           self.runAugmentedLagrangeAlgorithm()
        elif(self.optimizationSettings.optimization_algorithm == "penalized_projection"):
           self.runPenalizedProjectionAlgorithm()           
        else:
            sys.exit("Specified optimization_algorithm not implemented!")

    # --------------------------------------------------------------------------
    def runSteepestDescentAlgorithm( self ):

        # Flags to trigger proper function calls
        constraints_given = False

        # Get Id of objective
        only_F_id = None
        for F_id in self.objectives:
            only_F_id = F_id
            break

        # Initialize file where design evolution is recorded
        with open(self.optimizationSettings.design_history_file, 'w') as csvfile:
            historyWriter = csv.writer(csvfile, delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL)
            row = []
            row.append("itr\t")
            row.append("\tf\t")
            row.append("\tdf_absolute[%]\t")
            row.append("\tdf_relative[%]\t")
            row.append("\tstep_size[-]\t")
            row.append("\tt_iteration[s]\t")
            row.append("\tt_total[s]") 
            historyWriter.writerow(row)    
        
        # Define initial design (initial design corresponds to a zero shape update)
        # Note that we assume an incremental design variable in vertex morphing
        X = {}
        for node in self.optimizationModelPart.Nodes:
            X[node.Id] = [0.,0.,0.]

        # Miscellaneous working variables for data management
        initial_f = 0.0
        previous_f = 0.0

        # Start optimization loop
        for opt_itr in range(1,self.optimizationSettings.max_opt_iterations+1):

            # Some output
            print("\n>===================================================================")
            print("> ",time.ctime(),": Starting optimization iteration ",opt_itr)
            print(">===================================================================\n")

            timeAtStartOfCurrentOptimizationStep = time.time()

            self.communicator.deleteAllRequests()
            self.communicator.requestFunctionValueOf(only_F_id)
            self.communicator.requestGradientOf(only_F_id)   

            iterator = str(opt_itr) + ".0"
            self.analyzer(X, self.communicator.getControls(), iterator, response)

            self.storeGradientsOnNodes( self.communicator.getReportedGradientOf[only_F_id] )

            self.geom_utils.compute_unit_surface_normals()

            # Project received gradients on unit surface normal at each node to obtain normal gradients
            self.geom_utils.project_grad_on_unit_surface_normal( constraints_given )

            # Compute mapping matrix
            self.mapper.compute_mapping_matrix()

            # Map sensitivities to design space
            self.mapper.map_sensitivities_to_design_space( constraints_given )

            # Compute search direction
            self.opt_utils.compute_search_direction_steepest_descent()

            # # Adjustment of step size
            # if( opt_itr > 1 and response[only_F_id]["value"]>previous_f):
            #     self.optimizationSettings.step_size = self.optimizationSettings.step_size/2

            # Compute design variable update (do one step in the optimization algorithm)
            self.opt_utils.compute_design_update()

            # Map design update to geometry space
            self.mapper.map_design_update_to_geometry_space()

            # Compute and output some measures to track changes in the objective function
            delta_f_absolute = 0.0
            delta_f_relative = 0.0
            print("\n> Current value of objective function = ",response[only_F_id]["value"])
            if(opt_itr>1):
                delta_f_absolute = 100*(response[only_F_id]["value"]-initial_f)/initial_f
                delta_f_relative = 100*(response[only_F_id]["value"]-previous_f)/initial_f
                print("> Absolut change of objective function = ",round(delta_f_absolute,6)," [%]")
                print("> Relative change of objective function = ",round(delta_f_relative,6)," [%]")            

            # Take time needed for current optimization step
            timeAtEndOfCurrentOptimizationStep = time.time()
            time_current_step = round(timeAtEndOfCurrentOptimizationStep - timeAtStartOfCurrentOptimizationStep,2)
            time_optimization = round(timeAtEndOfCurrentOptimizationStep - self.timeAtStartOfOptimization,2)
            print("\n> Time needed for current optimization step = ",time_current_step,"s")
            print("> Time needed for total optimization so far = ",time_optimization,"s")     
            
            # Write design in GID format
            self.gidIO.write_results(opt_itr, self.optimizationModelPart, self.optimizationSettings.nodal_results, [])   

            # Write design history to file
            with open(self.optimizationSettings.design_history_directory+"/"+self.optimizationSettings.design_history_file, 'a') as csvfile:
                historyWriter = csv.writer(csvfile, delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL)
                row = []
                row.append(str(opt_itr)+"\t")
                row.append("\t"+str("%.12f"%(response[only_F_id]["value"]))+"\t")
                row.append("\t"+str("%.2f"%(delta_f_absolute))+"\t")
                row.append("\t"+str("%.6f"%(delta_f_relative))+"\t")
                row.append("\t"+str(self.optimizationSettings.step_size)+"\t")
                row.append("\t"+str("%.1f"%(time_current_step))+"\t")
                row.append("\t"+str("%.1f"%(time_optimization)))
                historyWriter.writerow(row)              

            # Check convergence
            if(opt_itr>1):

                # Check if maximum iterations were reached
                if(opt_itr==self.optimizationSettings.max_opt_iterations):
                    print("\n> Maximal iterations of optimization problem reached!")
                    break

                # Check for relative tolerance
                if(abs(delta_f_relative)<self.optimizationSettings.relative_tolerance_objective):
                    print("\n> Optimization problem converged within a relative objective tolerance of ",self.optimizationSettings.relative_tolerance_objective,"%.")
                    break

                # Check if value of objective increases
                if(response[only_F_id]["value"]>previous_f):
                    print("\n> Value of objective function increased!")
                    break

            # Update design
            X = self.getDesign()

            # Store values of for next iteration
            previous_f = response[only_F_id]["value"]

            # Store initial objective value
            if(opt_itr==1):
                initial_f = response[only_F_id]["value"]

    # --------------------------------------------------------------------------
    def runAugmentedLagrangeAlgorithm( self ):

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
                                                     self.optimizationSettings.penalty_fac_0, 
                                                     self.optimizationSettings.gamma, 
                                                     self.optimizationSettings.penalty_fac_max, 
                                                     self.optimizationSettings.lambda_0 )

        # Initialize file where design evolution is recorded
        with open(self.optimizationSettings.design_history_directory+"/"+self.optimizationSettings.design_history_file, 'w') as csvfile:
            historyWriter = csv.writer(csvfile, delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL)
            row = []
            row.append("itr\t")
            row.append("\tsub_itr\t")
            row.append("\tl\t")
            row.append("\tdl_relative[%]\t")
            row.append("\tf\t")
            row.append("\tdf_absolute[%]\t")
            row.append("\tpenalty_fac\t")
            row.append("\tC["+str(only_C_id)+"]:"+str(self.constraints[only_C_id]["type"])+"\t") 
            row.append("\tlambda["+str(only_C_id)+"]")
            historyWriter.writerow(row)  

        # Define initial design (initial design corresponds to a zero shape update)
        # Note that we assume an incremental design variable in vertex morphing
        X = {}
        for node in self.optimizationModelPart.Nodes:
            X[node.Id] = [0.,0.,0.]   

        # Miscellaneous working variables for data management
        initial_f = 0.0
        initial_l = 0.0
        previous_l = 0.0         

        # Start primary optimization loop
        for opt_itr in range(1,self.optimizationSettings.max_opt_iterations+1):

            # Some output
            print("\n>===================================================================")
            print("> ",time.ctime(),": Starting optimization iteration ",opt_itr)
            print(">===================================================================\n")

            timeAtStartOfCurrentOptimizationStep = time.time()

            # Solve optimization subproblem
            for sub_opt_itr in range(1,self.optimizationSettings.max_sub_opt_iterations+1):

                # Some output
                print("\n>===============================================")
                print("> Starting suboptimization iteration ",sub_opt_itr)
                print(">===============================================\n")

                # Start measuring time needed for current suboptimization step
                subopt_start_time = time.time()

                # Set communicator to evaluate objective and constraint
                self.communicator.initializeControls()
                self.communicator.getControls()[only_F_id]["calculateValue"] = 1
                self.communicator.getControls()[only_C_id]["calculateValue"] = 1

                # Set communicator to evaluate objective and constraint gradient
                self.communicator.getControls()[only_F_id]["calculateGradient"] = 1
                self.communicator.getControls()[only_C_id]["calculateGradient"] = 1                

                # Initialize response container
                response = self.communicator.create_response_container()          

                # Start analyzer according to specified controls
                iterator = str(opt_itr) + "." + str(sub_opt_itr)
                self.analyzer( X, self.communicator.getControls(), iterator, response )

                # Evaluate Lagrange function
                l = self.opt_utils.get_value_of_augmented_lagrangian( only_F_id, self.constraints, response )

                # Store gradients on the nodes of the model_part
                self.storeGradientsOnNodes( response[only_F_id]["gradient"], response[only_C_id]["gradient"] )

                # Compute unit surface normals at each node of current design
                self.geom_utils.compute_unit_surface_normals()

                # Project received gradients on unit surface normal at each node to obtain normal gradients
                self.geom_utils.project_grad_on_unit_surface_normal( constraints_given )

                # Compute mapping matrix
                self.mapper.compute_mapping_matrix()

                # Map sensitivities to design space
                self.mapper.map_sensitivities_to_design_space( constraints_given )

                # Compute search direction
                self.opt_utils.compute_search_direction_augmented_lagrange( self.constraints, response )

                # Compute design variable update (do one step in the optimization algorithm)
                self.opt_utils.compute_design_update()
    
                # Map design update to geometry space
                self.mapper.map_design_update_to_geometry_space()

                # Compute and output some measures to track changes in the objective function
                delta_f_absolute = 0.0
                delta_l_relative = 0.0
                print("\n> Current value of Lagrange function = ",round(l,12))
                if(sub_opt_itr>1):
                    delta_f_absolute = 100*(response[only_F_id]["value"]-initial_f)/initial_f
                    delta_l_relative = 100*(l-previous_l)/initial_l
                    print("\n> Relative change of Lagrange function = ",round(delta_l_relative,6)," [%]")

                # We write every major and every suboptimization iteration in design history
                with open(self.optimizationSettings.design_history_directory+"/"+self.optimizationSettings.design_history_file, 'a') as csvfile:
                   historyWriter = csv.writer(csvfile, delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL)
                   row = []
                   row.append(str(opt_itr)+"\t")
                   row.append("\t"+str(sub_opt_itr)+"\t")
                   row.append("\t"+str("%.12f"%(l))+"\t")
                   row.append("\t"+str("%.6f"%(delta_l_relative))+"\t")
                   row.append("\t"+str("%.12f"%(response[only_F_id]["value"]))+"\t")
                   row.append("\t"+str("%.2f"%(delta_f_absolute))+"\t")
                   row.append("\t"+str("%.2f"%(self.vm_utils.get_penalty_fac()))+"\t")
                   row.append("\t"+str("%.12f"%(response[only_C_id]["value"]))+"\t")
                   row.append("\t"+str("%.12f"%(self.vm_utils.get_lambda(only_C_id))))
                   historyWriter.writerow(row)

                # Write design in GID format
                write_itr = float(iterator)
                self.gidIO.write_results(write_itr, self.optimizationModelPart, self.optimizationSettings.nodal_results, [])

                # Check convergence (We ensure that at least 2 subiterations are done)
                if(sub_opt_itr>3):

                    # Check if maximum iterations were reached
                    if(sub_opt_itr==self.optimizationSettings.max_sub_opt_iterations): 
                        print("\n> Maximal iterations of supoptimization problem reached!")
                        break

                    # Check for relative tolerance
                    if(abs(delta_l_relative)<self.optimizationSettings.relative_tolerance_sub_opt): 
                        print("\n> Optimization subproblem converged within a relative objective tolerance of ",self.optimizationSettings.relative_tolerance_sub_opt,"%.")
                        break

                    # Check if value of lagrangian increases
                    if(l>previous_l):
                        print("\n> Value of Lagrange function increased!")
                        break

                # Update design
                X = self.getDesign()   

                # Store value of Lagrange function for next iteration
                previous_l = l

                # Store initial objective value
                if(opt_itr==1 and sub_opt_itr==1):
                    initial_f = response[only_F_id]["value"]
                    initial_l = l

                # Take time needed for current suboptimization step as well as for the overall opt so far
                subopt_end_time = time.time()
                print("\n> Time needed for current suboptimization step = ",round(subopt_end_time - subopt_start_time,2),"s")
                print("\n> Time needed for total optimization so far = ",round(subopt_end_time - self.timeAtStartOfOptimization,2),"s")

            # Check Convergence (More convergence criterion for major optimization iteration to be implemented!)

            # Check if maximum iterations were reached
            if(opt_itr==self.optimizationSettings.max_opt_iterations):
                print("\n> Maximal iterations of optimization problem reached!")
                break
            
            # Update lagrange multipliers and penalty factor
            self.vm_utils.udpate_augmented_lagrange_parameters( self.constraints, response )  

            # Take time needed for current optimization step
            timeAtEndOfCurrentOptimizationStep = time.time()
            print("\n> Time needed for current optimization step = ",round(timeAtEndOfCurrentOptimizationStep - timeAtStartOfCurrentOptimizationStep,2),"s")

    # --------------------------------------------------------------------------
    def runPenalizedProjectionAlgorithm( self ):

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
        with open(self.optimizationSettings.design_history_directory+"/"+self.optimizationSettings.design_history_file, 'w') as csvfile:
            historyWriter = csv.writer(csvfile, delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL)
            row = []
            row.append("itr\t")
            row.append("\tf\t")
            row.append("\tdf_absolute[%]\t")
            row.append("\tdf_relative[%]\t")
            row.append("\tc["+str(only_C_id)+"]:"+str(self.constraints[only_C_id]["type"])+"\t")    
            row.append("\tc["+str(only_C_id)+"] / reference_value[%]"+"\t")        
            row.append("\tcorrection_scaling[-]\t")
            row.append("\tstep_size[-]\t")
            row.append("\tt_iteration[s]\t")
            row.append("\tt_total[s]") 
            historyWriter.writerow(row)    
        
        # Define initial design (initial design corresponds to a zero shape update)
        # Note that we assume an incremental design variable in vertex morphing
        X = {}
        for node in self.optimizationModelPart.Nodes:
            X[node.Id] = [0.,0.,0.]

        # Miscellaneous working variables for data management
        initial_f = 0.0
        previous_f = 0.0 
        previous_c = 0.0        

        # Start optimization loop
        for opt_itr in range(1,self.optimizationSettings.max_opt_iterations+1):

            # Some output
            print("\n>===================================================================")
            print("> ",time.ctime(),": Starting optimization iteration ",opt_itr)
            print(">===================================================================\n")

            timeAtStartOfCurrentOptimizationStep = time.time()

            # Set communicator to evaluate objective and constraint
            self.communicator.initializeControls()
            self.communicator.getControls()[only_F_id]["calculateValue"] = 1
            self.communicator.getControls()[only_C_id]["calculateValue"] = 1

            # Set communicator to evaluate objective and constraint gradient
            self.communicator.getControls()[only_F_id]["calculateGradient"] = 1
            self.communicator.getControls()[only_C_id]["calculateGradient"] = 1                

            # Initialize response container
            response = self.communicator.create_response_container()          

            # Start analyzer according to specified controls
            iterator = str(opt_itr) + ".0"
            self.analyzer( X, self.communicator.getControls(), iterator, response )

            # Check if constraint is active
            for func_id in self.optimizationSettings.constraints:
                if(self.optimizationSettings.constraints[func_id]["type"] == "eq"):
                    constraints_given = True
                elif(response[func_id]["value"]>0):
                    constraints_given = True
                else:
                    constraints_given = False
                break

            # Store gradients on the nodes of the model_part
            self.storeGradientsOnNodes( response[only_F_id]["gradient"], response[only_C_id]["gradient"] )

            # Compute unit surface normals at each node of current design
            self.geom_utils.compute_unit_surface_normals()

            # Project received gradients on unit surface normal at each node to obtain normal gradients
            self.geom_utils.project_grad_on_unit_surface_normal( constraints_given )

            # Compute mapping matrix
            self.mapper.compute_mapping_matrix()

            # Map sensitivities to design space
            self.mapper.map_sensitivities_to_design_space( constraints_given )

            # # Adjustment of step size
            # if( opt_itr > 1 and response[only_F_id]["value"]>previous_f):
            #     self.optimizationSettings.step_size = self.optimizationSettings.step_size/2

            correction_scaling = [False] 
            if(constraints_given):
                self.opt_utils.compute_projected_search_direction( response[only_C_id]["value"] )
                self.opt_utils.correct_projected_search_direction( response[only_C_id]["value"], previous_c, correction_scaling )
                self.opt_utils.compute_design_update()
            else:
                self.opt_utils.compute_search_direction_steepest_descent()
                self.opt_utils.compute_design_update()

            # Map design update to geometry space
            self.mapper.map_design_update_to_geometry_space()

            # Compute and output some measures to track objective function
            delta_f_absolute = 0.0
            delta_f_relative = 0.0
            print("\n> Current value of objective function = ",response[only_F_id]["value"])
            if(opt_itr>1):
                delta_f_absolute = 100*(response[only_F_id]["value"]-initial_f)/initial_f
                delta_f_relative = 100*(response[only_F_id]["value"]-previous_f)/initial_f
                print("\n> Absolut change of objective function = ",round(delta_f_absolute,6)," [%]")
                print("\n> Relative change of objective function = ",round(delta_f_relative,6)," [%]")  

            # Compute and output some measures to track function
            print("\n> Current value of constraint function = ",round(response[only_C_id]["value"],12))

            # Take time needed for current optimization step
            timeAtEndOfCurrentOptimizationStep = time.time()
            time_current_step = round(timeAtEndOfCurrentOptimizationStep - timeAtStartOfCurrentOptimizationStep,2)
            time_optimization = round(timeAtEndOfCurrentOptimizationStep - self.timeAtStartOfOptimization,2)
            print("\n> Time needed for current optimization step = ",time_current_step,"s")
            print("\n> Time needed for total optimization so far = ",time_optimization,"s")     
            
            # Write design in GID format
            self.gidIO.write_results(opt_itr, self.optimizationModelPart, self.optimizationSettings.nodal_results, [])   

            # Write design history to file
            with open(self.optimizationSettings.design_history_directory+"/"+self.optimizationSettings.design_history_file, 'a') as csvfile:
                historyWriter = csv.writer(csvfile, delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL)
                row = []
                row.append(str(opt_itr)+"\t")
                row.append("\t"+str("%.12f"%(response[only_F_id]["value"]))+"\t")
                row.append("\t"+str("%.2f"%(delta_f_absolute))+"\t")
                row.append("\t"+str("%.6f"%(delta_f_relative))+"\t")
                row.append("\t"+str("%.12f"%(response[only_C_id]["value"]))+"\t")
                if not response[only_C_id]["reference_value"]:
                    row.append("\t"+str("-\t"))
                else: 
                    percentage_of_reference = 100*(response[only_C_id]["value"]/response[only_C_id]["reference_value"])
                    row.append("\t"+str("%.6f"%(percentage_of_reference)))
                row.append("\t"+str("%.12f"%(correction_scaling[0]))+"\t")
                row.append("\t"+str(self.optimizationSettings.step_size)+"\t")
                row.append("\t"+str("%.1f"%(time_current_step))+"\t")
                row.append("\t"+str("%.1f"%(time_optimization)))
                historyWriter.writerow(row)       

            # Check convergence (Further convergence criterions to be implemented )
            if(opt_itr>1):

                # Check if maximum iterations were reached
                if(opt_itr==self.optimizationSettings.max_opt_iterations):
                    print("\n> Maximal iterations of optimization problem reached!")
                    break

                # Check for relative tolerance
                if(abs(delta_f_relative)<self.optimizationSettings.relative_tolerance_objective):
                    print("\n> Optimization problem converged within a relative objective tolerance of ",self.optimizationSettings.relative_tolerance_objective,"%.")
                    break

            # Update design
            X = self.getDesign()

            # Store values of for next iteration
            previous_f = response[only_F_id]["value"]
            previous_c = response[only_C_id]["value"]

            # Store initial objective value
            if(opt_itr==1):
                initial_f = response[only_F_id]["value"]            
    # --------------------------------------------------------------------------
    def storeGradientsOnNodes( self,objective_grads,constraint_grads={} ):

        # Read objective gradients
        eucledian_norm_obj_sens = 0.0
        for node_Id in objective_grads:

            # If deactivated, nodal sensitivities will not be assigned and hence remain zero
            if(self.optimizationModelPart.Nodes[node_Id].GetSolutionStepValue(SENSITIVITIES_DEACTIVATED)):
                continue

            # If not deactivated, nodal sensitivities will be assigned
            sens_i = Vector(3)
            sens_i[0] = objective_grads[node_Id][0]
            sens_i[1] = objective_grads[node_Id][1]
            sens_i[2] = objective_grads[node_Id][2]           
            self.optimizationModelPart.Nodes[node_Id].SetSolutionStepValue(OBJECTIVE_SENSITIVITY,0,sens_i)

        # When constraint_grads is defined also store constraint sensitivities (bool returns false if dictionary is empty)
        if(bool(constraint_grads)):
            eucledian_norm_cons_sens = 0.0
            for node_Id in constraint_grads:

                # If deactivated, nodal sensitivities will not be assigned and hence remain zero
                if(self.optimizationModelPart.Nodes[node_Id].GetSolutionStepValue(SENSITIVITIES_DEACTIVATED)):
                    continue

                # If not deactivated, nodal sensitivities will be assigned
                sens_i = Vector(3)
                sens_i[0] = constraint_grads[node_Id][0]
                sens_i[1] = constraint_grads[node_Id][1]
                sens_i[2] = constraint_grads[node_Id][2]           
                self.optimizationModelPart.Nodes[node_Id].SetSolutionStepValue(CONSTRAINT_SENSITIVITY,0,sens_i)                 

    # --------------------------------------------------------------------------
    def addVariablesNeededForOptimization( self ):
        self.inputModelPart.AddNodalSolutionStepVariable(NORMAL)
        self.inputModelPart.AddNodalSolutionStepVariable(NORMALIZED_SURFACE_NORMAL)
        self.inputModelPart.AddNodalSolutionStepVariable(OBJECTIVE_SENSITIVITY)
        self.inputModelPart.AddNodalSolutionStepVariable(OBJECTIVE_SURFACE_SENSITIVITY)
        self.inputModelPart.AddNodalSolutionStepVariable(MAPPED_OBJECTIVE_SENSITIVITY)
        self.inputModelPart.AddNodalSolutionStepVariable(CONSTRAINT_SENSITIVITY) 
        self.inputModelPart.AddNodalSolutionStepVariable(CONSTRAINT_SURFACE_SENSITIVITY)
        self.inputModelPart.AddNodalSolutionStepVariable(MAPPED_CONSTRAINT_SENSITIVITY) 
        self.inputModelPart.AddNodalSolutionStepVariable(DESIGN_UPDATE)
        self.inputModelPart.AddNodalSolutionStepVariable(DESIGN_CHANGE_ABSOLUTE)  
        self.inputModelPart.AddNodalSolutionStepVariable(SEARCH_DIRECTION) 
        self.inputModelPart.AddNodalSolutionStepVariable(SHAPE_UPDATE) 
        self.inputModelPart.AddNodalSolutionStepVariable(SHAPE_CHANGE_ABSOLUTE)
        self.inputModelPart.AddNodalSolutionStepVariable(IS_ON_BOUNDARY)
        self.inputModelPart.AddNodalSolutionStepVariable(BOUNDARY_PLANE) 
        self.inputModelPart.AddNodalSolutionStepVariable(SHAPE_UPDATES_DEACTIVATED) 
        self.inputModelPart.AddNodalSolutionStepVariable(SENSITIVITIES_DEACTIVATED) 
        
    # --------------------------------------------------------------------------
    def outputInformationAboutResponseFunctions( self ):
        print("> The following objectives are defined:\n")
        for func_id in self.optimizationSettings.objectives:
            print(func_id,":",self.optimizationSettings.objectives[func_id],"\n")

        if( len(self.optimizationSettings.constraints)!=0 ):
            print("> The following constraints are defined:\n")
            for func_id in self.optimizationSettings.constraints:
                print(func_id,":",self.optimizationSettings.constraints[func_id],"\n")
        else:
            print("> No constraints defined.\n")

    # --------------------------------------------------------------------------
    def getDesign( self ):

        # Read and return the current design in the corresponding mode
        X = {}

        if(self.optimizationSettings.design_output_mode=="relative"):
            for node in self.optimizationModelPart.Nodes:
                X[node.Id] = node.GetSolutionStepValue(SHAPE_CHANGE_ABSOLUTE)

        elif(self.optimizationSettings.design_output_mode=="absolute"):
            for node in self.optimizationModelPart.Nodes:
                X[node.Id] = [node.X,node.Y,node.Z]

        else:
            sys.exit("Wrong definition of design_output_mode!")

        return X

    # --------------------------------------------------------------------------

# ==============================================================================
class Communicator:

    # --------------------------------------------------------------------------
    def __init__( self, optimizationSettings ):
        self.initializeListOfRequests( optimizationSettings )
        self.initializeListOfResponses( optimizationSettings )

    # --------------------------------------------------------------------------
    def initializeListOfRequests( self, optimizationSettings ):
        self.listOfRequests = {}
        for functionId in optimizationSettings.objectives:
            self.listOfRequests[functionId] = {"calculateValue": False, "calculateGradient": False}
        for functionId in optimizationSettings.constraints:
            self.listOfRequests[functionId] = {"calculateValue": False, "calculateGradient": False}

    # --------------------------------------------------------------------------
    def initializeListOfResponses( self, optimizationSettings ):
        self.response_container = {}       
        for func_id in optimizationSettings.objectives:
            self.response_container[func_id] = {}
        for func_id in optimizationSettings.constraints:
            self.response_container[func_id] = {}          
        for func_id in self.response_container:
            self.response_container[func_id] = {"value": None, "reference_value": None, "gradient": None}  
    
    # --------------------------------------------------------------------------
    def deleteAllRequests( self ):
        for functionId in self.listOfRequests:
            self.listOfRequests[functionId] = {"calculateValue": False, "calculateGradient": False}

    # --------------------------------------------------------------------------
    def requestFunctionValueOf( self, functionId ):
        self.listOfRequests[functionId]["calculateValue"] = True

    # --------------------------------------------------------------------------
    def requestGradientOf( self, functionId ):
        self.listOfRequests[functionId]["calculateGradient"] = True
    # --------------------------------------------------------------------------
    def isRrequestingFunctionValueOf( self, functionId ):
        return self.listOfRequests[functionId]["calculateValue"]

    # --------------------------------------------------------------------------
    def isRrequestingGradientOf( self, functionId ):
        return self.listOfRequests[functionId]["calculateGradient"]

    # --------------------------------------------------------------------------
    def reportFunctionValue( self, functionId, functionValue ):
        self.response_container[func_id]["value"] = functionValue

    # --------------------------------------------------------------------------
    def reportGradient( self, functionId, gradient ):
        self.response_container[func_id]["gradient"] = gradient

    # --------------------------------------------------------------------------
    def getReportedFunctionValueOf( functionId ):
        return self.response_container[func_id]["value"]

    # --------------------------------------------------------------------------
    def getReportedGradientOf( functionId ):
        return self.response_container[func_id]["gradient"]

    # --------------------------------------------------------------------------

# ==============================================================================
class Analyzer:

    # --------------------------------------------------------------------------
    def __init__( self, importedFunctionForDesignAnalysis, optimizationSettings ):
        self.optimizationSettings = optimizationSettings
        self.importedFunctionForDesignAnalysis = importedFunctionForDesignAnalysis

    def analyzeDesignAndReportToCommunicator( currentDesign, optimizationIteration, communicator )
        self.importedFunctionForDesignAnalysis( currentDesign, optimizationIteration, communicator )

    # -------------------------------------------------------------------------

# ==============================================================================
def CreateOptimizer( inputModelPart, optimizationSettings ):

    # Create folder where all functions includig the optimizer may store their design history in
    os.system( "rm -rf " + optimizationSettings.design_history_directory )
    os.system( "mkdir -p " + optimizationSettings.design_history_directory )

    # Create optimizer according to selected optimization method
    if( optimizationSettings.design_control == "vertex_morphing" ):
        optimizer = VertexMorphingMethod( inputModelPart, optimizationSettings )
        return optimizer

    else:
        sys.exit( "Specified design_control not implemented" )

# ==============================================================================
