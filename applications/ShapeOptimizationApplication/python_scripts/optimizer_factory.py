# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Baumgärtner Daniel, https://github.com/dbaumgaertner
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
import json

# ==============================================================================
def CreateOptimizer( inputModelPart, optimizationSettings ):
    if  optimizationSettings["design_variables"]["design_control_type"].GetString() == "vertex_morphing":
        optimizer = VertexMorphingMethod( inputModelPart, optimizationSettings )
        return optimizer
    else:
        sys.exit( "Specified design_control not implemented" )

# ==============================================================================
class VertexMorphingMethod:
    # --------------------------------------------------------------------------
    def __init__( self, inputModelPart, optimizationSettings ):
        
        self.inputModelPart = inputModelPart
        self.optimizationSettings = optimizationSettings

        self.gidIO = self.createGiDIO( optimizationSettings )
        self.optimizationLogFile = self.createCompleteOptimizationLogFilename ( optimizationSettings )
        self.communicator = Communicator( optimizationSettings )

        self.addVariablesNeededForOptimization( inputModelPart )

    # --------------------------------------------------------------------------
    def createGiDIO( self, optimizationSettings ):
        resultsDirectory = optimizationSettings["ouptut"]["output_directory"].GetString()
        designHistoryFilename = optimizationSettings["ouptut"]["design_history_filename"].GetString()
        designHistoryFilenameWithPath =  resultsDirectory+"/"+designHistoryFilename
        gidIO = GiDOutput( designHistoryFilenameWithPath,
                           optimizationSettings["ouptut"]["output_type"]["VolumeOutput"].GetBool(),
                           optimizationSettings["ouptut"]["output_type"]["GiDPostMode"].GetString(),
                           optimizationSettings["ouptut"]["output_type"]["GiDMultiFileFlag"].GetString(),
                           optimizationSettings["ouptut"]["output_type"]["GiDWriteMeshFlag"].GetBool(),
                           optimizationSettings["ouptut"]["output_type"]["GiDWriteConditionsFlag"].GetBool() )
        return gidIO
    
    # --------------------------------------------------------------------------
    def createCompleteOptimizationLogFilename( self, optimizationSettings ):
        resultsDirectory = optimizationSettings["ouptut"]["output_directory"].GetString()
        optimizationLogFilename = optimizationSettings["ouptut"]["optimization_log_filename"].GetString()
        completeOptimizationLogFilename = resultsDirectory+"/"+optimizationLogFilename+".csv"
        return completeOptimizationLogFilename

    # --------------------------------------------------------------------------
    def addVariablesNeededForOptimization( self, inputModelPart ):
        inputModelPart.AddNodalSolutionStepVariable(NORMAL)
        inputModelPart.AddNodalSolutionStepVariable(NORMALIZED_SURFACE_NORMAL)
        inputModelPart.AddNodalSolutionStepVariable(OBJECTIVE_SENSITIVITY)
        inputModelPart.AddNodalSolutionStepVariable(OBJECTIVE_SURFACE_SENSITIVITY)
        inputModelPart.AddNodalSolutionStepVariable(MAPPED_OBJECTIVE_SENSITIVITY)
        inputModelPart.AddNodalSolutionStepVariable(CONSTRAINT_SENSITIVITY) 
        inputModelPart.AddNodalSolutionStepVariable(CONSTRAINT_SURFACE_SENSITIVITY)
        inputModelPart.AddNodalSolutionStepVariable(MAPPED_CONSTRAINT_SENSITIVITY) 
        inputModelPart.AddNodalSolutionStepVariable(DESIGN_UPDATE)
        inputModelPart.AddNodalSolutionStepVariable(DESIGN_CHANGE_ABSOLUTE)  
        inputModelPart.AddNodalSolutionStepVariable(SEARCH_DIRECTION) 
        inputModelPart.AddNodalSolutionStepVariable(SHAPE_UPDATE) 
        inputModelPart.AddNodalSolutionStepVariable(SHAPE_CHANGE_ABSOLUTE)
        inputModelPart.AddNodalSolutionStepVariable(IS_ON_BOUNDARY)
        inputModelPart.AddNodalSolutionStepVariable(BOUNDARY_PLANE) 
        inputModelPart.AddNodalSolutionStepVariable(SHAPE_UPDATES_DEACTIVATED) 
        inputModelPart.AddNodalSolutionStepVariable(SENSITIVITIES_DEACTIVATED) 

    # --------------------------------------------------------------------------
    def importModelPart( self ):
        model_part_io = ModelPartIO( self.optimizationSettings.input_model_part_name )
        model_part_io.ReadModelPart( self.inputModelPart )
        buffer_size = 1
        self.inputModelPart.SetBufferSize( buffer_size )
        self.inputModelPart.ProcessInfo.SetValue(DOMAIN_SIZE,self.optimizationSettings.domain_size)

    # --------------------------------------------------------------------------
    def importAnalyzer( self, newAnalyzer ): 
        self.analyzer = newAnalyzer

    # --------------------------------------------------------------------------
    def optimize( self ):
        
        print("\n> ==============================================================================================================")
        print("> Starting optimization using the following algorithm: ",self.optimizationSettings["optimization_algorithm"]["name"].GetString())
        print("> ==============================================================================================================\n")

        print(time.ctime())
        self.timeAtStartOfOptimization = time.time()

        self.prepareOptimization()

        iteratorForInitialDesign = 0
        self.gidIO.initialize_results(self.optimizationModelPart)
        self.gidIO.write_results(iteratorForInitialDesign, self.optimizationModelPart, self.optimizationSettings.nodal_results, [])   

        self.runSpecifiedOptimizationAlgorithm()

        self.gidIO.finalize_results()

        print("\n> ==============================================================================================================")
        print("> Finished optimization in ",round(time.time() - self.timeAtStartOfOptimization,2)," s!")
        print("> ==============================================================================================================\n")

    # --------------------------------------------------------------------------
    def prepareOptimization( self ):
        self.outputInformationAboutResponseFunctions()
        self.createFolderToStoreOptimizationHistory()
        self.extractDesignSurfaceAndDampingRegions()
        self.createToolsNecessaryForOptimization()

    # --------------------------------------------------------------------------
    def outputInformationAboutResponseFunctions( self ):
        print("\n> The following objectives are defined:\n")
        for func_id in self.optimizationSettings.objectives:
            print(func_id,":",self.optimizationSettings.objectives[func_id],"\n")
        if len(self.optimizationSettings.constraints) != 0:
            print("> The following constraints are defined:\n")
            for func_id in self.optimizationSettings.constraints:
                print(func_id,":",self.optimizationSettings.constraints[func_id],"\n")
        else:
            print("> No constraints defined.\n")

    # --------------------------------------------------------------------------
    def createFolderToStoreOptimizationHistory ( self ):
        os.system( "rm -rf " + self.optimizationSettings.design_history_directory )
        os.system( "mkdir -p " + self.optimizationSettings.design_history_directory )
   
    # --------------------------------------------------------------------------
    def extractDesignSurfaceAndDampingRegions( self ):
        self.optimizationModelPart = self.getDesignSurfaceFromInputModelPart()
        if self.optimizationSettings.perform_damping:
            self.dampingRegions = self.createContainerWithDampingRegionsAndSettings()
        else:
            self.dampingRegions = []

    # --------------------------------------------------------------------------
    def createToolsNecessaryForOptimization( self ):
        self.vertexMorphingMapper = VertexMorphingMapper( self.optimizationModelPart,
                                                          self.optimizationSettings.filter_function,
                                                          self.optimizationSettings.filter_size,
                                                          self.optimizationSettings.perform_damping,
                                                          self.dampingRegions ) 
        self.optimizationTools = OptimizationUtilities( self.optimizationModelPart,
                                                        self.objectives,
                                                        self.constraints,
                                                        self.optimizationSettings.step_size,
                                                        self.optimizationSettings.normalize_search_direction )
        self.geometryTools = GeometryUtilities( self.optimizationModelPart ) 
    
    # --------------------------------------------------------------------------
    def getDesignSurfaceFromInputModelPart ( self ):
        if self.inputModelPart.HasSubModelPart( self.optimizationSettings.design_surface_submodel_part_name):
            optimizationModel = self.inputModelPart.GetSubModelPart( self.optimizationSettings.design_surface_submodel_part_name )
            print("> The following design surface was defined:\n\n",optimizationModel)
            return optimizationModel
        else:
            raise RuntimeError("Sub-model part specified for optimization does not exist!")
    
    # --------------------------------------------------------------------------
    def createContainerWithDampingRegionsAndSettings ( self ):
        damping_regions = []
        print("> The following damping regions are defined: \n")
        for region in self.optimizationSettings.damping_regions:
            sub_mdpa_name = region[0]
            damp_in_X = region[1]
            damp_in_Y = region[2]
            damp_in_Z = region[3]
            damping_function = region[4]
            damping_radius = region[5]
            if self.inputModelPart.HasSubModelPart(sub_mdpa_name):
                print(region)
                damping_regions.append( [ self.inputModelPart.GetSubModelPart(sub_mdpa_name), 
                                        damp_in_X, 
                                        damp_in_Y, 
                                        damp_in_Z, 
                                        damping_function, 
                                        damping_radius] )
            else:
                raise ValueError("The following sub-model part specified for damping does not exist: ",sub_mdpa_name)         
        return damping_regions
    
    #---------------------------------------------------------------------------
    def runSpecifiedOptimizationAlgorithm( self ):
        if self.optimizationSettings.optimization_algorithm == "steepest_descent":
           self.runSteepestDescentAlgorithm()
        elif self.optimizationSettings.optimization_algorithm == "penalized_projection":
           self.runPenalizedProjectionAlgorithm()           
        else:
            sys.exit("Specified optimization_algorithm not implemented!")

    # --------------------------------------------------------------------------
    def runSteepestDescentAlgorithm( self ):

        # Flags to trigger proper function calls
        constraints_given = False

        # Get Id of objective
        onlyObjectiveFunction = None
        for objectiveFunctionId in self.objectives:
            onlyObjectiveFunction = objectiveFunctionId
            break

        # Initialize file where design evolution is recorded
        with open(self.optimizationLogFile, 'w') as csvfile:
            historyWriter = csv.writer(csvfile, delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL)
            row = []
            row.append("itr\t")
            row.append("\tf\t")
            row.append("\tdf_absolute[%]\t")
            row.append("\tdf_relative[%]\t")
            row.append("\tstep_size[-]\t")
            row.append("\tt_iteration[s]\t")
            row.append("\tt_total[s]") 
            row.append("\ttime_stamp") 
            historyWriter.writerow(row)    

        # Miscellaneous working variables for data management
        initialValueOfObjectiveFunction = 0.0
        previousValueOfObjectiveFunction = 0.0

        # Start optimization loop
        for optimizationIteration in range(1,self.optimizationSettings.max_opt_iterations+1):

            # Some output
            print("\n>===================================================================")
            print("> ",time.ctime(),": Starting optimization iteration ",optimizationIteration)
            print(">===================================================================\n")

            timeAtStartOfCurrentOptimizationStep = time.time()

            self.communicator.initializeCommunication()
            self.communicator.requestFunctionValueOf( onlyObjectiveFunction )
            self.communicator.requestGradientOf( onlyObjectiveFunction )

            self.analyzer.analyzeDesignAndReportToCommunicator( self.optimizationModelPart, optimizationIteration, self.communicator )

            valueOfObjectiveFunction = self.communicator.getReportedFunctionValueOf ( onlyObjectiveFunction )
            gradientOfObjectiveFunction = self.communicator.getReportedGradientOf ( onlyObjectiveFunction )  
                        
            self.storeGradientOnNodes( gradientOfObjectiveFunction )

            self.geometryTools.compute_unit_surface_normals()
            self.geometryTools.project_grad_on_unit_surface_normal( constraints_given )

            self.vertexMorphingMapper.compute_mapping_matrix()
            self.vertexMorphingMapper.map_sensitivities_to_design_space( constraints_given )

            self.optimizationTools.compute_search_direction_steepest_descent()
            self.optimizationTools.compute_design_update()

            self.vertexMorphingMapper.map_design_update_to_geometry_space()         
            
            # Compute and output some measures to track changes in the objective function
            absoluteChangeOfObjectiveValue = 0.0
            relativeChangeOfObjectiveValue = 0.0
            print("\n> Current value of objective function = ",valueOfObjectiveFunction)
            if optimizationIteration > 1:
                absoluteChangeOfObjectiveValue = 100*(valueOfObjectiveFunction-initialValueOfObjectiveFunction) / initialValueOfObjectiveFunction
                relativeChangeOfObjectiveValue = 100*(valueOfObjectiveFunction-previousValueOfObjectiveFunction) / initialValueOfObjectiveFunction
                print("> Absolut change of objective function = ",round(absoluteChangeOfObjectiveValue,6)," [%]")
                print("> Relative change of objective function = ",round(relativeChangeOfObjectiveValue,6)," [%]") 

            # Write design in GID format
            self.gidIO.write_results(optimizationIteration, self.optimizationModelPart, self.optimizationSettings.nodal_results, [])   

            # Write design history to file
            runTimeOptimizationStep = round(time.time() - timeAtStartOfCurrentOptimizationStep,2)
            runTimeOptimization = round(time.time() - self.timeAtStartOfOptimization,2)
            with open(self.optimizationLogFile, 'a') as csvfile:
                historyWriter = csv.writer(csvfile, delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL)
                row = []
                row.append(str(optimizationIteration)+"\t")
                row.append("\t"+str("%.12f"%(valueOfObjectiveFunction))+"\t")
                row.append("\t"+str("%.2f"%(absoluteChangeOfObjectiveValue))+"\t")
                row.append("\t"+str("%.6f"%(relativeChangeOfObjectiveValue))+"\t")
                row.append("\t"+str(self.optimizationSettings.step_size)+"\t")
                row.append("\t"+str("%.1f"%(runTimeOptimizationStep))+"\t")
                row.append("\t"+str("%.1f"%(runTimeOptimization))+"\t")
                row.append("\t"+str(time.ctime()))
                historyWriter.writerow(row)     

            # Take time needed for current optimization step
            print("\n> Time needed for current optimization step = ",runTimeOptimizationStep,"s")
            print("> Time needed for total optimization so far = ",runTimeOptimization,"s")                              

            # Check convergence
            if optimizationIteration > 1 :

                # Check if maximum iterations were reached
                if optimizationIteration == self.optimizationSettings.max_opt_iterations:
                    print("\n> Maximal iterations of optimization problem reached!")
                    break

                # Check for relative tolerance
                if abs(relativeChangeOfObjectiveValue) < self.optimizationSettings.relative_tolerance_objective:
                    print("\n> Optimization problem converged within a relative objective tolerance of ",self.optimizationSettings.relative_tolerance_objective,"%.")
                    break

                # Check if value of objective increases
                if valueOfObjectiveFunction > previousValueOfObjectiveFunction:
                    print("\n> Value of objective function increased!")

            # Store values of for next iteration
            previousValueOfObjectiveFunction = valueOfObjectiveFunction

            # Store initial objective value
            if optimizationIteration == 1:
                initialValueOfObjectiveFunction = valueOfObjectiveFunction

    # --------------------------------------------------------------------------
    def runPenalizedProjectionAlgorithm( self ):

        # README!!!
        # Note that the current implementation assumes that only one scalar objective & constraint is present

        # Flags to trigger proper function calls
        constraints_given = True

        # Get Id of objective
        onlyObjectiveFunction = None
        for objectiveFunctionId in self.objectives:
            onlyObjectiveFunction = objectiveFunctionId
            break

        # Get Id of constraint
        onlyConstraintFunction = None
        for constraintFunctionId in self.constraints:
            onlyConstraintFunction = constraintFunctionId
            break            

        # Initialize file where design evolution is recorded
        with open(self.optimizationLogFile, 'w') as csvfile:
            historyWriter = csv.writer(csvfile, delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL)
            row = []
            row.append("itr\t")
            row.append("\tf\t")
            row.append("\tdf_absolute[%]\t")
            row.append("\tdf_relative[%]\t")
            row.append("\tc["+str(onlyConstraintFunction)+"]:"+str(self.constraints[onlyConstraintFunction]["type"])+"\t")    
            row.append("\tc["+str(onlyConstraintFunction)+"] / reference_value[%]"+"\t")        
            row.append("\tcorrection_scaling[-]\t")
            row.append("\tstep_size[-]\t")
            row.append("\tt_iteration[s]\t")
            row.append("\tt_total[s]") 
            row.append("\ttime_stamp") 
            historyWriter.writerow(row)    

        # Miscellaneous working variables for data management
        initialValueOfObjectiveFunction = 0.0
        previousValueOfObjectiveFunction = 0.0 
        previousValueOfConstraintFunction = 0.0        

        # Start optimization loop
        for optimizationIteration in range(1,self.optimizationSettings.max_opt_iterations+1):

            # Some output
            print("\n>===================================================================")
            print("> ",time.ctime(),": Starting optimization iteration ",optimizationIteration)
            print(">===================================================================\n")

            timeAtStartOfCurrentOptimizationStep = time.time()   

            self.communicator.initializeCommunication()
            self.communicator.requestFunctionValueOf( onlyObjectiveFunction )
            self.communicator.requestFunctionValueOf( onlyConstraintFunction )
            self.communicator.requestGradientOf( onlyObjectiveFunction )
            self.communicator.requestGradientOf( onlyConstraintFunction )

            self.analyzer.analyzeDesignAndReportToCommunicator( self.optimizationModelPart, optimizationIteration, self.communicator )

            valueOfObjectiveFunction = self.communicator.getReportedFunctionValueOf ( onlyObjectiveFunction )
            valueOfConstraintFunction = self.communicator.getReportedFunctionValueOf ( onlyConstraintFunction )
            referenceValueOfConstraintFunction = self.communicator.getReportedFunctionReferenceValueOf ( onlyConstraintFunction )
            gradientOfObjectiveFunction = self.communicator.getReportedGradientOf ( onlyObjectiveFunction )
            gradientOfConstraintFunction = self.communicator.getReportedGradientOf ( onlyConstraintFunction )

            self.storeGradientOnNodes( gradientOfObjectiveFunction, gradientOfConstraintFunction )
            
            # Check if constraint is active
            for functionId in self.constraints:
                if self.constraints[functionId]["type"] == "eq":
                    constraints_given = True
                elif valueOfConstraintFunction > 0:
                    constraints_given = True
                else:
                    constraints_given = False
                break           

            self.geometryTools.compute_unit_surface_normals()
            self.geometryTools.project_grad_on_unit_surface_normal( constraints_given )

            self.vertexMorphingMapper.compute_mapping_matrix()
            self.vertexMorphingMapper.map_sensitivities_to_design_space( constraints_given )

            correctionScaling = [False] 
            if constraints_given:
                self.optimizationTools.compute_projected_search_direction( valueOfConstraintFunction )
                self.optimizationTools.correct_projected_search_direction( valueOfConstraintFunction, previousValueOfConstraintFunction, correctionScaling )
                self.optimizationTools.compute_design_update()
            else:
                self.optimizationTools.compute_search_direction_steepest_descent()
                self.optimizationTools.compute_design_update()

            self.vertexMorphingMapper.map_design_update_to_geometry_space()  
            
            # Compute and output some measures to track response functions
            absoluteChangeOfObjectiveValue = 0.0
            relativeChangeOfObjectiveValue = 0.0
            print("\n> Current value of objective function = ",valueOfObjectiveFunction)
            if optimizationIteration > 1:
                absoluteChangeOfObjectiveValue = 100*(valueOfObjectiveFunction-initialValueOfObjectiveFunction) / initialValueOfObjectiveFunction
                relativeChangeOfObjectiveValue = 100*(valueOfObjectiveFunction-previousValueOfObjectiveFunction) / initialValueOfObjectiveFunction
                print("> Absolut change of objective function = ",round(absoluteChangeOfObjectiveValue,6)," [%]")
                print("> Relative change of objective function = ",round(relativeChangeOfObjectiveValue,6)," [%]")  
            print("\n> Current value of constraint function = ",round(valueOfConstraintFunction,12))     

            # Write design in GID format
            self.gidIO.write_results(optimizationIteration, self.optimizationModelPart, self.optimizationSettings.nodal_results, [])   

            # Write design history to file
            runTimeOptimizationStep = round(time.time() - timeAtStartOfCurrentOptimizationStep,2)
            runTimeOptimization = round(time.time() - self.timeAtStartOfOptimization,2)            
            with open(self.optimizationLogFile, 'a') as csvfile:
                historyWriter = csv.writer(csvfile, delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL)
                row = []
                row.append(str(optimizationIteration)+"\t")
                row.append("\t"+str("%.12f"%(valueOfObjectiveFunction))+"\t")
                row.append("\t"+str("%.2f"%(absoluteChangeOfObjectiveValue))+"\t")
                row.append("\t"+str("%.6f"%(relativeChangeOfObjectiveValue))+"\t")
                row.append("\t"+str("%.12f"%(valueOfConstraintFunction))+"\t")
                if not referenceValueOfConstraintFunction:
                    row.append("\t"+str("-\t"))
                else: 
                    percentageOfReference = 100*(valueOfConstraintFunction / referenceValueOfConstraintFunction)
                    row.append("\t"+str("%.6f"%(percentageOfReference)))
                row.append("\t"+str("%.12f"%(correctionScaling[0]))+"\t")
                row.append("\t"+str(self.optimizationSettings.step_size)+"\t")
                row.append("\t"+str("%.1f"%(runTimeOptimizationStep))+"\t")
                row.append("\t"+str("%.1f"%(runTimeOptimization))+"\t")
                row.append("\t"+str(time.ctime()))
                historyWriter.writerow(row)   

            # Take time needed for current optimization step
            print("\n> Time needed for current optimization step = ",runTimeOptimizationStep,"s")
            print("\n> Time needed for total optimization so far = ",runTimeOptimization,"s")                       

            # Check convergence (Further convergence criterions to be implemented )
            if optimizationIteration>1:

                # Check if maximum iterations were reached
                if optimizationIteration==self.optimizationSettings.max_opt_iterations:
                    print("\n> Maximal iterations of optimization problem reached!")
                    break

                # Check for relative tolerance
                if abs(relativeChangeOfObjectiveValue)<self.optimizationSettings.relative_tolerance_objective:
                    print("\n> Optimization problem converged within a relative objective tolerance of ",self.optimizationSettings.relative_tolerance_objective,"%.")
                    break

                # Check if value of objective increases
                if  valueOfObjectiveFunction > previousValueOfObjectiveFunction:
                    print("\n> Value of objective function increased!")                    

            # Store values of for next iteration
            previousValueOfObjectiveFunction = valueOfObjectiveFunction
            previousValueOfConstraintFunction = valueOfConstraintFunction

            # Store initial objective value
            if optimizationIteration == 1:
                initialValueOfObjectiveFunction = valueOfObjectiveFunction      

    # --------------------------------------------------------------------------
    def storeGradientOnNodes( self,objective_grads,constraint_grads={} ):

        # Read objective gradients
        eucledian_norm_obj_sens = 0.0
        for node_Id in objective_grads:

            # If deactivated, nodal sensitivities will not be assigned and hence remain zero
            if self.optimizationModelPart.Nodes[node_Id].GetSolutionStepValue(SENSITIVITIES_DEACTIVATED):
                continue

            # If not deactivated, nodal sensitivities will be assigned
            sens_i = Vector(3)
            sens_i[0] = objective_grads[node_Id][0]
            sens_i[1] = objective_grads[node_Id][1]
            sens_i[2] = objective_grads[node_Id][2]           
            self.optimizationModelPart.Nodes[node_Id].SetSolutionStepValue(OBJECTIVE_SENSITIVITY,0,sens_i)

        # When constraint_grads is defined also store constraint sensitivities (bool returns false if dictionary is empty)
        if bool(constraint_grads):
            eucledian_norm_cons_sens = 0.0
            for node_Id in constraint_grads:

                # If deactivated, nodal sensitivities will not be assigned and hence remain zero
                if self.optimizationModelPart.Nodes[node_Id].GetSolutionStepValue(SENSITIVITIES_DEACTIVATED):
                    continue

                # If not deactivated, nodal sensitivities will be assigned
                sens_i = Vector(3)
                sens_i[0] = constraint_grads[node_Id][0]
                sens_i[1] = constraint_grads[node_Id][1]
                sens_i[2] = constraint_grads[node_Id][2]           
                self.optimizationModelPart.Nodes[node_Id].SetSolutionStepValue(CONSTRAINT_SENSITIVITY,0,sens_i)                 

    # --------------------------------------------------------------------------

# ==============================================================================
class Communicator:
    # --------------------------------------------------------------------------
    def __init__( self, optimizationSettings ):
        self.listOfRequests = self.initializeListOfRequests( optimizationSettings )
        self.listOfResponses = self.initializeListOfResponses( optimizationSettings )

    # --------------------------------------------------------------------------
    def initializeListOfRequests( self, optimizationSettings ):

        listOfRequests = {}
        numberOfObjectives = optimizationSettings["objectives"].size()
        numberOfConstraints = optimizationSettings["constraints"].size()

        for objectiveNumber in range(numberOfObjectives):
            objectiveId = optimizationSettings["objectives"][objectiveNumber]["identifier"].GetString()
            listOfRequests[objectiveId] = {"calculateValue": False, "calculateGradient": False}

        for constraintNumber in range(numberOfConstraints):
            constraintId = optimizationSettings["constraints"][constraintNumber]["identifier"].GetString()
            listOfRequests[constraintId] = {"calculateValue": False, "calculateGradient": False}

        return listOfRequests

    # --------------------------------------------------------------------------
    def initializeListOfResponses( self, optimizationSettings ):
        
        listOfResponses = {}  
        numberOfObjectives = optimizationSettings["objectives"].size()
        numberOfConstraints = optimizationSettings["constraints"].size()

        for objectiveNumber in range(numberOfObjectives):
            objectiveId = optimizationSettings["objectives"][objectiveNumber]["identifier"].GetString()
            listOfResponses[objectiveId] = {}

        for constraintNumber in range(numberOfConstraints):
            constraintId = optimizationSettings["constraints"][constraintNumber]["identifier"].GetString()
            listOfResponses[constraintId] = {}  

        for responseId in listOfResponses:
            listOfResponses[responseId] = {"value": None, "referenceValue": None, "gradient": None} 

        return listOfResponses
    
    # --------------------------------------------------------------------------
    def initializeCommunication( self ):
        self.deleteAllRequests()
        self.deleteAllReportedValues()

    # --------------------------------------------------------------------------
    def deleteAllRequests( self ):
        for responseId in self.listOfRequests:
            self.listOfRequests[responseId] = {"calculateValue": False, "calculateGradient": False}

    # --------------------------------------------------------------------------
    def deleteAllReportedValues( self ):
        for responseId in self.listOfResponses:
            self.listOfResponses[responseId] = {"value": None, "referenceValue": None, "gradient": None}  

    # --------------------------------------------------------------------------
    def requestFunctionValueOf( self, responseId ):
        self.listOfRequests[responseId]["calculateValue"] = True

    # --------------------------------------------------------------------------
    def requestGradientOf( self, responseId ):
        self.listOfRequests[responseId]["calculateGradient"] = True

    # --------------------------------------------------------------------------
    def isRequestingFunctionValueOf( self, responseId ):
        return self.listOfRequests[responseId]["calculateValue"]

    # --------------------------------------------------------------------------
    def isRequestingGradientOf( self, responseId ):
        return self.listOfRequests[responseId]["calculateGradient"]

    # --------------------------------------------------------------------------
    def reportFunctionValue( self, responseId, functionValue ):
        if responseId in self.listOfResponses.keys():
            self.listOfResponses[responseId]["value"] = functionValue
        else:
            raise NameError("Reported function is not specified: " + responseId)

    # --------------------------------------------------------------------------
    def reportFunctionReferenceValue( self, responseId, functionReferenceValue ):
        if responseId in self.listOfResponses.keys():
            self.listOfResponses[responseId]["referenceValue"] = functionReferenceValue
        else:
            raise NameError("Reported function is not specified: " + responseId)

    # --------------------------------------------------------------------------
    def reportGradient( self, responseId, gradient ):
        if responseId in self.listOfResponses.keys():        
            self.listOfResponses[responseId]["gradient"] = gradient
        else:
            raise NameError("Reported function is not specified: " + responseId)

    # --------------------------------------------------------------------------
    def getReportedFunctionValueOf( self, responseId ):
        return self.listOfResponses[responseId]["value"]

    # --------------------------------------------------------------------------
    def getReportedFunctionReferenceValueOf( self, responseId ):
        return self.listOfResponses[responseId]["referenceValue"]

    # --------------------------------------------------------------------------    
    def getReportedGradientOf( self, responseId ):
        return self.listOfResponses[responseId]["gradient"]

    # --------------------------------------------------------------------------

# ==============================================================================
class analyzerBaseClass:
    def analyzeDesignAndReportToCommunicator( self, currentDesign, optimizationIteration, communicator ):
        raise RuntimeError("Analyzer base class is called. Please check your implementation of the function >> analyzeDesignAndReportToCommunicator << .")

# ==============================================================================
