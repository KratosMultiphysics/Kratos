# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Baumgärtner Daniel, https://github.com/dbaumgaertner
#
# ==============================================================================

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
    if  optimizationSettings["design_variables"]["variable_type"].GetString() == "vertex_morphing":
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
        self.addVariablesNeededForOptimization( inputModelPart )
    
    # --------------------------------------------------------------------------
    def createFolderToStoreOptimizationResults ( self, optimizationSettings ):
        resultsDirectory = optimizationSettings["output"]["output_directory"].GetString()
        os.system( "rm -rf " + resultsDirectory )
        os.system( "mkdir -p " + resultsDirectory )

    # --------------------------------------------------------------------------
    def createGiDIO( self, optimizationSettings ):
        resultsDirectory = optimizationSettings["output"]["output_directory"].GetString()
        designHistoryFilename = optimizationSettings["output"]["design_history_filename"].GetString()
        designHistoryFilenameWithPath =  resultsDirectory+"/"+designHistoryFilename
        gidIO = GiDOutput( designHistoryFilenameWithPath,
                           optimizationSettings["output"]["output_type"]["VolumeOutput"].GetBool(),
                           optimizationSettings["output"]["output_type"]["GiDPostMode"].GetString(),
                           optimizationSettings["output"]["output_type"]["GiDMultiFileFlag"].GetString(),
                           optimizationSettings["output"]["output_type"]["GiDWriteMeshFlag"].GetBool(),
                           optimizationSettings["output"]["output_type"]["GiDWriteConditionsFlag"].GetBool() )
        return gidIO
    
    # --------------------------------------------------------------------------
    def createCompleteOptimizationLogFilename( self, optimizationSettings ):
        resultsDirectory = optimizationSettings["output"]["output_directory"].GetString()
        optimizationLogFilename = optimizationSettings["output"]["optimization_log_filename"].GetString()
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
        model_part_io = ModelPartIO( self.optimizationSettings["design_variables"]["input_model_part_name"].GetString() )
        model_part_io.ReadModelPart( self.inputModelPart )
        buffer_size = 1
        self.inputModelPart.SetBufferSize( buffer_size )
        self.inputModelPart.ProcessInfo.SetValue( DOMAIN_SIZE, self.optimizationSettings["design_variables"]["domain_size"].GetInt() )

    # --------------------------------------------------------------------------
    def importAnalyzer( self, newAnalyzer ): 
        self.analyzer = newAnalyzer

    # --------------------------------------------------------------------------
    def optimize( self ):
        
        print("\n> ==============================================================================================================")
        print("> ",time.ctime(),": Starting optimization using the following algorithm: ",self.optimizationSettings["optimization_algorithm"]["name"].GetString())
        print("> ==============================================================================================================\n")

        
        self.outputInformationAboutResponseFunctions()

        # timer
        self.timeAtStartOfOptimization = time.time()

        # model data
        self.designSurface = self.getDesignSurfaceFromInputModelPart()

        # communicator
        self.communicator = Communicator( self.optimizationSettings )        

        # mapper
        listOfDampingRegions = self.getListOfDampingRegionsFromInputModelPart()
        self.vertexMorphingMapper = VertexMorphingMapper( self.designSurface, listOfDampingRegions, self.optimizationSettings["design_variables"] ) 

        # dataWriter
        optimizatoinDataWriterFactory = __import__("optimization_data_writer_factory")
        optimizationDataWriter = optimizatoinDataWriterFactory.CreateDataWriter( self.designSurface, self.optimizationSettings )

        
        # optimizationDataWriter.initializeWriting()


        



        # outputWriter
        # self.createFolderToStoreOptimizationResults( self.optimizationSettings )
        # self.gidIO = self.createGiDIO( self.optimizationSettings )
        # self.nodalResults = self.generateListOfNodalResults( self.optimizationSettings["output"] )
        # self.optimizationLogFile = self.createCompleteOptimizationLogFilename ( self.optimizationSettings )   

        # algorithm_driver
        # specifiedAlgorithm = OptimizationAlgorithems( self.design_surface, communicator, mapper, outputWriter, timer)

        # outputWriter.initializeOutput()

        # specifiedAlgorithm.run()

         # outputWriter.finalizeOutput()


        # Old format

        
        # iteratorForInitialDesign = 0
        # self.gidIO.initialize_results( self.designSurface )
        # self.gidIO.write_results(iteratorForInitialDesign, self.designSurface, self.nodalResults, [])   

        # self.optimizationTools = OptimizationUtilities( self.designSurface, self.optimizationSettings )
        # self.runSpecifiedOptimizationAlgorithm()

        # self.gidIO.finalize_results()

        print("\n> ==============================================================================================================")
        print("> Finished optimization in ",round(time.time() - self.timeAtStartOfOptimization,2)," s!")
        print("> ==============================================================================================================\n")

    # --------------------------------------------------------------------------
    def outputInformationAboutResponseFunctions( self ):

        numberOfObjectives = self.optimizationSettings["objectives"].size()
        numberOfConstraints = self.optimizationSettings["constraints"].size()

        print("\n> The following objectives are defined:\n")
        for objectiveNumber in range(numberOfObjectives):
            print(self.optimizationSettings["objectives"][objectiveNumber])

        if numberOfConstraints != 0:
            print("> The following constraints are defined:\n")
            for constraintNumber in range(numberOfConstraints):
                print(self.optimizationSettings["constraints"][constraintNumber],"\n")
        else:
            print("> No constraints defined.\n")     
    
    # --------------------------------------------------------------------------
    def getDesignSurfaceFromInputModelPart( self ):
        nameOfDesingSurface = self.optimizationSettings["design_variables"]["design_submodel_part_name"].GetString()
        if self.inputModelPart.HasSubModelPart( nameOfDesingSurface ):
            optimizationModel = self.inputModelPart.GetSubModelPart( nameOfDesingSurface )
            print("> The following design surface was defined:\n\n",optimizationModel)
            return optimizationModel
        else:
            raise RuntimeError("Sub-model part specified for optimization does not exist!")

    # --------------------------------------------------------------------------
    def getListOfDampingRegionsFromInputModelPart( self ):
        listOfDampingRegions = {}
        print("> The following damping regions are defined: \n")
        for regionNumber in range(self.optimizationSettings["design_variables"]["damping"]["damping_regions"].size()):
            regionName = self.optimizationSettings["design_variables"]["damping"]["damping_regions"][regionNumber]["sub_model_part_name"].GetString()
            if self.inputModelPart.HasSubModelPart(regionName):
                print(regionName)
                listOfDampingRegions[regionName] = self.inputModelPart.GetSubModelPart(regionName)
            else:
                raise ValueError("The following sub-model part specified for damping does not exist: ",regionName)         
        return listOfDampingRegions

    # ---------------------------------------------------------------------------
    def generateListOfNodalResults( self, outputSettings ):
        listOfNodalResults = []
        for nodalResultNumber in range(outputSettings["nodal_results"].size()):
            listOfNodalResults.append( outputSettings["nodal_results"][nodalResultNumber].GetString() )
        return listOfNodalResults

    # --------------------------------------------------------------------------
    def runSpecifiedOptimizationAlgorithm( self ):
        if self.optimizationSettings["optimization_algorithm"]["name"].GetString() == "steepest_descent":
           self.runSteepestDescentAlgorithm()
        elif self.optimizationSettings["optimization_algorithm"]["name"].GetString() == "penalized_projection":
           self.runPenalizedProjectionAlgorithm()           
        else:
            sys.exit("Specified optimization_algorithm not implemented!")

    # --------------------------------------------------------------------------
    def runSteepestDescentAlgorithm( self ):

        # README!!!
        # Note that the current implementation assumes that only one scalar objective is present

        # Flags to trigger proper function calls
        constraints_given = False

        # Tools to perform geometric operations
        geometryTools = GeometryUtilities( self.designSurface )

        # Get Id of objective (right now only one objective is considered!!)
        onlyObjectiveFunction = self.optimizationSettings["objectives"][0]["identifier"].GetString()

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
        maxOptimizationIterations = self.optimizationSettings["optimization_algorithm"]["max_iterations"].GetInt() + 1 
        for optimizationIteration in range(1,maxOptimizationIterations):

            # Some output
            print("\n>===================================================================")
            print("> ",time.ctime(),": Starting optimization iteration ",optimizationIteration)
            print(">===================================================================\n")

            timeAtStartOfCurrentOptimizationStep = time.time()

            self.communicator.initializeCommunication()
            self.communicator.requestFunctionValueOf( onlyObjectiveFunction )
            self.communicator.requestGradientOf( onlyObjectiveFunction )

            self.analyzer.analyzeDesignAndReportToCommunicator( self.designSurface, optimizationIteration, self.communicator )

            valueOfObjectiveFunction = self.communicator.getReportedFunctionValueOf ( onlyObjectiveFunction )
            gradientOfObjectiveFunction = self.communicator.getReportedGradientOf ( onlyObjectiveFunction )  
                        
            self.storeGradientOnNodes( gradientOfObjectiveFunction )

            geometryTools.compute_unit_surface_normals()
            geometryTools.project_grad_on_unit_surface_normal( constraints_given )

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
            self.gidIO.write_results(optimizationIteration, self.designSurface, self.nodalResults, [])   

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
                row.append("\t"+str(self.optimizationSettings["line_search"]["step_size"].GetDouble())+"\t")
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
                if optimizationIteration == maxOptimizationIterations:
                    print("\n> Maximal iterations of optimization problem reached!")
                    break

                # Check for relative tolerance
                relativeTolerance = self.optimizationSettings["optimization_algorithm"]["relative_tolerance"].GetDouble()
                if abs(relativeChangeOfObjectiveValue) < relativeTolerance:
                    print("\n> Optimization problem converged within a relative objective tolerance of ",relativeTolerance,"%.")
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

        # Tools to perform geometric operations
        geometryTools = GeometryUtilities( self.designSurface )        

        # Get Id of objective
        onlyObjectiveFunction = self.optimizationSettings["objectives"][0]["identifier"].GetString()

        # Get Id of constraint
        onlyConstraintFunction = self.optimizationSettings["constraints"][0]["identifier"].GetString()
        typOfOnlyConstraintFunction = self.optimizationSettings["constraints"][0]["type"].GetString()      

        # Initialize file where design evolution is recorded
        with open(self.optimizationLogFile, 'w') as csvfile:
            historyWriter = csv.writer(csvfile, delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL)
            row = []
            row.append("itr\t")
            row.append("\tf\t")
            row.append("\tdf_absolute[%]\t")
            row.append("\tdf_relative[%]\t")
            row.append("\tc["+onlyConstraintFunction+"]: "+typOfOnlyConstraintFunction+"\t")    
            row.append("\tc["+onlyConstraintFunction+"] / reference_value[%]"+"\t")        
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
        maxOptimizationIterations = self.optimizationSettings["optimization_algorithm"]["max_iterations"].GetInt() + 1 
        for optimizationIteration in range(1,maxOptimizationIterations):

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

            self.analyzer.analyzeDesignAndReportToCommunicator( self.designSurface, optimizationIteration, self.communicator )

            valueOfObjectiveFunction = self.communicator.getReportedFunctionValueOf ( onlyObjectiveFunction )
            valueOfConstraintFunction = self.communicator.getReportedFunctionValueOf ( onlyConstraintFunction )
            referenceValueOfConstraintFunction = self.communicator.getReportedFunctionReferenceValueOf ( onlyConstraintFunction )
            gradientOfObjectiveFunction = self.communicator.getReportedGradientOf ( onlyObjectiveFunction )
            gradientOfConstraintFunction = self.communicator.getReportedGradientOf ( onlyConstraintFunction )

            self.storeGradientOnNodes( gradientOfObjectiveFunction, gradientOfConstraintFunction )
            
            # Check if constraint is active
            if typOfOnlyConstraintFunction == "equality":
                constraints_given = True
            elif valueOfConstraintFunction > 0:
                constraints_given = True
            else:
                constraints_given = False       

            geometryTools.compute_unit_surface_normals()
            geometryTools.project_grad_on_unit_surface_normal( constraints_given )

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
            self.gidIO.write_results(optimizationIteration, self.designSurface, self.nodalResults, [])   

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
                row.append("\t"+str(self.optimizationSettings["line_search"]["step_size"].GetDouble())+"\t")
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
                if optimizationIteration == maxOptimizationIterations:
                    print("\n> Maximal iterations of optimization problem reached!")
                    break

                # Check for relative tolerance
                relativeTolerance = self.optimizationSettings["optimization_algorithm"]["relative_tolerance"].GetDouble()
                if abs(relativeChangeOfObjectiveValue) < relativeTolerance:
                    print("\n> Optimization problem converged within a relative objective tolerance of ",relativeTolerance,"%.")
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
            if self.designSurface.Nodes[node_Id].GetSolutionStepValue(SENSITIVITIES_DEACTIVATED):
                continue

            # If not deactivated, nodal sensitivities will be assigned
            sens_i = Vector(3)
            sens_i[0] = objective_grads[node_Id][0]
            sens_i[1] = objective_grads[node_Id][1]
            sens_i[2] = objective_grads[node_Id][2]           
            self.designSurface.Nodes[node_Id].SetSolutionStepValue(OBJECTIVE_SENSITIVITY,0,sens_i)

        # When constraint_grads is defined also store constraint sensitivities (bool returns false if dictionary is empty)
        if bool(constraint_grads):
            eucledian_norm_cons_sens = 0.0
            for node_Id in constraint_grads:

                # If deactivated, nodal sensitivities will not be assigned and hence remain zero
                if self.designSurface.Nodes[node_Id].GetSolutionStepValue(SENSITIVITIES_DEACTIVATED):
                    continue

                # If not deactivated, nodal sensitivities will be assigned
                sens_i = Vector(3)
                sens_i[0] = constraint_grads[node_Id][0]
                sens_i[1] = constraint_grads[node_Id][1]
                sens_i[2] = constraint_grads[node_Id][2]           
                self.designSurface.Nodes[node_Id].SetSolutionStepValue(CONSTRAINT_SENSITIVITY,0,sens_i)                 

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
