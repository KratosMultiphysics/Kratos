# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Suneth Warnakulasuriya, https://github.com/sunethwarna
#
# ==============================================================================

# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.ALEApplication import *
from KratosMultiphysics.ShapeOptimizationApplication import *
import KratosMultiphysics.AdjointFluidApplication as AdjointFluidApplication

import time as timer


class kratosAdjointFluidAnalyzer( (__import__("analyzer_base")).analyzerBaseClass ):
    # --------------------------------------------------------------------------   
    def __init__ ( self, inputModelPart, optimizationSettings ):
        self.main__optimization_model_part = inputModelPart
        self.analyzer_settings = optimizationSettings
        self.history_path = self.analyzer_settings["adjoint_fluid_analyzer"]["history_path"].GetString()

        if not os.path.isdir(self.history_path):
            os.mkdir(self.history_path)
            
        setup_factory = __import__("setup_factory")

        self.primal = setup_factory.OptimizationSetup(self.main__optimization_model_part, self.analyzer_settings["adjoint_fluid_analyzer"]["primal_setup"])
        self.adjoint = setup_factory.OptimizationSetup(self.main__optimization_model_part, self.analyzer_settings["adjoint_fluid_analyzer"]["adjoint_setup"])            

    def Initialize( self, Model ):
        self.Model = Model
        self.initializeGIDOutput()
        self.initializeProcesses(Model)
        self.initializeSolutionLoop()   

    # --------------------------------------------------------------------------   
    def SetMeshSolver( self, inputMeshSolver ):
        self.mesh_solver = inputMeshSolver

    def ImportModelPart(self):
    
        if(self.analyzer_settings["adjoint_fluid_analyzer"]["model_import_settings"]["input_type"].GetString() == "mdpa"):
            # here it would be the place to import restart data if required
            ModelPartIO(self.analyzer_settings["adjoint_fluid_analyzer"]["model_import_settings"]["input_filename"].GetString()).ReadModelPart(self.main__optimization_model_part)

            # here we shall check that the input read has the shape we like
            aux_params = Parameters("{}")
            aux_params.AddValue("volume_model_part_name",self.analyzer_settings["adjoint_fluid_analyzer"]["model_import_settings"]["volume_model_part_name"])
            aux_params.AddValue("skin_parts",self.analyzer_settings["adjoint_fluid_analyzer"]["sub_model_part_list"])

            # here we replace the dummy elements we read with proper elements
            self.analyzer_settings.AddEmptyValue("element_replace_settings")
            if(self.main__optimization_model_part.ProcessInfo[DOMAIN_SIZE] == 3):
                self.analyzer_settings["element_replace_settings"] = Parameters("""
                    {
                    "element_name": "VMSAdjointElement3D",
                    "condition_name": "SurfaceCondition3D3N"
                    }
                    """)
            elif(self.main__optimization_model_part.ProcessInfo[DOMAIN_SIZE] == 2):
                self.analyzer_settings["element_replace_settings"] = Parameters("""
                    {
                    "element_name": "VMSAdjointElement2D",
                    "condition_name": "LineCondition2D2N"
                    }
                    """)
            else:
                raise Exception("domain size is not 2 or 3")

            ReplaceElementsAndConditionsProcess(self.main__optimization_model_part, self.analyzer_settings["element_replace_settings"]).Execute()

            # import check_and_prepare_model_process_fluid
            # check_and_prepare_model_process_fluid.CheckAndPrepareModelProcess(self.main__optimization_model_part, aux_params).Execute()

            #here we read the KINEMATIC VISCOSITY and DENSITY and we apply it to the nodes
            # for el in self.main__optimization_model_part.Elements:
            #     rho = el.Properties.GetValue(KratosMultiphysics.DENSITY)
            #     kin_viscosity = el.Properties.GetValue(KratosMultiphysics.VISCOSITY)
            #     break

            # VariableUtils().SetScalarVar(KratosMultiphysics.DENSITY, rho, self.main__optimization_model_part.Nodes)
            # VariableUtils().SetScalarVar(KratosMultiphysics.VISCOSITY, kin_viscosity, self.main__optimization_model_part.Nodes)

        else:
            raise Exception("Other input options are not yet implemented.")

        current_buffer_size = self.main__optimization_model_part.GetBufferSize()
        if(self.GetMinimumBufferSize() > current_buffer_size):
            self.main__optimization_model_part.SetBufferSize( self.GetMinimumBufferSize() )

        print ("Model reading finished.")

    def GetMinimumBufferSize(self):
        return 2
    
    # --------------------------------------------------------------------------
    def initializeProcesses( self, Model ):
        import process_factory
        #the process order of execution is important
        self.list_of_processes  = process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( self.analyzer_settings["optimization_settings"]["boundary_conditions_process_list"] )
        # if(ProjectParameters.Has("problem_process_list")):
        #     self.list_of_processes += process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( self.analyzer_settings["problem_process_list"] )
        # if(ProjectParameters.Has("output_process_list")):
        #     self.list_of_processes += process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( self.analyzer_settings["output_process_list"] )
                    
        #print list of constructed processes
        # if(echo_level>1):
        #     for process in self.list_of_processes:
        #         print(process)

        for process in self.list_of_processes:
            process.ExecuteInitialize()

    # --------------------------------------------------------------------------
    def initializeGIDOutput( self ):
        computing_model_part = self.main__optimization_model_part.GetSubModelPart(self.analyzer_settings["adjoint_fluid_analyzer"]["model_import_settings"]["volume_model_part_name"].GetString())
        problem_name = self.analyzer_settings["problem_data"]["problem_name"].GetString()

        from gid_output_process import GiDOutputProcess
        output_settings = self.analyzer_settings["output_configuration"]
        self.gid_output = GiDOutputProcess(computing_model_part, problem_name, output_settings)

        self.gid_output.ExecuteInitialize()

    # --------------------------------------------------------------------------
    def initializeSolutionLoop( self ):

        # print(self.adjoint.setup_parameters)

        self.design_surface = self.adjoint.setup_parameters["solver_settings"]["response_function_settings"]["sensitivity_model_part_name"]

        hdf5_output_process = {
                                    "kratos_module"    : "KratosMultiphysics.AdjointFluidApplication",
                                    "python_module"    : "output_primal_solution_process",
                                    "help"             : "",
                                    "process_name"     : "OutputPrimalSolutionProcess",
                                    "Parameters"          : {
                                        "model_part_name" : self.analyzer_settings["problem_data"]["model_part_name"].GetString(),
                                        "file_name"       : "primal_hdf5",
                                        "variable_list"   : ["VELOCITY", "ACCELERATION", "PRESSURE", "REACTION"],
                                        "alpha_bossak"    : self.adjoint.setup_parameters["solver_settings"]["scheme_settings"]["alpha_bossak"]
                                    }
                                }
        self.primal.AddCustomProcess(hdf5_output_process)

        hdf5_input_process = {
                                "kratos_module"    : "KratosMultiphysics.AdjointFluidApplication",
                                "python_module"    : "input_primal_solution_process",
                                "help"             : "",
                                "process_name"     : "InputPrimalSolutionProcess",
                                "Parameters"          : {
                                    "model_part_name" : self.analyzer_settings["problem_data"]["model_part_name"].GetString(),
                                    "file_name"       : "../primal/primal_hdf5",
                                    "variable_list"   : ["VELOCITY", "ACCELERATION", "PRESSURE", "REACTION"]
                                }
                            }
        
        self.adjoint.AddCustomProcess(hdf5_input_process, overwrite = True, existingProcessPythonModuleName = "input_primal_solution_process")

        hdf5_output_process = {
                                    "kratos_module"    : "KratosMultiphysics.AdjointFluidApplication",
                                    "python_module"    : "output_primal_solution_process",
                                    "help"             : "",
                                    "process_name"     : "OutputPrimalSolutionProcess",
                                    "Parameters"          : {
                                        "model_part_name" : self.analyzer_settings["problem_data"]["model_part_name"].GetString(),
                                        "file_name"       : "adjoint_hdf5",
                                        "variable_list"   : ["SHAPE_SENSITIVITY"],
                                        "alpha_bossak"    : self.adjoint.setup_parameters["solver_settings"]["scheme_settings"]["alpha_bossak"]
                                    }
                                }
        self.adjoint.AddCustomProcess(hdf5_output_process)

        self.mesh_solver.Initialize()
        # self.mesh_solver.SetEchoLevel(echo_level)

        # Start process
        for process in self.list_of_processes:
            process.ExecuteBeforeSolutionLoop()

        ## Set results when are written in a single file
        self.gid_output.ExecuteBeforeSolutionLoop()

    def initializeSolutionStep( self ):
        # processes to be executed at the begining of the solution step
        for process in self.list_of_processes:
            process.ExecuteInitializeSolutionStep()

        self.gid_output.ExecuteInitializeSolutionStep()
            
    def finalizeSolutionStep( self ):
        for process in self.list_of_processes:
            process.ExecuteFinalizeSolutionStep()
        
        self.gid_output.ExecuteFinalizeSolutionStep()

        # processes to be executed at the end of the solution step
        for process in self.list_of_processes:
            process.ExecuteFinalizeSolutionStep()

        # processes to be executed before witting the output
        for process in self.list_of_processes:
            process.ExecuteBeforeOutputStep()
        
        # write output results GiD: (frequency writing is controlled internally)
        if(self.gid_output.IsOutputStep()):
            self.gid_output.PrintOutput()
                        
        # processes to be executed after witting the output
        for process in self.list_of_processes:
            process.ExecuteAfterOutputStep()                  
            

    def importOptimizer( self, optimizer ):
        self.optimizer = optimizer

    # --------------------------------------------------------------------------
    def analyzeDesignAndReportToCommunicator( self, currentDesign, optimizationIteration, communicator ):

         # Calculation of objective function
        self.initializeNewSolutionStep( optimizationIteration )

        if optimizationIteration>1:
            print("\n> Starting ALEApplication to update the mesh")
            startTime = timer.time()
            self.updateMeshForAnalysis()
            print("> Time needed for updating the mesh = ",round(timer.time() - startTime,2),"s")
        
        self.initializeSolutionStep()

        current_path = "%s/%d" % (self.history_path, optimizationIteration)

        if not os.path.isdir(current_path):
            os.mkdir(current_path)

        print("\n> Starting primal problem solver")
        startTime = timer.time()
        self.primal.SetCurrentPath( "%s/primal" % current_path )
        self.primal.SetIterationPostFix( "%d_primal" % optimizationIteration )
        # for node in self.main__optimization_model_part.Nodes:
        #     print("Objective sensitivity: ", node.SolutionStepsDataHas(OBJECTIVE_SENSITIVITY))
        #     break
        self.primal.Execute()
        # for node in self.main__optimization_model_part.Nodes:
        #     print("Objective sensitivity: ", node.SolutionStepsDataHas(OBJECTIVE_SENSITIVITY))
        #     break
        print("> Time needed for solving the primal problem = ",round(timer.time() - startTime,2),"s")

        print("\n> Starting adjoint problem solver")
        startTime = timer.time()
        self.adjoint.SetCurrentPath( "%s/adjoint" % current_path )
        self.adjoint.SetIterationPostFix( "%d_adjoint" % optimizationIteration )
        # self.optimizer.AddVariables()
        # for node in self.main__optimization_model_part.Nodes:
        #     print("Objective sensitivity: ", node.SolutionStepsDataHas(OBJECTIVE_SENSITIVITY))
        #     break
        self.adjoint.Execute()
        # for node in self.main__optimization_model_part.Nodes:
        #     print("Objective sensitivity: ", node.SolutionStepsDataHas(OBJECTIVE_SENSITIVITY))
        #     break
        print("> Time needed for solving the adjoint problem = ",round(timer.time() - startTime,2),"s")

        print("\n> Starting calculation of response value")
        startTime = timer.time()                    
        response_value = self.calculateResponseValue()
        print("> Time needed for calculation of response value = ",round(timer.time() - startTime,2),"s")

        self.finalizeSolutionStep()

        communicator.reportFunctionValue("adjoint_fluid", response_value)  
        gradients = self.calculateGradientsOnDesignSurface("%s/adjoint" % current_path)
        communicator.reportGradient("adjoint_fluid", gradients) 

    # --------------------------------------------------------------------------
    def initializeNewSolutionStep( self, optimizationIteration ):
        self.main__optimization_model_part.CloneTimeStep( optimizationIteration )

    # --------------------------------------------------------------------------
    def updateMeshForAnalysis( self ):

        # Apply shape update as boundary condition for computation of mesh displacement 
        for node in self.main__optimization_model_part.GetSubModelPart(self.design_surface).Nodes:
            node.Fix(MESH_DISPLACEMENT_X)
            node.Fix(MESH_DISPLACEMENT_Y)
            node.Fix(MESH_DISPLACEMENT_Z)              
            disp = Vector(3)
            disp[0] = node.GetSolutionStepValue(SHAPE_UPDATE_X, 0)
            disp[1] = node.GetSolutionStepValue(SHAPE_UPDATE_Y, 0)
            disp[2] = node.GetSolutionStepValue(SHAPE_UPDATE_Z, 0)
            node.SetSolutionStepValue(MESH_DISPLACEMENT,0,disp)

        # Solve for mesh-update
        self.mesh_solver.Solve()

        # Update reference mesh (Since shape updates are imposed as incremental quantities)
        self.mesh_solver.get_mesh_motion_solver().UpdateReferenceMesh()

    # --------------------------------------------------------------------------
    def finalizeSolutionLoop( self ):
        for process in self.list_of_processes:
            process.ExecuteFinalize()
        self.gid_output.ExecuteFinalize()

    # --------------------------------------------------------------------------
    def calculateResponseValue( self ):
        objective_file = self.adjoint.setup_parameters["solver_settings"]["response_function_settings"]["output_file"]
        with open("%s/%s.data" % (self.adjoint.path, objective_file), "r") as file_input:
            lines = file_input.readlines()
        file_input.close()

        lines = lines[1:]

        if (lines[-1].strip()==""):
            lines = lines[:-1]

        if len(lines) > 1:
            components_0 = lines[0].strip().split()
            components_1 = lines[1].strip().split()
            dt = float(components_0[0]) - float(components_1[0])
            components_1 = lines[-1].strip().split()
            T = float(components_0[0]) - float(components_1[0])
        else:
            dt = 0.0
            T = 1.0

        response_value = 0.0
        for line in lines:
            components = line.strip().split()
            response_value += float(components[1])
        
        return response_value*dt/T
    
    # --------------------------------------------------------------------------
    def calculateGradientsOnDesignSurface( self, path ):
        design_surface_model_part = self.main__optimization_model_part.GetSubModelPart(self.design_surface)
        file_list = os.listdir( path )

        adjoint_sensitivity_files = []
        for file in file_list:
            if file[:12] == "adjoint_hdf5":
                adjoint_sensitivity_files.append( "%s/%s" % (path, file))

        designSurfaceNodeList = []
        for node in design_surface_model_part.Nodes:
            designSurfaceNodeList.append(node.Id)

        gradientOnDesignSurface = OptimizationUtilities(self.main__optimization_model_part, self.analyzer_settings).get_adjoint_design_surface_sensitivities (\
                                    adjoint_sensitivity_files,
                                    designSurfaceNodeList                               
                                    )
        return gradientOnDesignSurface

    # --------------------------------------------------------------------------




