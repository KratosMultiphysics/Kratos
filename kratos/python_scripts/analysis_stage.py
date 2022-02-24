# Importing Kratos
import KratosMultiphysics
from KratosMultiphysics.process_factory import KratosProcessFactory
from KratosMultiphysics.kratos_utilities import IssueDeprecationWarning
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

class AnalysisStage(object):
    """The base class for the AnalysisStage-classes in the applications
    Changes to this BaseClass have to be discussed first!
    """
    def __init__(self, model, project_parameters):
        """The constructor of the AnalysisStage-Object.

        It is intended to be called from the constructor
        of deriving classes:
        super(DerivedAnalysis, self).__init__(project_parameters)

        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- The Model to be used
        project_parameters -- The ProjectParameters used
        """
        if not isinstance(model, KratosMultiphysics.Model):
            raise Exception("Input is expected to be provided as a Kratos Model object")

        if not isinstance(project_parameters, KratosMultiphysics.Parameters):
            raise Exception("Input is expected to be provided as a Kratos Parameters object")

        self.model = model
        self.project_parameters = project_parameters

        ## Get echo level and parallel type
        self.echo_level = self.project_parameters["problem_data"]["echo_level"].GetInt()
        self.parallel_type = self.project_parameters["problem_data"]["parallel_type"].GetString()
        is_distributed_run = KratosMultiphysics.IsDistributedRun()

        if self.parallel_type == "OpenMP" and is_distributed_run:
            KratosMultiphysics.Logger.PrintWarning("Parallel Type", '"OpenMP" is specified as "parallel_type", but Kratos is running distributed!')
        if self.parallel_type == "MPI" and not is_distributed_run:
            KratosMultiphysics.Logger.PrintWarning("Parallel Type", '"MPI" is specified as "parallel_type", but Kratos is not running distributed!')

        self._GetSolver().AddVariables() # this creates the solver and adds the variables

    def Run(self):
        """This function executes the entire AnalysisStage
        It can be overridden by derived classes
        """
        self.Initialize()
        self.RunSolutionLoop()
        self.Finalize()

    def KeepAdvancingSolutionLoop(self):
        """This function specifies the stopping criteria for breaking the solution loop
        It can be overridden by derived classes
        """
        return self.time < self.end_time

    def RunSolutionLoop(self):
        """This function executes the solution loop of the AnalysisStage
        It can be overridden by derived classes
        """
        file1 = open('Results.txt', 'w')
        L = ["Displacement      Reaction\n"]
        file1.writelines(L)    
            
        file1.close()
        while self.KeepAdvancingSolutionLoop():
        
        # ##################################################
        
        #     # print("*********************************")
        #     # print("*********************************")
        #     # print("*********************************")
            displ = 0.0
            reac  = 0.0
             
            file1 = open('Results.txt', 'a') 
             
            for node in self._GetSolver().main_model_part.GetSubModelPart("DISPLACEMENT_Displacement_Auto2").Nodes:
        #         # print("Id-> ", node.Id)
        #         # print("Displ-> ", node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT))
        #         # print("Reaction-> ", node.GetSolutionStepValue(KratosMultiphysics.REACTION))
                displ = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
                reac += node.GetSolutionStepValue(KratosMultiphysics.REACTION_X)
               
            file1.writelines("%.5f      %.2f\n" % (displ,reac))
        #     # print("Displ-> ", displ)
        #     # print("Reaction-> ", reac)
            
                         
        #     # for elem in self._GetSolver().main_model_part.Elements:
        #     #         print(elem.Id)
        #     # print("*********************************")
        #     # print("*********************************")
        #     # print("*********************************")
                                    
            file1.close()
                       
        #     ##################################################


        

        ################################################## INICIO DE LA PRUEBA .mesh
       # file1 = open('Prueba.post.msh', 'w')

        # file1.close()
        # while self.KeepAdvancingSolutionLoop():
        #     coord_x = 0.0
        #     coord_y = 0.0
        #     coord_z = 0.0     
                
        #     file1 = open('Prueba.post.msh', 'a') 
        #     #file1.writelines("Coordinates\n")  
        
        #     if self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.STEP] == 0:

        #         dimension =  self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        #         Nnode = len (self._GetSolver().main_model_part.Nodes) 

                                              
        #         if self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
        #             file1.writelines("2D\n")
        #             # falta completar con los otros tipos de elementos: Point, Line, Triangle, Quadrilateral, Tetrahedra, Hexahedra, Prism, Pyramid, Sphere, Circle

        #         if self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 3:
            #         if Nnode ==4:
            #             file1.writelines("MESH \"Kratos_Hexahedra3D8_Mesh_1\" dimension %d ElemType Tetrahedra Nnode %.0d\n" % (dimension, Nnode))
            #         if Nnode ==8:
            #             file1.writelines("MESH \"Kratos_Hexahedra3D8_Mesh_1\" dimension %d ElemType Hexahedra Nnode %.0d\n" % (dimension, Nnode))
                
            #             #KRATOS_ERROR_IF(dimension != 2 && dimension !=3) << "Dimension has to be either 2 or 3! Current dimension: " << dimension << std::endl;
                
            #     file1.writelines("Coordinates\n")

            # for node in self._GetSolver().main_model_part.Nodes:
                
            #     id = node.Id
            #     coord_x = node.X
            #     coord_y = node.Y
            #     coord_z = node.Z
                
            #     if self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.STEP] == 0:
            #         file1.writelines("%.0f %.6f %.6f %.1f\n" % (id, coord_x, coord_y, coord_z))    
                                             
            # if self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.STEP] == 10:
            #     file1.writelines("End Coordinates\n") 
            #     file1.writelines("Elements\n")
            #     for elem in self._GetSolver().main_model_part.Elements:
            #         #node_1 = elem.GetGeometry()[0]
            #         elem_1_id = elem.Id
            #         node_1 = elem.GetGeometry()[0].Id
            #         node_2 = elem.GetGeometry()[1].Id
            #         node_3 = elem.GetGeometry()[2].Id
            #         node_4 = elem.GetGeometry()[3].Id
            #         node_5 = elem.GetGeometry()[4].Id
            #         node_6 = elem.GetGeometry()[5].Id
            #         node_7 = elem.GetGeometry()[6].Id
            #         node_8 = elem.GetGeometry()[7].Id
            #         color = 2
                    
            #     file1.writelines("%.0d %.0d %.0d %.0d %.0d %.0d %.0d %.0d %.0d %.0d\n" % (elem_1_id, node_1, node_2,  node_3, node_4, node_5, node_6, node_7, node_8, color))
            #     file1.writelines("End Elements\n")
                                    
             #   file1.close()
          
        ################################################## FIN DE LA PRUEBA .mesh

            ################################################## INICIO DE LA PRUEBA .res
                # file2 = open('Prueba.post.res', 'w')
                # file2.close()
            
                # file2 = open('Prueba.post.res', 'a')
                # file2.writelines("GiD Post Results File 1.0\n")
                # if Nnode == 4:
                #     file2.writelines("GaussPoints \"hex8_element_gp\" ElemType Tetrahedra\n")
                # if Nnode == 8:
                #     file2.writelines("GaussPoints \"hex8_element_gp\" ElemType Hexahedra\n")
                # file2.writelines("Number Of Gauss Points: %d\n" % (Nnode))
                # file2.writelines("Natural Coordinates: Internal\n")
                # file2.writelines("End GaussPoints\n")
                # time = self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME]
                # file2.writelines("Result \"DISPLACEMENT\" \"Kratos\" %.1f Vector OnNodes\n" % (time))
                # file2.writelines("Values\n")
                            
                # for node in self._GetSolver().main_model_part.Nodes:
                #     id = node.Id
                #     displ_x = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
                #     displ_y = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y)
                #     displ_z = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z)
                #     #file2.writelines("%d %f %f %f\n" % (id, displ_x, displ_y, displ_z))
                #     file2.writelines("%d %s %s %s\n" % (id, "{0:.4e}".format(displ_x), "{0:.4e}".format(displ_y), "{0:.4e}".format(displ_z)))
                            
                # file2.writelines("End Values\n")
                # file2.writelines("Result \"REACTION\" \"Kratos\" %.1f Vector OnNodes\n" % (time))
                # file2.writelines("Values\n")

                
                # for node in self._GetSolver().main_model_part.Nodes:
                #     id = node.Id
                #     reac_x = node.GetSolutionStepValue(KratosMultiphysics.REACTION_X)
                #     reac_y = node.GetSolutionStepValue(KratosMultiphysics.REACTION_Y)
                #     reac_z = node.GetSolutionStepValue(KratosMultiphysics.REACTION_Z)
                #     file2.writelines("%d %s %s %s \n" % (id, "{0:.4e}".format(reac_x), "{0:.4e}".format(reac_y), "{0:.4e}".format(reac_z)))
                   
                # file2.writelines("End Values\n")
                # #file2.writelines("Result \"VON_MISES_STRESS\" \"Kratos\" %.1f Scalar OnGaussPoints \"hex8_element_gp\"\n" % (time))
                # file2.writelines("Result \"VON_MISES_STRESS//a\" \"Kratos\" %.1f Scalar OnGaussPoints \"hex8_element_gp\"\n" % (time))
                # file2.writelines("Values\n")

                
                # for elem in self._GetSolver().main_model_part.Elements:
                #     elem_id=elem.Id
                #     Von_Mises = elem.CalculateOnIntegrationPoints(StructuralMechanicsApplication.VON_MISES_STRESS, self._GetSolver().GetComputingModelPart().ProcessInfo)
                #     Von_Mises_1 = Von_Mises [0]
                #     Von_Mises_2 = Von_Mises [1]
                #     Von_Mises_3 = Von_Mises [2]
                #     Von_Mises_4 = Von_Mises [3]
                #     Von_Mises_5 = Von_Mises [4]
                #     Von_Mises_6 = Von_Mises [5]
                #     Von_Mises_7 = Von_Mises [6]
                #     Von_Mises_8 = Von_Mises [7]
                    
                #     file2.writelines("%d %s\n %s\n %s\n %s\n %s\n %s\n %s\n %s \n" % (elem_id, "{0:.4e}".format(Von_Mises_1), "{0:.4e}".format(Von_Mises_2), "{0:.4e}".format(Von_Mises_3), "{0:.4e}".format(Von_Mises_4), "{0:.4e}".format(Von_Mises_5), "{0:.4e}".format(Von_Mises_6), "{0:.4e}".format(Von_Mises_7), "{0:.4e}".format(Von_Mises_8)))
                                                        
                # file2.writelines("End Values\n")
                
                # file2.close()

                

            ################################################## FIN DE LA PRUEBA .res

            self.time = self._GetSolver().AdvanceInTime(self.time)
            self.InitializeSolutionStep()
            self._GetSolver().Predict()
            is_converged = self._GetSolver().SolveSolutionStep()
            self.__CheckIfSolveSolutionStepReturnsAValue(is_converged)
            self.FinalizeSolutionStep()
            self.OutputSolutionStep()
        
    
    def Initialize(self):
        """This function initializes the AnalysisStage
        Usage: It is designed to be called ONCE, BEFORE the execution of the solution-loop
        This function has to be implemented in deriving classes!
        """
        # Modelers:
        self._CreateModelers()
        self._ModelersSetupGeometryModel()
        self._ModelersPrepareGeometryModel()
        self._ModelersSetupModelPart()

        self._GetSolver().ImportModelPart()
        self._GetSolver().PrepareModelPart()
        self._GetSolver().AddDofs()

        self.ModifyInitialProperties()
        self.ModifyInitialGeometry()

        ##here we initialize user-provided processes
        self.__CreateListOfProcesses() # has to be done after importing and preparing the ModelPart
        for process in self._GetListOfProcesses():
            process.ExecuteInitialize()

        self._GetSolver().Initialize()
        self.Check()

        self.ModifyAfterSolverInitialize()

        for process in self._GetListOfProcesses():
            process.ExecuteBeforeSolutionLoop()

        ## Stepping and time settings
        self.end_time = self.project_parameters["problem_data"]["end_time"].GetDouble()

        if self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
            self.time = self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME]
        else:
            self.time = self.project_parameters["problem_data"]["start_time"].GetDouble()

        ## If the echo level is high enough, print the complete list of settings used to run the simualtion
        if self.echo_level > 1:
            with open("ProjectParametersOutput.json", 'w') as parameter_output_file:
                parameter_output_file.write(self.project_parameters.PrettyPrintJsonString())

        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Analysis -START- ")

    def Finalize(self):
        """This function finalizes the AnalysisStage
        Usage: It is designed to be called ONCE, AFTER the execution of the solution-loop
        """
        for process in self._GetListOfProcesses():
            process.ExecuteFinalize()

        self._GetSolver().Finalize()

        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Analysis -END- ")

    def InitializeSolutionStep(self):
        """This function performs all the required operations that should be executed
        (for each step) BEFORE solving the solution step.
        """
        self.PrintAnalysisStageProgressInformation()

        self.ApplyBoundaryConditions() #here the processes are called
        self.ChangeMaterialProperties() #this is normally empty
        self._GetSolver().InitializeSolutionStep()


    def PrintAnalysisStageProgressInformation(self):
        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "STEP: ", self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.STEP])
        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "TIME: ", self.time)

    def FinalizeSolutionStep(self):
        """This function performs all the required operations that should be executed
        (for each step) AFTER solving the solution step.
        """
        self._GetSolver().FinalizeSolutionStep()

        for process in self._GetListOfProcesses():
            process.ExecuteFinalizeSolutionStep()

    def OutputSolutionStep(self):
        """This function printed / writes output files after the solution of a step
        """
        execute_was_called = False
        for output_process in self._GetListOfOutputProcesses():
            if output_process.IsOutputStep():
                if not execute_was_called:
                    for process in self._GetListOfProcesses():
                        process.ExecuteBeforeOutputStep()
                    execute_was_called = True

                output_process.PrintOutput()

        if execute_was_called:
            for process in self._GetListOfProcesses():
                process.ExecuteAfterOutputStep()

    def Check(self):
        """This function checks the AnalysisStage
        """
        self._GetSolver().Check()
        for process in self._GetListOfProcesses():
            process.Check()

    def Clear(self):
        """This function clears the AnalysisStage
        """
        self._GetSolver().Clear()

    def ModifyInitialProperties(self):
        """this is the place to eventually modify material properties in the stage """
        pass

    def ModifyInitialGeometry(self):
        """this is the place to eventually modify geometry (for example moving nodes) in the stage """
        pass

    def ModifyAfterSolverInitialize(self):
        """this is the place to eventually do any modification that requires the solver to be initialized """
        pass

    def ApplyBoundaryConditions(self):
        """here the boundary conditions is applied, by calling the InitializeSolutionStep function of the processes"""

        for process in self._GetListOfProcesses():
            process.ExecuteInitializeSolutionStep()

        #other operations as needed

    def ChangeMaterialProperties(self):
        """this function is where the user could change material parameters as a part of the solution step """
        pass

    def _GetSolver(self):
        if not hasattr(self, '_solver'):
            self._solver = self._CreateSolver()
        return self._solver

    def _CreateSolver(self):
        """Create the solver
        """
        raise Exception("Creation of the solver must be implemented in the derived class.")

    ### Modelers
    def _ModelersSetupGeometryModel(self):
        # Import or generate geometry models from external input.
        for modeler in self._GetListOfModelers():
            if self.echo_level > 1:
                KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Modeler: ", str(modeler), " Setup Geometry Model started.")
            modeler.SetupGeometryModel()
            if self.echo_level > 1:
                KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Modeler: ", str(modeler), " Setup Geometry Model finished.")

    def _ModelersPrepareGeometryModel(self):
        # Prepare or update the geometry model_part.
        for modeler in self._GetListOfModelers():
            if self.echo_level > 1:
                KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Modeler: ", str(modeler), " Prepare Geometry Model started.")
            modeler.PrepareGeometryModel()
            if self.echo_level > 1:
                KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Modeler: ", str(modeler), " Prepare Geometry Model finished.")

    def _ModelersSetupModelPart(self):
        # Convert the geometry model or import analysis suitable models.
        for modeler in self._GetListOfModelers():
            if self.echo_level > 1:
                KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Modeler: ", str(modeler), " Setup ModelPart started.")
            modeler.SetupModelPart()
            if self.echo_level > 1:
                KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Modeler: ", str(modeler), " Setup ModelPart finished.")

    ### Modelers
    def _GetListOfModelers(self):
        """ This function returns the list of modelers
        """
        if not hasattr(self, '_list_of_modelers'):
            raise Exception("The list of modelers was not yet created!")
        return self._list_of_modelers

    def _CreateModelers(self):
        """ List of modelers in following format:
        "modelers" : [{
            "modeler_name" : "geometry_import":
            "parameters" : {
                "echo_level" : 0:
                // settings for this modeler
            }
        },{ ... }]
        """
        self._list_of_modelers = []

        if self.project_parameters.Has("modelers"):
            from KratosMultiphysics.modeler_factory import KratosModelerFactory
            factory = KratosModelerFactory()

            modelers_list = self.project_parameters["modelers"]
            self._list_of_modelers = factory.ConstructListOfModelers(self.model, modelers_list)

    ### Processes
    def _GetListOfProcesses(self):
        """This function returns the list of processes involved in this Analysis
        """
        if not hasattr(self, '_list_of_processes'):
            raise Exception("The list of processes was not yet created!")
        return self._list_of_processes

    def _GetListOfOutputProcesses(self):
        """This function returns the list of output processes involved in this Analysis
        """
        if not hasattr(self, '_list_of_output_processes'):
            raise Exception("The list of output-processes was not yet created!")
        return self._list_of_output_processes

    def _CreateProcesses(self, parameter_name, initialization_order):
        """Create a list of processes
        Format:
        "processes" : {
            initial_processes : [
                { proces_specific_params },
                { proces_specific_params }
            ],
            boundary_processes : [
                { proces_specific_params },
                { proces_specific_params }
            ]
        }
        The order of intialization can be specified by setting it in "initialization_order"
        if e.g. the "boundary_processes" should be constructed before the "initial_processes", then
        initialization_order should be a list containing ["boundary_processes", "initial_processes"]
        see the functions _GetOrderOfProcessesInitialization and _GetOrderOfOutputProcessesInitialization
        """
        list_of_processes = []

        factory = KratosProcessFactory(self.model)

        if self.project_parameters.Has(parameter_name):
            processes_params = self.project_parameters[parameter_name]

            # first initialize the processes that depend on the order
            for processes_names in initialization_order:
                if processes_params.Has(processes_names):
                    list_of_processes += factory.ConstructListOfProcesses(processes_params[processes_names])

            # then initialize the processes that don't depend on the order
            for name, value in processes_params.items():
                if not name in initialization_order:
                    list_of_processes += factory.ConstructListOfProcesses(value) # Does this work? or should it be processes[name]

        return list_of_processes

    def _GetOrderOfProcessesInitialization(self):
        """This function can be overridden in derived classes if the order of
        initialization for the processes matters
        """
        return []

    def _GetOrderOfOutputProcessesInitialization(self):
        """This function can be overridden in derived classes if the order of
        initialization for the output-processes matters
        """
        return []

    def _CheckDeprecatedOutputProcesses(self, list_of_processes):
        deprecated_output_processes = []
        for process in list_of_processes:
            if issubclass(type(process), KratosMultiphysics.OutputProcess):
                deprecated_output_processes.append(process)
                msg  = "{} is an OutputProcess. However, it has been constructed as a regular process.\n"
                msg += "Please, define it as an 'output_processes' in the ProjectParameters."
                IssueDeprecationWarning("AnalysisStage", msg.format(process.__class__.__name__))
        return deprecated_output_processes

    def _GetSimulationName(self):
        """Returns the name of the Simulation
        """
        return "Analysis"

    def __CreateListOfProcesses(self):
        """This function creates the processes and the output-processes
        """
        order_processes_initialization = self._GetOrderOfProcessesInitialization()
        self._list_of_processes        = self._CreateProcesses("processes", order_processes_initialization)
        deprecated_output_processes    = self._CheckDeprecatedOutputProcesses(self._list_of_processes)
        order_processes_initialization = self._GetOrderOfOutputProcessesInitialization()
        self._list_of_output_processes = self._CreateProcesses("output_processes", order_processes_initialization)
        self._list_of_processes.extend(self._list_of_output_processes) # Adding the output processes to the regular processes
        self._list_of_output_processes.extend(deprecated_output_processes)

    def __CheckIfSolveSolutionStepReturnsAValue(self, is_converged):
        """In case the solver does not return the state of convergence
        (same as the SolvingStrategy does) then issue ONCE a deprecation-warning

        """
        if is_converged is None:
            if not hasattr(self, '_map_ret_val_depr_warnings'):
                self._map_ret_val_depr_warnings = []
            solver_class_name = self._GetSolver().__class__.__name__
            # used to only print the deprecation-warning once
            if not solver_class_name in self._map_ret_val_depr_warnings:
                self._map_ret_val_depr_warnings.append(solver_class_name)
                warn_msg  = 'Solver "{}" does not return '.format(solver_class_name)
                warn_msg += 'the state of convergence from "SolveSolutionStep"'
                IssueDeprecationWarning("AnalysisStage", warn_msg)
