from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.ShapeOptimizationApplication as ShapeOptimizationApplication
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication # needed for POINT_LOAD variable

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

def CreateSolver(main_model_part, custom_settings):
    return AdjointStructuralSolver(main_model_part, custom_settings)

class AdjointStructuralSolver:

    def __init__(self, main_model_part, custom_settings):

        self.main_model_part = main_model_part

        # default settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type" : "adjoint_structural_solver",
            "scheme_settings" : {
                "scheme_type" : "structural"
            },
          
            "response_function_settings" : {
                "response_type" : "unknown_name"
            },
            "model_import_settings" : {
                "input_type"     : "mdpa",
                "input_filename" : "unknown_name"
            },
            "linear_solver_settings" : {
                "solver_type" : "AMGCL"
            },
            "bodies_list": [],
            "problem_domain_sub_model_part_list": ["solid"],
            "processes_sub_model_part_list": [""],
            "computing_model_part_name" : "computing_domain",
            "rotation_dofs": false,
            "echo_level"  : 0
        }""")

        # overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        # construct the linear solver
        import linear_solver_factory
        self.linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])

        print("Construction of AdjointStructuralSolver finished")

    def GetMinimumBufferSize(self):
        return 2

    def AddVariables(self):
        
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)  
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ADJOINT_DISPLACEMENT) 
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ADJOINT_ROTATION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SHAPE_SENSITIVITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.POINT_LOAD_SENSITIVITY)
        self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_LOAD) #Is there a better solution for Neumann BCs? 

        #---> is reaction and torque necessary????????????
        #self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ADJOINT_REACTION)
        #self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ADJOINT_TORQUE)
        
        print("Variables for the adjoint structural solver added correctly")

    def ImportModelPart(self):

        if(self.settings["model_import_settings"]["input_type"].GetString() == "mdpa"):
            # here it would be the place to import restart data if required
            KratosMultiphysics.ModelPartIO(self.settings["model_import_settings"]["input_filename"].GetString()).ReadModelPart(self.main_model_part)

            #--------------------------> which aux_params are really necessary? <------------------------------------
            # here we shall check that the input read has the shape we like
            aux_params = KratosMultiphysics.Parameters("{}")
            aux_params.AddValue("computing_model_part_name",self.settings["computing_model_part_name"])
            aux_params.AddValue("problem_domain_sub_model_part_list",self.settings["problem_domain_sub_model_part_list"])
            aux_params.AddValue("processes_sub_model_part_list",self.settings["processes_sub_model_part_list"])
            if( self.settings.Has("bodies_list") ):
                aux_params.AddValue("bodies_list",self.settings["bodies_list"])
            #aux_params.AddValue("problem_model_part_name",self.settings["problem_model_part_name"])
            #aux_params.AddValue("skin_parts",self.settings["skin_parts"])    
           

            # here we replace the dummy elements we read with proper elements
            self.settings.AddEmptyValue("element_replace_settings")
            if(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 3):
                self.settings["element_replace_settings"] = KratosMultiphysics.Parameters("""
                    {
                    "element_name": "CrBeamElement3D2NForSA",
                    "condition_name": "PointLoadCondition3D1NForSA" 
                    }
                    """)
            # ---------------> the condition here is used as dummy. In this case only the 
            # elements have to be replaced
            elif(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
                raise Exception("there is currently no 2D adjoint element")
            else:
                raise Exception("domain size is not 2 or 3")

            # TODO: the replacement function has to be reworked. Since only the elements and 
            # not the conditions have to be replaced. Also it has to condsidered that a model can consist 
            # out of more than one element type (e.g. model with beam shell elements). The current 
            # replacement function replaced all elements with the same element which is chosen 
            # by "element_name": some_element (see above)
            KratosMultiphysics.ReplaceElementsAndConditionsProcess(self.main_model_part, self.settings["element_replace_settings"]).Execute()

            import check_and_prepare_model_process_structural
            check_and_prepare_model_process_structural.CheckAndPrepareModelProcess(self.main_model_part, aux_params).Execute()

        else:
            raise Exception("Other input options are not yet implemented.")

        current_buffer_size = self.main_model_part.GetBufferSize()
        if(self.GetMinimumBufferSize() > current_buffer_size):
            self.main_model_part.SetBufferSize( self.GetMinimumBufferSize() )

        print ("Model reading finished.")


    def AddDofs(self):

        for node in self.main_model_part.Nodes:
            # adding dofs
            node.AddDof(KratosMultiphysics.ADJOINT_DISPLACEMENT_X)
            node.AddDof(KratosMultiphysics.ADJOINT_DISPLACEMENT_Y)
            node.AddDof(KratosMultiphysics.ADJOINT_DISPLACEMENT_Z)
            if self.settings["rotation_dofs"].GetBool():
                node.AddDof(KratosMultiphysics.ADJOINT_ROTATION_X)
                node.AddDof(KratosMultiphysics.ADJOINT_ROTATION_Y)
                node.AddDof(KratosMultiphysics.ADJOINT_ROTATION_Z)

        #---> is reaction and torque necessary????????????        

        print("DOFs for the structural adjoint solver added correctly.")


    def Initialize(self):

        self.computing_model_part = self.GetComputingModelPart()

        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        if self.settings["response_function_settings"]["response_type"].GetString() == "local_stress":
            if (domain_size == 2):
                raise Exception("Currently only availible for 3D. Your choice is 2D")
            elif (domain_size == 3):
                self.response_function = ShapeOptimizationApplication.LocalStressResponseFunction(self.main_model_part, self.settings["response_function_settings"])
            else:
                raise Exception("Invalid DOMAIN_SIZE: " + str(domain_size))
        elif self.settings["response_function_settings"]["response_type"].GetString() == "nodal_displacement":
            if (domain_size == 2):
                raise Exception("Currently only availible for 3D. Your choice is 2D")
            elif (domain_size == 3):
                self.response_function = ShapeOptimizationApplication.NodalDisplacementResponseFunction(self.main_model_part, self.settings["response_function_settings"])
            else:
                raise Exception("Invalid DOMAIN_SIZE: " + str(domain_size))        
        else:
            raise Exception("invalid response_type: " + self.settings["response_function_settings"]["response_type"].GetString())

        # -------------------------------->   add here all structural response functions   <-------------------------------------------


        if self.settings["scheme_settings"]["scheme_type"].GetString() == "structural":
            self.time_scheme = ShapeOptimizationApplication.AdjointStructuralScheme(self.settings["scheme_settings"], self.response_function)
        else:
            raise Exception("invalid scheme_type: " + self.settings["scheme_settings"]["scheme_type"].GetString())

        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(self.linear_solver)

        self.solver = KratosMultiphysics.ResidualBasedLinearStrategy(self.main_model_part,
                                                                     self.time_scheme,
                                                                     self.linear_solver,
                                                                     builder_and_solver,
                                                                     False,
                                                                     False,
                                                                     False,
                                                                     False)

        (self.solver).SetEchoLevel(self.settings["echo_level"].GetInt())

        self.solver.Check()
        
        print ("Adjoint solver initialization finished.")

    def GetComputingModelPart(self):
        # get the submodelpart generated in CheckAndPrepareModelProcess
        return self.main_model_part.GetSubModelPart(self.settings["computing_model_part_name"].GetString())

    def GetOutputVariables(self):
        pass

    def ComputeDeltaTime(self):
        pass

    def SaveRestart(self):
        pass #one should write the restart file here

    def DivergenceClearance(self):
        pass
        
    def SolverInitialize(self):
        self.solver.Initialize()

    def SolverInitializeSolutionStep(self):
        self.solver.InitializeSolutionStep()

    def SolverPredict(self):
        self.solver.Predict()

    def SolverSolveSolutionStep(self):
        self.solver.SolveSolutionStep()

    def SolverFinalizeSolutionStep(self):
        self.solver.FinalizeSolutionStep()

    def Solve(self):
        self.solver.Solve()

    def SetEchoLevel(self, level):
        self.solver.SetEchoLevel(level)

    def Clear(self):
        self.solver.Clear()

    def Check(self):
        self.solver.Check()
