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
            "material_import_settings" :{
                "materials_filename": "material_name.json"
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
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION) 
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

            #--------------------------Is this really needed again? the same import is also performed durning the primal solution!
            # Import constitutive laws.
            materials_imported = self.import_constitutive_laws()
            if materials_imported:
               print("    Constitutive law was successfully imported.")
            else:
                print("    Constitutive law was not imported.") 
            #----------------------------------------------------------------------------------------------------------

            # here we replace the dummy elements we read with proper elements
            self.settings.AddEmptyValue("element_replace_settings")
            if(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 3):
                self.settings["element_replace_settings"] = KratosMultiphysics.Parameters("""
                    {
                    "Add_string": "Adjoint",
                    "Add_before_in_element_name": "Element", 
                    "Add_before_in_condition_name": "Condition"
                    }
                    """)

            elif(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
                raise Exception("there is currently no 2D adjoint element")
            else:
                raise Exception("domain size is not 2 or 3")

            
            #KratosMultiphysics.ReplaceElementsAndConditionsProcess(self.main_model_part, self.settings["element_replace_settings"]).Execute()
            #print("I replace now elements!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            ShapeOptimizationApplication.ReplaceElementsAndConditionsForAdjointProblemProcess(self.main_model_part, self.settings["element_replace_settings"]).Execute()
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
        elif self.settings["response_function_settings"]["response_type"].GetString() == "rework_strain_energy":
            if (domain_size == 2):
                raise Exception("Currently only availible for 3D. Your choice is 2D")
            elif (domain_size == 3):
                self.response_function = ShapeOptimizationApplication.ReworkStrainEnergyResponseFunction(self.main_model_part, self.settings["response_function_settings"])
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
        if self.settings["response_function_settings"]["response_type"].GetString() == "rework_strain_energy":
            self._SolveSpecialStrainEnergy()
        else:
            self.solver.Solve()

    def SetEchoLevel(self, level):
        self.solver.SetEchoLevel(level)

    def Clear(self):
        self.solver.Clear()

    def Check(self):
        self.solver.Check()

    def _SolveSpecialStrainEnergy(self):

        self.response_function.Initialize()

        for node in self.main_model_part.Nodes:
            disp_x = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X,0)
            disp_y = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y,0)
            disp_z = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z,0)
            node.SetSolutionStepValue(KratosMultiphysics.ADJOINT_DISPLACEMENT_X,0,disp_x * 0.5)
            node.SetSolutionStepValue(KratosMultiphysics.ADJOINT_DISPLACEMENT_Y,0,disp_y * 0.5)
            node.SetSolutionStepValue(KratosMultiphysics.ADJOINT_DISPLACEMENT_Z,0,disp_z * 0.5)
            if self.settings["rotation_dofs"].GetBool():
                rot_x = node.GetSolutionStepValue(KratosMultiphysics.ROTATION_X,0)
                rot_y = node.GetSolutionStepValue(KratosMultiphysics.ROTATION_Y,0)
                rot_z = node.GetSolutionStepValue(KratosMultiphysics.ROTATION_Z,0)
                node.SetSolutionStepValue(KratosMultiphysics.ADJOINT_ROTATION_X,0,rot_x * 0.5)
                node.SetSolutionStepValue(KratosMultiphysics.ADJOINT_ROTATION_Y,0,rot_y * 0.5)
                node.SetSolutionStepValue(KratosMultiphysics.ADJOINT_ROTATION_Z,0,rot_z * 0.5)
                
        self.response_function.FinalizeSolutionStep()        


    def import_constitutive_laws(self): # copied from structural_mechanics_solver.py
        materials_filename = self.settings["material_import_settings"]["materials_filename"].GetString()
        if (materials_filename != ""):
            import read_materials_process
            # Create a dictionary of model parts.
            Model = {self.main_model_part.Name : self.main_model_part}
            for i in range(self.settings["problem_domain_sub_model_part_list"].size()):
                part_name = self.settings["problem_domain_sub_model_part_list"][i].GetString()
                Model.update({part_name: self.main_model_part.GetSubModelPart(part_name)})
            for i in range(self.settings["processes_sub_model_part_list"].size()):
                part_name = self.settings["processes_sub_model_part_list"][i].GetString()
                Model.update({part_name: self.main_model_part.GetSubModelPart(part_name)})
            # Add constitutive laws and material properties from json file to model parts.
            read_materials_process.ReadMaterialsProcess(Model, self.settings["material_import_settings"])
            materials_imported = True
        else:
            materials_imported = False
        return materials_imported     
            
