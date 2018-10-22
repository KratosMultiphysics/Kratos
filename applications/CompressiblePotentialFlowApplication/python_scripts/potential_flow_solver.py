from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.ConvectionDiffusionApplication
import KratosMultiphysics.CompressiblePotentialFlowApplication
import eigen_solver_factory
KratosMultiphysics.CheckForPreviousImport()


def CreateSolver(main_model_part, custom_settings):
    return LaplacianSolver(main_model_part, custom_settings)

class LaplacianSolver:
    def __init__(self, model_part, custom_settings):
        self.MoveMeshFlag = False

        #TODO: shall obtain the compute_model_part from the MODEL once the object is implemented
        self.main_model_part = model_part    
        
        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "potential_flow_solver",
            "echo_level": 1,
            "relative_tolerance": 1e-5,
            "absolute_tolerance": 1e-9,
            "maximum_iterations": 10,
            "compute_reactions": false,
            "reform_dofs_at_each_step": false,
            "calculate_solution_norm" : false,
            "volume_model_part_name" : "volume_model_part",
            "skin_parts":[],
            "no_skin_parts"                : [],
            "model_import_settings": {
                    "input_type": "mdpa",
                    "input_filename": "unknown_name"
            },
            "linear_solver_settings": {
                    "solver_type": "AMGCL",
                    "max_iteration": 400,
                    "gmres_krylov_space_dimension": 100,
                    "coarse_enough" : 1000,
                    "smoother_type":"ilu0",
                    "coarsening_type":"ruge_stuben",
                    "krylov_type": "lgmres",
                    "tolerance": 1e-9,
                    "verbosity": 2,
                    "scaling": false
            }
        }""")
                    
        ##overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)
        
        #construct the linear solvers
        import linear_solver_factory
        self.linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])
        
        ##settings for the condition number
        settings_max = KratosMultiphysics.Parameters("""
        {
            "solver_type"             : "power_iteration_highest_eigenvalue_solver",
            "max_iteration"           : 10000,
            "tolerance"               : 1e-9,
            "required_eigen_number"   : 1,
            "verbosity"               : 0,
            "linear_solver_settings"  : {
                "solver_type"             : "SuperLUSolver",
                "max_iteration"           : 500,
                "tolerance"               : 1e-9,
                "scaling"                 : false,
                "verbosity"               : 0
            }
        }
        """)
        self.eigen_solver_max = eigen_solver_factory.ConstructSolver(settings_max)
        settings_min = KratosMultiphysics.Parameters("""
        {
            "solver_type"             : "power_iteration_eigenvalue_solver",
            "max_iteration"           : 10000,
            "tolerance"               : 1e-9,
            "required_eigen_number"   : 1,
            "verbosity"               : 0,
            "linear_solver_settings"  : {
                "solver_type"             : "SuperLUSolver",
                "max_iteration"           : 500,
                "tolerance"               : 1e-9,
                "scaling"                 : false,
                "verbosity"               : 0
            }
        }
        """)
        self.eigen_solver_min = eigen_solver_factory.ConstructSolver(settings_min)

        print("Construction of LaplacianSolver finished")

    def AddVariables(self):
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.POSITIVE_FACE_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NEGATIVE_FACE_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.WATER_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.CompressiblePotentialFlowApplication.UPPER_SURFACE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.CompressiblePotentialFlowApplication.LOWER_SURFACE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.CompressiblePotentialFlowApplication.WAKE_NORMAL)

        
        
    def AddDofs(self):
        for node in self.main_model_part.Nodes:
            # adding dofs
            node.AddDof(KratosMultiphysics.POSITIVE_FACE_PRESSURE)
            node.AddDof(KratosMultiphysics.NEGATIVE_FACE_PRESSURE)
        
    def Initialize(self):
        time_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        move_mesh_flag = False #USER SHOULD NOT CHANGE THIS
        
        #self.solver = KratosMultiphysics.ResidualBasedLinearStrategy(
            #self.main_model_part, 
            #time_scheme, 
            #self.linear_solver,
            #self.settings["compute_reactions"].GetBool(), 
            #self.settings["reform_dofs_at_each_step"].GetBool(), 
            #self.settings["calculate_solution_norm"].GetBool(), 
            #move_mesh_flag)
            
        conv_criteria = KratosMultiphysics.ResidualCriteria(
            self.settings["relative_tolerance"].GetDouble(), 
            self.settings["absolute_tolerance"].GetDouble())
        max_iterations = self.settings["maximum_iterations"].GetInt()
                
        self.solver = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(
            self.main_model_part, 
            time_scheme, 
            self.linear_solver,
            conv_criteria,
            max_iterations,
            self.settings["compute_reactions"].GetBool(), 
            self.settings["reform_dofs_at_each_step"].GetBool(), 
            move_mesh_flag)
        
        (self.solver).SetEchoLevel(self.settings["echo_level"].GetInt())
        self.solver.Check()
        
        #'''
        for node in self.main_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.POSITIVE_FACE_PRESSURE,0)
            node.SetSolutionStepValue(KratosMultiphysics.NEGATIVE_FACE_PRESSURE,0)
        #'''    
        
    def ImportModelPart(self):
        
        if(self.settings["model_import_settings"]["input_type"].GetString() == "mdpa"):
            #here it would be the place to import restart data if required
            print(self.settings["model_import_settings"]["input_filename"].GetString())
            KratosMultiphysics.ModelPartIO(self.settings["model_import_settings"]["input_filename"].GetString()).ReadModelPart(self.main_model_part)
                     
            throw_errors = False
            KratosMultiphysics.TetrahedralMeshOrientationCheck(self.main_model_part,throw_errors).Execute()
            #here we replace the dummy elements we read with proper elements
            self.settings.AddEmptyValue("element_replace_settings")
            if(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 3):
                self.settings["element_replace_settings"] = KratosMultiphysics.Parameters("""
                    {
                    "element_name":"CompressiblePotentialFlowElement3D4N",
                    "condition_name": "PotentialWallCondition3D3N"
                    }
                    """)
            elif(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
                self.settings["element_replace_settings"] = KratosMultiphysics.Parameters("""
                    {
                    "element_name":"CompressiblePotentialFlowElement2D3N",
                    "condition_name": "PotentialWallCondition2D2N"
                    }
                    """)
            else:
                raise Exception("Domain size is not 2 or 3!!")
            
            KratosMultiphysics.ReplaceElementsAndConditionsProcess(self.main_model_part, self.settings["element_replace_settings"]).Execute()
            
        else:
            raise Exception("other input options are not yet implemented")
        
        current_buffer_size = self.main_model_part.GetBufferSize()
        if(self.GetMinimumBufferSize() > current_buffer_size):
            self.main_model_part.SetBufferSize( self.GetMinimumBufferSize() )
                
        print ("model reading finished")
        
    def GetMinimumBufferSize(self):
        return 2;
    
    def GetComputingModelPart(self):
        return self.main_model_part
        
    def GetOutputVariables(self):
        pass
        
    def ComputeDeltaTime(self):
        pass
        
    def SaveRestart(self):
        pass #one should write the restart file here
        
    def Solve(self):
        (self.solver).Solve()
        #'''
        # Solve condition number
        condition_number_utility = KratosMultiphysics.ConditionNumberUtility()
        condition_number = condition_number_utility.GetConditionNumber(self.solver.GetSystemMatrix(), self.eigen_solver_max, self.eigen_solver_min)
        print('Condition number = {:.6e}'.format(condition_number))
        condition_number_file = open("condition_number.txt", "a")
        condition_number_file.write('{0:.6e}\t'.format(condition_number))
        #'''

    #
    def SetEchoLevel(self, level):
        (self.solver).SetEchoLevel(level)

    #
    def Clear(self):
        (self.solver).Clear()

