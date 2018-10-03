from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
KratosMultiphysics.CheckForPreviousImport()
from python_solver import PythonSolver

def CreateSolver(model, custom_settings):
    return LaplacianSolver(model, custom_settings["solver_settings"])

class LaplacianSolver(PythonSolver):
    def __init__(self, model, custom_settings):
        self.MoveMeshFlag = False

        #TODO: shall obtain the compute_model_part from the MODEL once the object is implemented
        
        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "model_part_name"        : "model",                   
             "domain_size"            : 2,
            "solver_type": "potential_flow_solver",
            "echo_level": 1,
            "relative_tolerance": 1e-5,
            "absolute_tolerance": 1e-9,
            "maximum_iterations": 1,
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
                    "smoother_type":"ilu0",
                    "coarsening_type":"ruge_stuben",
                    "coarse_enough" : 5000,
                    "krylov_type": "lgmres",
                    "tolerance": 1e-9,
                    "verbosity": 3,
                    "scaling": false
            }


        }""")

            # "linear_solver_settings"       : {
            #      "solver_type"     : "SkylineLUFactorizationSolver"
            #   }
         
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        model_part_name = self.settings["model_part_name"].GetString()
        super(LaplacianSolver,self).__init__(model, self.settings)

        if model_part_name == "":
            raise Exception('Please provide the model part name as the "model_part_name" (string) parameter!')

        if self.model.HasModelPart(model_part_name):
            self.main_model_part = self.model.GetModelPart(model_part_name)
        else:
            self.main_model_part = KratosMultiphysics.ModelPart(model_part_name)
        
        self.domain_size = custom_settings["domain_size"].GetInt()
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, self.domain_size)
                    
        #construct the linear solvers
        import linear_solver_factory
        self.linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])

        print("Construction of LaplacianSolver finished")

    def AddVariables(self):
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.POSITIVE_FACE_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NEGATIVE_FACE_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.CompressiblePotentialFlowApplication.VELOCITY_INFINITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_H)
        
    def AddDofs(self):
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.POSITIVE_FACE_PRESSURE, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.NEGATIVE_FACE_PRESSURE, self.main_model_part)
        
    def Initialize(self):
        time_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        move_mesh_flag = False #USER SHOULD NOT CHANGE THIS
        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(self.linear_solver)
        self.solver = KratosMultiphysics.ResidualBasedLinearStrategy(
            self.main_model_part, 
            time_scheme, 
            self.linear_solver,
            builder_and_solver,
            self.settings["compute_reactions"].GetBool(), 
            self.settings["reform_dofs_at_each_step"].GetBool(), 
            self.settings["calculate_solution_norm"].GetBool(), 
            move_mesh_flag)
        
        (self.solver).SetEchoLevel(self.settings["echo_level"].GetInt())
        self.solver.Check()
    def PrepareModelPart(self):
        if not self.model.HasModelPart(self.settings["model_part_name"].GetString()):
            self.model.AddModelPart(self.main_model_part)


    def ImportModelPart(self):
        
        if(self.settings["model_import_settings"]["input_type"].GetString() == "mdpa"):
            #here it would be the place to import restart data if required
            print(self.settings["model_import_settings"]["input_filename"].GetString())
            KratosMultiphysics.ModelPartIO(self.settings["model_import_settings"]["input_filename"].GetString()).ReadModelPart(self.main_model_part)
            
            throw_errors = False
            KratosMultiphysics.TetrahedralMeshOrientationCheck(self.main_model_part,throw_errors).Execute()
            #here we replace the dummy elements we read with proper elements
            self.settings.AddEmptyValue("element_replace_settings")
            if(self.domain_size == 3):
                self.settings["element_replace_settings"] = KratosMultiphysics.Parameters("""
                    {
                    "element_name":"CompressiblePotentialFlowElement3D4N",
                    "condition_name": "PotentialWallCondition3D3N"
                    }
                    """)
            elif(self.domain_size == 2):
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
        
    def SolveSolutionStep(self):
        (self.solver).Solve()
    
    #
    def SetEchoLevel(self, level):
        (self.solver).SetEchoLevel(level)

    #
    def Clear(self):
        (self.solver).Clear()

