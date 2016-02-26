from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics as kratoscore
import KratosMultiphysics.SolidMechanicsApplication as solid_application
import KratosMultiphysics.StructuralMechanicsApplication as structural_application

# Check that KratosMultiphysics was imported in the main script
kratoscore.CheckForPreviousImport()



def AddVariables(model_part, settings=default_settings):
    model_part.AddNodalSolutionStepVariable(kratoscore.DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(kratoscore.REACTION)
    
    if(setting.RotationDofs == True):
        model_part.AddNodalSolutionStepVariable(kratoscore.ROTATION)
        model_part.AddNodalSolutionStepVariable(kratoscore.TORQUE)
 
    if(setting.PressureDofs == True):
        model_part.AddNodalSolutionStepVariable(kratoscore.PRESSURE)
        model_part.AddNodalSolutionStepVariable(kratoscore.PRESSURE_REACTION)

    print("variables for the  Structural Mechanics Static solver symbolic added correctly")


def AddDofs(model_part, settings=default_settings):
    for node in model_part.Nodes:
        # adding dofs
        node.AddDof(kratoscore.DISPLACEMENT_X, kratoscore.REACTION_X)
        node.AddDof(kratoscore.DISPLACEMENT_Y, kratoscore.REACTION_Y)
        node.AddDof(kratoscore.DISPLACEMENT_Z, kratoscore.REACTION_Z)
        
        if(setting.RotationDofs == True):
            node.AddDof(kratoscore.ROTATION_X, kratoscore.TORQUE_X)
            node.AddDof(kratoscore.ROTATION_Y, kratoscore.TORQUE_Y)
            node.AddDof(kratoscore.ROTATION_Z, kratoscore.TORQUE_Z)
            
        if(setting.PressureDofs == True):
            node.AddDof(kratoscore.PRESSURE, kratoscore.PRESSURE_REACTION)

    print("dofs for the  Structural Mechanics Static solver added correctly")

def CreateSolver(model_part, settings=default_settings):
    static_solver = StaticSolver(model_part, settings) 
    return static_solver

class StaticSolver:
    
    ##here set the defaults
    self.settings = {
        "solver_type" : "static_structural_solver",
        "domain_size" :  3,
        "echo_level" :  1,
        "RotationDofs" :  False,          
        "PressureDofs" :  False,          
        "problem_is_linear" : = False,
        "displacement_relative_tolerance" :  1e-4,
        "displacement_relative_tolerance" :  1e-9,
        "residual_relative_tolerance" :  1e-4,
        "residual_absolute_tolerance" :  1e-9,
        "max_iteration" :  30,
        "linear_solver_settings" : 
            {
            "solver_type" : "Super LU" 
            "max_iteration" : 500
            "tolerance" : 1e-9
            "scaling" : False  #false->False
            "verbosity" : 1 
            }
        }  
            
    
    #TODO: this shall be in the base class
    def apply_settings(self,custom_settings,settings):
        
        unexpected_parameter = []
        for param in self.custom_settings:
            if self.settings.has(param.Name()) == False:
                unexpected_parameter.append(param )
            else:
                self.settings[param.Name()] = param.Value()
                
        if(len(unexpected_parameter) != 0):
            for wrong_param in unexpected_parameter:
                print ("unexpected parameter : ", wrong_param.Name(), " with value ", wrong_param.Value())
            
            print("   setting currently employed : "
            print(self.settings)
            
            
    
    def __init__(self, model_part, custom_settings=None): 
        
        if(custom_settings != None):
            apply_settings(custom_settings, self.settings)
        
        #TODO: shall obtain the model_part from the MODEL once the object is implemented
        self.model_part = model_part
        
        #construct the linear solver
        import linear_solver_factory
        linear_solver = linear_solver_factory.ConstructSolver(settings.linear_solver_settings)
                
        time_scheme = kratoscore.ResidualBasedIncrementalUpdateStaticScheme() 
        
        builder_and_solver = kratoscore.ResidualBasedBlockBuilderAndSolver(self.linear_solver)
        
        move_mesh_flag = False #user should NOT configure this
        
        if(problem_is_linear == True):
            self.solver = kratoscore.ResidualBasedLinearStrategy(self.model_part, 
                                                                 time_scheme, 
                                                                 linear_solver, 
                                                                 builder_and_solver, 
                                                                 self.settings["compute_reactions"], 
                                                                 self.settings["reform_dofs_at_each_step"], 
                                                                 move_mesh_flag)
        else: #nonlinear case
            
            #TODO: we shall construct a factory for the convergence criteria, like what is done for the linear solver
            conv_criteria = kratoscore.DisplacementCriteria(
                                                                 self.settings["displacement_relative_tolerance"],
                                                                 self.settings["displacement_absoulute_tolerance"]
                                                                 )
                    
            self.solver = kratoscore.ResidualBasedNewtonRaphsonStrategy(self.model_part, 
                                                                 time_scheme, 
                                                                 linear_solver, 
                                                                 conv_criteria,
                                                                 builder_and_solver, 
                                                                 self.settings["max_iter"], 
                                                                 self.settings["compute_reactions"], 
                                                                 self.settings["reform_dofs_at_each_step"], 
                                                                 move_mesh_flag)
            
        (self.solver).SetEchoLevel(self.settings["echo_level"])
        self.solver.Check()

        print("Construction Static Solver finished")
        
    def GetMinimumBufferSize(self):
        return 1;

    def ImportModelPart(self):
        #here it would be the place to import restart data if required
        ModelPartIO(self.settings.input_filename).ReadModelPart(self.model_part)
        
        # set the constitutive law
        import constitutive_law_python_utility as constitutive_law_utils
        constitutive_law = constitutive_law_utils.ConstitutiveLawUtility(model_part, self.settings["domain_size"]);
        constitutive_law.Initialize();
        
        print ("model reading finished")
        
    def GetOutputVariables(self):
        pass
        
    def ComputeDeltaTime(self):
        pass
        
    def SaveRestart(self):
        pass #one should write the restart file here
        
    def Initialize(self):
        print ("Initialization stokes solver finished")
        
    def (self):
    
    def Solve(self):
        self.solver.Solve()

    def SetEchoLevel(self, level):
        self.solver.SetEchoLevel(level)

    def Clear(self):
        self.solver.Clear()
        
    def Check(self):
        self.solver.Check()


