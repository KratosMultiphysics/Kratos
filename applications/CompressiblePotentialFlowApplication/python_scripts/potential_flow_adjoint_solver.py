from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
KratosMultiphysics.CheckForPreviousImport()
from python_solver import PythonSolver
from potential_flow_solver import PotentialSolver

def CreateSolver(model, custom_settings):
    return PotentialAdjointSolver(model, custom_settings["solver_settings"])

class PotentialAdjointSolver(PotentialSolver):
    def __init__(self, model, custom_settings):
        adjoint_settings = KratosMultiphysics.Parameters("""
        {
            "scheme_settings" : {
                "scheme_type": "adjoint_potential"
            },
            "element_replace_settings" : {
                "element_name":"IncompressibleAdjointPotentialFlowElement",
                "condition_name": "IncompressibleAdjointPotentialWallCondition"
            }
        }
        """)

        self.validate_and_transfer_matching_settings(custom_settings, adjoint_settings)
        self.scheme_settings = adjoint_settings["scheme_settings"]
        self.element_replace_settings=adjoint_settings["element_replace_settings"]

        self.response_function_settings = custom_settings["response_function_settings"].Clone()
        self.sensitivity_settings = custom_settings["sensitivity_settings"].Clone()
        custom_settings.RemoveValue("response_function_settings")
        custom_settings.RemoveValue("sensitivity_settings")
        # Construct the base solver.
        super(PotentialAdjointSolver, self).__init__(model, custom_settings)
        self.print_on_rank_zero("::[PotentialAdjointSolver]:: ", "Construction finished")

    def AddVariables(self):
        super(PotentialAdjointSolver, self).AddVariables()
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.CompressiblePotentialFlowApplication.ADJOINT_POSITIVE_POTENTIAL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.CompressiblePotentialFlowApplication.ADJOINT_NEGATIVE_POTENTIAL)
        # self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SHAPE_SENSITIVITY)
        self.print_on_rank_zero("::[PotentialAdjointSolver]:: ", "Variables ADDED")
        
    def AddDofs(self):
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.CompressiblePotentialFlowApplication.ADJOINT_POSITIVE_POTENTIAL, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.CompressiblePotentialFlowApplication.ADJOINT_NEGATIVE_POTENTIAL, self.main_model_part)
        
    def Initialize(self):
        time_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        move_mesh_flag = False #USER SHOULD NOT CHANGE THIS

        if self.settings["problem_type"].GetString() == "compressible":
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
        else:
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
        super(PotentialAdjointSolver, self).PrepareModelPart()
        # TODO Why does replacement need to happen after reading materials?
        KratosMultiphysics.ReplaceElementsAndConditionsProcess(self.main_model_part,self.element_replace_settings).Execute()
        self.print_on_rank_zero("::[PotentialAdjointSolver]:: ", "ModelPart prepared for Solver.")

    def ImportModelPart(self):

        if(self.settings["model_import_settings"]["input_type"].GetString() == "mdpa"):
            #here it would be the place to import restart data if required
            print(self.settings["model_import_settings"]["input_filename"].GetString())
            IOdir=self.settings["model_import_settings"]["input_filename"].GetString()
            KratosMultiphysics.ModelPartIO(IOdir).ReadModelPart(self.main_model_part)
            
            throw_errors = False
            KratosMultiphysics.TetrahedralMeshOrientationCheck(self.main_model_part,throw_errors).Execute()
            #here we replace the dummy elements we read with proper elements
            self.settings.AddEmptyValue("element_replace_settings")
            if (self.settings["problem_type"].GetString() == "incompressible"):
                if(self.domain_size == 3):
                    self.settings["element_replace_settings"] = KratosMultiphysics.Parameters("""
                        {
                        "element_name":"IncompressiblePotentialFlowElement3D4N",
                        "condition_name": "IncompressiblePotentialWallCondition3D3N"
                        }
                        """)
                elif(self.domain_size == 2):
                    self.settings["element_replace_settings"] = KratosMultiphysics.Parameters("""
                        {
                        "element_name":"IncompressiblePotentialFlowElement2D3N",
                        "condition_name": "IncompressiblePotentialWallCondition2D2N"
                        }
                        """)
                else:
                    raise Exception("Domain size is not 2 or 3!!")
            elif (self.settings["problem_type"].GetString() == "compressible"):
                if(self.domain_size == 3):
                    self.settings["element_replace_settings"] = KratosMultiphysics.Parameters("""
                        {
                        "element_name":"CompressiblePotentialFlowElement3D4N",
                        "condition_name": "CompressiblePotentialWallCondition3D3N"
                        }
                        """)
                elif(self.domain_size == 2):
                    self.settings["element_replace_settings"] = KratosMultiphysics.Parameters("""
                        {
                        "element_name":"CompressiblePotentialFlowElement2D3N",
                        "condition_name": "CompressiblePotentialWallCondition2D2N"
                        }
                        """)
                else:
                    raise Exception("Domain size is not 2 or 3!!")
            elif (self.settings["problem_type"].GetString() == "incompressible_stresses"):
                if(self.domain_size == 2):
                    self.settings["element_replace_settings"] = KratosMultiphysics.Parameters("""
                        {
                        "element_name":"IncompressibleStressesPotentialFlowElement2D3N",
                        "condition_name": "IncompressibleStressesPotentialWallCondition2D2N"
                        }
                        """)
                else:
                    raise Exception("Domain size is not 2!!")
            elif (self.settings["problem_type"].GetString() == "incompressible_alpha"):
                if(self.domain_size == 2):
                    self.settings["element_replace_settings"] = KratosMultiphysics.Parameters("""
                        {
                        "element_name":"IncompressibleAlphaPotentialFlowElement2D3N",
                        "condition_name": "IncompressibleStressesPotentialWallCondition2D2N"
                        }
                        """)
                else:
                    raise Exception("Domain size is not 2!!")
            elif (self.settings["problem_type"].GetString() == "incompressible_alpha_full"):
                if(self.domain_size == 2):
                    self.settings["element_replace_settings"] = KratosMultiphysics.Parameters("""
                        {
                        "element_name":"IncompressibleAlphaFullPotentialFlowElement2D3N",
                        "condition_name": "IncompressiblePotentialWallCondition2D2N"
                        }
                        """)
                else:
                    raise Exception("Domain size is not 2!!")
            elif (self.settings["problem_type"].GetString() == "incompressible_stresses_mix"):
                if(self.domain_size == 2):
                    self.settings["element_replace_settings"] = KratosMultiphysics.Parameters("""
                        {
                        "element_name":"IncompressibleStressesMixPotentialFlowElement2D3N",
                        "condition_name": "IncompressibleStressesPotentialWallCondition2D2N"
                        }
                        """)
                else:
                    raise Exception("Domain size is not 2!!")
            else:
                raise Exception("Problem type not defined!!")
            
            KratosMultiphysics.ReplaceElementsAndConditionsProcess(self.main_model_part, self.settings["element_replace_settings"]).Execute()

        else:
            raise Exception("other input options are not yet implemented")
        print("Solving",self.settings["problem_type"].GetString() ,"case")
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
    
    def AdvanceInTime(self, current_time):
        dt = 1 #self._ComputeDeltaTime()
        new_time = current_time + dt

        # self.main_model_part.CloneTimeStep(new_time)
        self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] += 1

        return new_time

    def InitializeSolutionStep(self):        
        self.solver.InitializeSolutionStep()


    def SolveSolutionStep(self):
        (self.solver).Solve() 

    def FinalizeSolutionStep(self):        
        self.solver.FinalizeSolutionStep()

    def Predict(self):
        self.solver.Predict()

    #
    def SetEchoLevel(self, level):
        (self.solver).SetEchoLevel(level)

    #
    def Clear(self):
        (self.solver).Clear()

