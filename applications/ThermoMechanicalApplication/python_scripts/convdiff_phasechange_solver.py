from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.ThermoMechanicalApplication import *
from KratosMultiphysics.IncompressibleFluidApplication import *
#from KratosMultiphysics.ExternalSolversApplication import *
import KratosMultiphysics.ExternalSolversApplication as ExternalSolver

#from KratosMultiphysics.StructuralApplication import *
#from KratosMultiphysics.SolidMechanicsApplication import *
import KratosMultiphysics.SolidMechanicsApplication as SolidMechanics

# Check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()


def AddVariables(model_part, settings):
    #my_settings.SetDensityVariable(DENSITY)
    #my_settings.SetDiffusionVariable(CONDUCTIVITY)
    #my_settings.SetUnknownVariable(TEMPERATURE)
    #my_settings.SetVolumeSourceVariable(HEAT_FLUX)
    #my_settings.SetSurfaceSourceVariable(FACE_HEAT_FLUX)
    #my_settings.SetMeshVelocityVariable(MESH_VELOCITY)
    #my_settings.SetConvectionVariable(VELOCITY)
    #my_settings.SetTransferCoefficientVariable(HTC)
    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part.AddNodalSolutionStepVariable(NODAL_AREA);
    model_part.AddNodalSolutionStepVariable(settings.GetMeshVelocityVariable());
    model_part.AddNodalSolutionStepVariable(settings.GetUnknownVariable());
    model_part.AddNodalSolutionStepVariable(SPECIFIC_HEAT);
    model_part.AddNodalSolutionStepVariable(settings.GetVolumeSourceVariable());
    model_part.AddNodalSolutionStepVariable(settings.GetDensityVariable());
    model_part.AddNodalSolutionStepVariable(settings.GetDiffusionVariable());
    model_part.AddNodalSolutionStepVariable(settings.GetSurfaceSourceVariable());
    #model_part.AddNodalSolutionStepVariable(settings.GetTransferCoefficientVariable());
    #model_part.AddNodalSolutionStepVariable(HTC); #Mirar para quitarla de enmedio
    model_part.AddNodalSolutionStepVariable(ENTHALPY);
    model_part.AddNodalSolutionStepVariable(SOLIDFRACTION);
    model_part.AddNodalSolutionStepVariable(SOLIDFRACTION_RATE);
    model_part.AddNodalSolutionStepVariable(DISTANCE);


    #TODO: remove this variables!!
    model_part.AddNodalSolutionStepVariable(NODAL_PAUX);
    #model_part.AddNodalSolutionStepVariable(NODAL_VOLUME);
    #model_part.AddNodalSolutionStepVariable(SOLIDIF_TIME);
    #model_part.AddNodalSolutionStepVariable(SOLIDIF_MODULUS);
    model_part.AddNodalSolutionStepVariable(IS_VISITED);
    #model_part.AddNodalSolutionStepVariable(MOLD_INNER_TEMPERATURE);

    print("variables for the THERMAL_SOLVER added correctly")


def AddDofs(model_part, settings):
    for node in model_part.Nodes:
        node.AddDof(settings.GetUnknownVariable());

    print("dofs for the THERMAL_SOLVER added correctly")


class Solver:
    #

    def __init__(self, model_part, domain_size, my_settings, SolverSettings):
        Timer.Start("ConvDiffPhaseChange_Constructor")
        self.model_part = model_part
        self.settings = my_settings
        self.domain_size = domain_size
        self.max_distance = 1e6

        self.echo_level = 0
        self.CalculateReactionFlag = False
        self.ReformDofSetAtEachStep = False
        self.CalculateNormDxFlag = True
        self.MoveMeshFlag = False

        #functions to skip one of the stages
        self.skip_stage0 = SolverSettings.SkipStage0
        self.skip_stage1 = SolverSettings.SkipStage1
        self.use_linesearch_in_step1 = SolverSettings.UseLineSearch

        ##utility to effective store a copy of the variables
        self.variable_utils = VariableUtils()
        self.activation_utils = ActivationUtilities()


        ##strategy to be used in step0 - pure convection step - a non-symmetric linear solver is required
        gmres_size = SolverSettings.LinearSolverGmresSize
        tol = SolverSettings.LinearSolverTolerance#BEFORE 1e-4
        verbosity =SolverSettings.Verbosity
        niter=SolverSettings.LinearSolverIteration
        if(SolverSettings.NonSymmetricLinearSolver=="BICGSTAB_WITH_GMRES_FALLBACK"):
            self.non_symmetric_linear_solver = ExternalSolver.AMGCLSolver(ExternalSolver.AMGCLSmoother.ILU0,ExternalSolver.AMGCLIterativeSolverType.BICGSTAB_WITH_GMRES_FALLBACK, tol, niter, verbosity, gmres_size) #AMGCLIterativeSolverType.BICGSTAB_WITH_GMRES_FALLBACK, tol, 200, verbosity, gmres_size)#
        else:
            self.non_symmetric_linear_solver = ExternalSolver.AMGCLSolver(ExternalSolver.AMGCLSmoother.ILU0,ExternalSolver.AMGCLIterativeSolverType.GMRES, tol, niter, verbosity, gmres_size) #self.non_symmetric_linear_solver = AMGCLSolver(AMGCLSmoother.ILU0, AMGCLIterativeSolverType.GMRES, tol, 400, verbosity, gmres_size) #AMGCLIterativeSolverType.BICGSTAB_WITH_GMRES_FALLBACK, tol, 200, verbosity, gmres_size)#

        self.stage0_time_scheme = ResidualBasedIncrementalUpdateStaticScheme()
        self.stage0_conv_criteria = IncrementalDisplacementCriteria(SolverSettings.Stage0ConvCriteria1, SolverSettings.Stage0ConvCriteria2)
        self.stage0_max_iterations = SolverSettings.Stage0MaxIterations

        ##strategy to be used in step1 - diffusion + phase change - a symmetric linear solver is sufficient
        #tol = 1e-6
        #verbosity = 0
        self.symmetric_linear_solver = ExternalSolver.AMGCLSolver(ExternalSolver.AMGCLSmoother.ILU0, ExternalSolver.AMGCLIterativeSolverType.CG, tol, niter, verbosity, gmres_size)
        self.stage1_time_scheme = ResidualBasedIncrementalUpdateStaticVariablePropertyScheme()
        self.stage1_conv_criteria = IncrementalDisplacementCriteria(SolverSettings.Stage1ConvCriteria1, SolverSettings.Stage1ConvCriteria2) #IncrementalDisplacementCriteria(1e-4, 1e-6)
        self.stage1_max_iterations = SolverSettings.Stage1MaxIterations

        ##strategy to be used in step3 - diffusion + phase change - a symmetric linear solver is sufficient
        #self.stage2_time_scheme = ResidualBasedIncrementalUpdateStaticScheme()
        #self.stage2_conv_criteria = IncrementalDisplacementCriteria(1e-4,1e-6) #IncrementalDisplacementCriteria(1e-4, 1e-6)
        #self.stage2_max_iterations = 1
        self.MaxLineSearchIterations = SolverSettings.MaxLineSearchIterations
        Timer.Stop("ConvDiffPhaseChange_Constructor")
    def Initialize(self):
        Timer.Start("ConvDiffPhaseChange_Initialize")
        (self.model_part.ProcessInfo).SetValue(CONVECTION_DIFFUSION_SETTINGS, self.settings)

        self.stage0_solver = ResidualBasedNewtonRaphsonStrategy(self.model_part,
                                                                self.stage0_time_scheme,
                                                                self.non_symmetric_linear_solver,
                                                                self.stage0_conv_criteria,
                                                                self.stage0_max_iterations,
                                                                self.CalculateReactionFlag,
                                                                self.ReformDofSetAtEachStep,
                                                                self.MoveMeshFlag)
        if(self.use_linesearch_in_step1 == False):
            self.stage1_solver = ResidualBasedNewtonRaphsonStrategy(self.model_part,
                                                                    self.stage1_time_scheme,
                                                                    self.symmetric_linear_solver,
                                                                    self.stage1_conv_criteria,
                                                                    self.stage1_max_iterations,
                                                                    self.CalculateReactionFlag,
                                                                    self.ReformDofSetAtEachStep,
                                                                    self.MoveMeshFlag)
        else:

            self.tolls = 0.8           # energy tolerance factor on LineSearch (0.8 is ok)
            self.amp = 1.618         # maximum amplification factor
            self.etmxa = 3.0           # maximum allowed step length
            self.etmna = 0.01           # minimum allowed step length
            self.toler = 1.0E-9
            self.norm = 1.0E-6
            #self.stage1_solver = ResidualBasedNewtonRaphsonLineSearchesStrategy(self.model_part,
            #                                                                    self.stage1_time_scheme,
            #                                                                    self.symmetric_linear_solver,
            #                                                                    self.stage1_conv_criteria,
            #                                                                    self.stage1_max_iterations,
            #                                                                    self.MaxLineSearchIterations,
            #                                                                    self.tolls,
            #                                                                    self.amp,
            #                                                                    self.etmxa,
            #                                                                    self.etmna,
            #                                                                    self.CalculateReactionFlag,
            #                                                                    self.ReformDofSetAtEachStep,
            #                                                                    self.MoveMeshFlag,
            #                                                                    self.use_linesearch_in_step1)
            builder_and_solver=ResidualBasedBlockBuilderAndSolver(self.symmetric_linear_solver)
            self.stage1_solver = SolidMechanics.ResidualBasedNewtonRaphsonLineSearchStrategy(self.model_part,
                                                                                self.stage1_time_scheme,
                                                                                self.symmetric_linear_solver,
                                                                                self.stage1_conv_criteria,
                                                                                builder_and_solver,
                                                                                self.stage1_max_iterations,
                                                                                self.CalculateReactionFlag,
                                                                                self.ReformDofSetAtEachStep,
                                                                                self.MoveMeshFlag)

        self.Check()
        self.SetEchoLevel(self.echo_level)
        Timer.Stop("ConvDiffPhaseChange_Initialize")
    def Solve(self):
        Timer.Start("ConvDiffPhaseChange_Solve")
        print("")
        print("******************************")
        print("*  Solving for Thermal+P.C.  *")
        print("******************************")
        (self.model_part.ProcessInfo).SetValue(CONVECTION_DIFFUSION_SETTINGS, self.settings)

        #compute BDF2_coefficients
        ComputeBDFCoefficientsProcess(self.model_part,2).Execute()

        self.max_distance = 1e6
        #solve stage0 - pure convection step
        Timer.Start("ConvDiffPhaseChange_Solve_Stage0Solve")
        if(self.skip_stage0 == False):
            print("* Stage 0: Convection")
            self.activation_utils.ActivateElementsAndConditions( self.model_part, DISTANCE, self.max_distance, True)
            self.model_part.ProcessInfo.SetValue(FRACTIONAL_STEP,0)
            self.stage0_solver.Solve()
        else:
            BDFVector = self.model_part.ProcessInfo.GetValue(BDF_COEFFICIENTS)
            aux1 = BDFVector[1]/BDFVector[0]
            aux2 = BDFVector[2]/BDFVector[0]

            for node in self.model_part.Nodes:
                Told = node.GetSolutionStepValue(TEMPERATURE,1)
                Toldold = node.GetSolutionStepValue(TEMPERATURE,2)
                Tstar = -aux1*Told - aux2*Toldold
                node.SetSolutionStepValue(TEMPERATURE,0,Tstar)
        Timer.Stop("ConvDiffPhaseChange_Solve_Stage0Solve")
        Timer.Start("ConvDiffPhaseChange_Solve_Stage1Solve")
        if(self.skip_stage1 == False):
            #solve stage1 - diffusion + phase change - a symmetric linear solver is sufficient
            print("* Stage 1: Diffusion + P.C")
            #store "unknown" in the database without history
            self.variable_utils.SaveScalarVar(TEMPERATURE,TEMPERATURE,self.model_part.Nodes)
            self.activation_utils.ActivateElementsAndConditions( self.model_part, DISTANCE,  self.max_distance, True)  #0.0, True)
            self.model_part.ProcessInfo.SetValue(FRACTIONAL_STEP,1)
            self.stage1_solver.Solve()
        Timer.Stop("ConvDiffPhaseChange_Solve_Stage1Solve")
        Timer.Stop("ConvDiffPhaseChange_Solve")
    def SetEchoLevel(self, level):
        (self.stage0_solver).SetEchoLevel(level)
        (self.stage1_solver).SetEchoLevel(level)

    def Check(self):
        self.stage0_solver.Check()
        self.stage1_solver.Check()
        pass
