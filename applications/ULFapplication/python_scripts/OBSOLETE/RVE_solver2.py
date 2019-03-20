from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.ULFApplication as KratosULFapp

## Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

## Import base class file
import navier_stokes_solver_vmsmonolithic

def CreateSolver(main_model_part, custom_settings):
    return RVEMonolithicSolver2(main_model_part, custom_settings)


BaseAlgorithm = navier_stokes_solver_vmsmonolithic.NavierStokesSolverMonolithic


class RVEMonolithicSolver2(BaseAlgorithm):

    def __init__(self, main_model_part, custom_settings):
        super(RVEMonolithicSolver2,self).__init__(main_model_part, custom_settings)
        self.rel_vel_tol = 1e-3
        self.abs_vel_tol = 1e-6
        self.rel_pres_tol = 1e-3
        self.abs_pres_tol = 1e-6
        self.dynamic_tau = 0.0
        self.oss_switch = 0
        self.alpha = -0.3
        #2 Lagr 0 Eul #1 ALE
        #self.move_mesh_strategy = 2
        self.move_mesh_strategy = 0
        
        self.regularization_coef = 1000
        self.max_iter = 30
        

        # default settings
        self.echo_level = 0
        self.compute_reactions = False
        self.ReformDofSetAtEachStep = False
        self.CalculateNormDxFlag = False
        self.MoveMeshFlag = False
        self.use_slip_conditions = False
        
        #pDiagPrecond = DiagonalPreconditioner()
        #self.linear_solver = BICGSTABSolver(1e-6, 5000, pDiagPrecond)
        self.linear_solver=SuperLUSolver()
        
        


    def Initialize(self):

        self.computing_model_part = self.GetComputingModelPart()      
        # Creating the solution strategy
        self.conv_criteria = VelPrCriteria(self.rel_vel_tol, self.abs_vel_tol, self.rel_pres_tol, self.abs_pres_tol)
        self.time_scheme = FluidDynamicsApplication.ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent( self.alpha,                                        
                                        self.domain_size,
                                        FluidDynamicsApplication.PATCH_INDEX
                                        )


        builder_and_solver = KratosULFapp.ResidualBasedEliminationBuilderAndSolverPeriodicNormalOnly(self.linear_solver,
                                                                                KratosCFD.PATCH_INDEX)

        self.solver = ResidualBasedNewtonRaphsonStrategy(
            self.model_part, self.time_scheme, self.linear_solver, self.conv_criteria,
            builder_and_solver, self.max_iter, self.compute_reactions, self.ReformDofSetAtEachStep, self.MoveMeshFlag)

        self.solver = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(self.main_model_part,
                                                                            self.time_scheme,
                                                                            self.linear_solver,
                                                                            self.conv_criteria,
                                                                            builder_and_solver,
                                                                            self.max_iter,
                                                                            self.compute_reactions,
                                                                            self.ReformDofSetAtEachStep,
                                                                            self.MoveMeshFlag)

        (self.solver).SetEchoLevel(self.echo_level)
        self.solver.Check()   

        ##########################################################
        #self.Remesh()        
        #(self.append_model_part_process).AppendPart(self.model_part, solid_model_part);        
        (self.neigh_finder).Execute();       
        #we need normals to prescribe the inlet velocity
        #self.normal_util = NormalCalculationUtils()
        
        #self.normal_util.CalculateOnSimplex(self.model_part.Conditions, self.domain_size)
        normal_tools = BodyNormalCalculationUtils()
        normal_tools.CalculateBodyNormals(self.model_part, self.domain_size)  
        
        #self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DYNAMIC_TAU, self.settings["dynamic_tau"].GetDouble())
        #self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.OSS_SWITCH, self.settings["oss_switch"].GetInt())
        #self.process_data[KratosMultiphysics.DYNAMIC_TAU] = self.settings["dynamic_tau"].GetDouble()
        (self.solver).Initialize()

        print ("Monolithic solver initialization finished.")

