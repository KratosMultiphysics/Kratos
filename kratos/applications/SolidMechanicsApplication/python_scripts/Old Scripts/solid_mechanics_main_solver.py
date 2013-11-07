#importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
CheckForPreviousImport()


def AddVariables(model_part,rotation_dofs):
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part.AddNodalSolutionStepVariable(PRESSURE);
    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(ACCELERATION);
    model_part.AddNodalSolutionStepVariable(REACTION);
    model_part.AddNodalSolutionStepVariable(FORCE_INTERNAL);
    model_part.AddNodalSolutionStepVariable(FORCE_EXTERNAL);
    model_part.AddNodalSolutionStepVariable(FORCE_DYNAMIC);
    model_part.AddNodalSolutionStepVariable(REACTION_PRESSURE);

    AddExtraVariables(model_part,rotation_dofs);

    print "GENERAL VARIABLES ADDED CORRECTLY"

def AddExtraVariables(model_part,rotation_dofs):
    AddConditionVariables(model_part,rotation_dofs);
    #AddRigidWallVariables(model_part);
    #model_part.AddNodalSolutionStepVariable(MEAN_ERROR);
    #model_part.AddNodalSolutionStepVariable(OFFSET);
    #model_part.AddNodalSolutionStepVariable(NODAL_H);    
    #model_part.AddNodalSolutionStepVariable(NORMAL);
    #model_part.AddNodalSolutionStepVariable(FORCE_CONTACT_NORMAL);
    #model_part.AddNodalSolutionStepVariable(FORCE_CONTACT_TANGENT);
    
    print "EXTRA VARIABLES ADDED CORRECTLY"

def AddConditionVariables(model_part,rotation_dofs):
    #add specific variables for the problem (conditions)
    model_part.AddNodalSolutionStepVariable(IMPOSED_DISPLACEMENT);
    model_part.AddNodalSolutionStepVariable(POSITIVE_FACE_PRESSURE);
    model_part.AddNodalSolutionStepVariable(NEGATIVE_FACE_PRESSURE);
    model_part.AddNodalSolutionStepVariable(FORCE);
    model_part.AddNodalSolutionStepVariable(FACE_LOAD);
    model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);
    if(rotation_dofs == "True"):
        model_part.AddNodalSolutionStepVariable(ROTATION);
        model_part.AddNodalSolutionStepVariable(MOMENTUM);

    print "CONDITION VARIABLES ADDED CORRECTLY"

def AddRigidWallVariables(model_part):
    #add specific variables for the problem (rigid walls)
    model_part.AddNodalSolutionStepVariable(RIGID_WALL);
    model_part.AddNodalSolutionStepVariable(WALL_TIP_RADIUS);
    model_part.AddNodalSolutionStepVariable(WALL_REFERENCE_POINT);
    model_part.AddNodalSolutionStepVariable(WALL_VELOCITY);

    print "RIGID WALL VARIABLES ADDED CORRECTLY"
        
def AddDofs(model_part,problemtype,rotation_dofs):
    if(problemtype == "Mechanical"):
        if(rotation_dofs == "True"):
            for node in model_part.Nodes:
                node.AddDof(DISPLACEMENT_X,REACTION_X);
                node.AddDof(DISPLACEMENT_Y,REACTION_Y);
                node.AddDof(DISPLACEMENT_Z,REACTION_Z);
                node.AddDof(PRESSURE,REACTION_PRESSURE);
                node.AddDof(ROTATION_X,MOMENTUM_X);
                node.AddDof(ROTATION_Y,MOMENTUM_Y);
                node.AddDof(ROTATION_Z,MOMENTUM_Z);
        else:
            for node in model_part.Nodes:
                node.AddDof(DISPLACEMENT_X,REACTION_X);
                node.AddDof(DISPLACEMENT_Y,REACTION_Y);
                node.AddDof(DISPLACEMENT_Z,REACTION_Z);
                node.AddDof(PRESSURE,REACTION_PRESSURE);

    print "DOF'S ADDED CORRECTLY"


class SolidMechanicsSolver:
    #######################################################################
    def __init__(self,model_part,domain_size,echo_level):

        self.model_part = model_part
        
        #set the echo level
        self.echo_level = echo_level

       
        #definition of the default solver: -can be set using SetLinearSolver() method-
        self.linear_solver = SkylineLUFactorizationSolver()

        #definition of the default convergence criterion: -can be set using SetConvergenceCriterion() method-
        self.rel_tol   = 1e-4
        self.abs_tol   = 1e-9
        self.max_iters = 30;

        #self.mechanical_convergence_criterion  = DisplacementCriteria(self.rel_tol,self.abs_tol)
        self.mechanical_convergence_criterion  = ResidualCriteria(self.rel_tol,self.abs_tol)
        #self.mechanical_convergence_criterion = ResidualConvergenceCriteria(self.rel_tol,self.abs_tol)
        #self.mechanical_convergence_criterion = DisplacementConvergenceCriteria(self.rel_tol,self.abs_tol)

        #definition of the default builder_and_solver: 
        #for normal execution
        self.builder_and_solver = ResidualBasedBuilderAndSolver(self.linear_solver)
        #to conserve matrix blocks AMG solving
        #self.builder_and_solver = BlockResidualBasedBuilderAndSolver(self.linear_solver)
  
        #definition of computing flags
        self.CalculateReactionsFlag = True
        #Xref(n+1) = Xref(n)+dX
        self.MoveMeshFlag           = True


    #######################################################################
    def Initialize(self,load_restart):

        #creating the solution SCHEME:

        #type of solver (static=0, dynamic=1)
        if(self.solver_type == 0): 
            self.mechanical_scheme = ResidualBasedStaticScheme()
        elif(self.solver_type == 1):
            #definition of time scheme
            self.damp_factor_f  =  0.00; 
            self.damp_factor_m  = -0.01; 
            self.dynamic_factor =  1;
            self.mechanical_scheme = ResidualBasedBossakScheme(self.damp_factor_m,self.dynamic_factor)
        elif(self.solver_type == 2):
            #definition of time scheme
            self.damp_factor_f  =  0.00; 
            self.damp_factor_m  = -0.01; 
            self.dynamic_factor =  0;
            self.mechanical_scheme = ResidualBasedBossakScheme(self.damp_factor_m,self.dynamic_factor)
        elif(self.solver_type == 3):
            self.mechanical_scheme = ResidualBasedRelaxationScheme(-0.3,10.0)


        #creating the solution STRATEGY:
        self.ReformDofSetAtEachStep = False
                
        if(self.ComputeMechanicsFlag == True):
            self.ReformDofSetAtEachStep = True

        # option 1.- Application Strategy:
        self.mechanical_solver = ResidualBasedNewtonRaphsonStrategy(self.model_part,self.mechanical_scheme,self.linear_solver,self.mechanical_convergence_criterion,self.builder_and_solver,self.max_iters,self.CalculateReactionsFlag,self.ReformDofSetAtEachStep,self.MoveMeshFlag)
        #self.mechanical_solver = ResidualBasedNewtonRaphsonStrategy(self.model_part,self.mechanical_scheme,self.linear_solver,self.mechanical_convergence_criterion,self.max_iters,self.CalculateReactionsFlag,self.ReformDofSetAtEachStep,self.MoveMeshFlag

        # option 2.- Phyton Strategy:

        #import solid_mechanics_python_strategy
        
        #self.mechanical_solver = solid_mechanics_python_strategy.ResidualStrategy(self.model_part,self.mechanical_scheme,self.solid_linear_solver,self.mechanical_convergence_criterion,self.CalculateReactionsFlag,self.ReformDofSetAtEachStep,self.LineSearch)
           

        #check if everything is assigned correctly
        if(self.ComputeMechanicsFlag == True):
            self.mechanical_scheme.Check(self.model_part)
            self.mechanical_convergence_criterion.Check(self.model_part)
            self.mechanical_solver.Check();


        if(load_restart == "True"):
            self.mechanical_solver.SetInitialized()
            

                 
    #######################################################################   
    def Solve(self):

        if(self.ComputeMechanicsFlag == True):
            print " MECHANICAL SOLUTION START "
            (self.mechanical_solver).Solve()
            print " MECHANICAL SOLUTION PERFORMED "
             
                
        #move the mesh as needed
        #if(self.MoveMeshFlag == True):
        #    self.mechanical_scheme.MoveMesh(self.model_part.Nodes);

        ## MoveMesh is not a function of the scheme is a function of add strategies to python



    #######################################################################   
    def Check(self):
        self.builder_and_solver.Check(self.model_part)
        self.mechanical_scheme.Check(self.model_part)
        self.mechanical_convergence_criterion.Check(self.model_part)


    #######################################################################   
    def SetProblemType(self,problem_type,solver_type,line_search_type):
        
        self.LineSearchFlag = False
        self.ComputeMechanicsFlag = False;

        #type of problem (mechanical, thermal, thermo-mechanical)
        if(problem_type == "Mechanical"):
            self.ComputeMechanicsFlag = True;
        
        #type of solver (static=0, dynamic=1, psedo-static=2)
        if(solver_type == "StaticSolver"):
            print " Static Solver "
            self.solver_type = 0; 
        elif(solver_type == "DynamicSolver"):
            print " Dynamic Solver "
            self.solver_type = 1;
        elif(solver_type == "RelaxedDynamicSolver"):
            print " Pseudo-Static Solver "
            self.solver_type = 2;
        elif(solver_type == "RelaxationSolver"):
            print " Relaxation Solver "
            self.solver_type = 3;
            
            
        #line_search flag and type
        if(line_search_type  == "True"):
            self.LineSearchFlag = True
            

    #######################################################################   
    def SetLinearSolver(self,linear_solver_type,solver_tolerance,max_iters):

        class linear_solver_config:
            solver_type     = linear_solver_type
            scaling         = False
            tolerance       = solver_tolerance #1e-7
            max_iteration   = max_iters  #300     
            verbosity       = 0
            is_symmetric    = False  

            #Pastix Iterative Solver:
            gmres_krylov_space_dimension = 100
            ilu_level_of_fill            = 3 #5

            #GMRES or CG:
            preconditioner_type          = "None"

            #Deflated CG:
            assume_constant_structure    = True
            max_reduced_size             = 1000

            #AMG: (requires ResidualBasedBlockBuilderAndSolver )
            smoother_type  = "ILU0" #"DAMPED_JACOBI"
            krylov_type    = "GMRES"

        import linear_solver_factory

        self.linear_solver = linear_solver_factory.ConstructSolver(linear_solver_config)

        self.SetBuilderAndSolver(linear_solver_type)
                
    #######################################################################   
    def SetBuilderAndSolver(self,linear_solver_type):
        if(linear_solver_type == "AMGCL"):
            #to conserve matrix blocks AMG solving
            self.builder_and_solver = BlockResidualBasedBuilderAndSolver(self.linear_solver)      
        else:
            #for normal execution
            self.builder_and_solver = ResidualBasedBuilderAndSolver(self.linear_solver)
            

    #######################################################################   
    def SetConvergenceCriterion(self,convergence_criterion_type,convergence_tol,absolute_tol,max_iters, rotation_dofs):
        
        if(max_iters > 1):
            self.max_iters = max_iters
        #mechanical convergence criteria
        CT = convergence_tol;
        AT = absolute_tol; 

        if(rotation_dofs == "True"):
            if(convergence_criterion_type == "Displacement_criteria"):
                self.mechanical_convergence_criterion  =  DisplacementCriteria(CT,AT)
            elif(convergence_criterion_type == "Residual_criteria"):
                self.mechanical_convergence_criterion  =  ResidualCriteria(CT,AT)
            elif(convergence_criterion_type == "And_criteria"):
                Displacement   =   DisplacementCriteria(CT,AT)
                Residual       =   ResidualCriteria(CT,AT)
                self.mechanical_convergence_criterion  = AndCriteria(Residual, Displacement)
            elif(convergence_criterion_type == "Or_criteria"):
                Displacement   =   DisplacementCriteria(CT,AT)
                Residual       =   ResidualCriteria(CT,AT)
                self.mechanical_convergence_criterion  = OrCriteria(Residual, Displacement)
            elif(convergence_criterion_type == "Mixed_criteria"):
                Displacement   =   MixedElementCriteria(CT,AT)
                Residual       =   ResidualCriteria(CT,AT)
                self.mechanical_convergence_criterion  = AndCriteria(Residual, Displacement)
        else:
            if(convergence_criterion_type == "Displacement_criteria"):
                self.mechanical_convergence_criterion  =  DisplacementConvergenceCriteria(CT,AT)
            elif(convergence_criterion_type == "Residual_criteria"):
                self.mechanical_convergence_criterion  =  ResidualConvergenceCriteria(CT,AT)
            elif(convergence_criterion_type == "And_criteria"):
                Displacement   =   DisplacementConvergenceCriteria(CT,AT)
                Residual       =   ResidualConvergenceCriteria(CT,AT)
                self.mechanical_convergence_criterion  = AndCriteria(Residual, Displacement)
            elif(convergence_criterion_type == "Or_criteria"):
                Displacement   =   DisplacementConvergenceCriteria(CT,AT)
                Residual       =   ResidualConvergenceCriteria(CT,AT)
                self.mechanical_convergence_criterion  = OrCriteria(Residual, Displacement)
            elif(convergence_criterion_type == "Mixed_criteria"):
                Displacement   =   MixedElementConvergenceCriteria(CT,AT)
                Residual       =   ResidualConvergenceCriteria(CT,AT)
                self.mechanical_convergence_criterion  = AndCriteria(Residual, Displacement)

        #move the mesh as needed
        #if(self.MoveMeshFlag == True):
        #    self.mechanical_scheme.MoveMesh(self.model_part.Nodes);

        ## MoveMesh is not a function of the scheme is a function of add strategies to python
