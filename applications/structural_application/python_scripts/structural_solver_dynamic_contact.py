#importing the Kratos Library
from Kratos import *
from KratosStructuralApplication import *
from KratosIncompressibleFluidApplication import *
from KratosMeshingApplication import *
#from KratosExternalSolversApplication import *


def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_OLD);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_NULL);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_EINS);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_DT);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_NULL_DT);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_EINS_DT);
    model_part.AddNodalSolutionStepVariable(ACCELERATION_NULL);
    model_part.AddNodalSolutionStepVariable(ACCELERATION_EINS);
    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(ACCELERATION);
    model_part.AddNodalSolutionStepVariable(REACTION);
    model_part.AddNodalSolutionStepVariable(NEGATIVE_FACE_PRESSURE);
    model_part.AddNodalSolutionStepVariable(POSITIVE_FACE_PRESSURE);
    model_part.AddNodalSolutionStepVariable(INSITU_STRESS);
    model_part.AddNodalSolutionStepVariable(FACE_LOAD);
    model_part.AddNodalSolutionStepVariable(IS_STRUCTURE);
    model_part.AddNodalSolutionStepVariable(IS_FLUID);

    print "variables for the dynamic structural solution added correctly"
        
def AddDofs(model_part):
    for node in model_part.Nodes:
        #adding dofs
##        node.AddDof(DISPLACEMENT_X);
##        node.AddDof(DISPLACEMENT_Y);
##        node.AddDof(DISPLACEMENT_Z);
        node.AddDof(DISPLACEMENT_X,REACTION_X);
        node.AddDof(DISPLACEMENT_Y,REACTION_Y);
        node.AddDof(DISPLACEMENT_Z,REACTION_Z);
    print "dofs for the dynamic structural solution added correctly"


class DynamicStructuralSolver:
    #######################################################################
    def __init__(self,model_part,domain_size):

        self.model_part = model_part
        self.echo_level = 0
        
        self.damp_factor = -0.1
        self.toll = 1e-6
        self.absolute_tol = 1e-9

        self.Ide                        = 20
	self.factor_delta_lmax          = 1.00
	self.toler                      = 1.0E-9
        self.norm                       = 1.0E-5 
	self.MaxIterations              = 100
        self.damp_factor                = 0; 


        self.CalculateReactionFlag  = True
        self.ReformDofSetAtEachStep = True
        self.MoveMeshFlag           = True
	self.ApplyBodyForce         = True
        

        #definition of the solvers
        self.structure_linear_solver =  SkylineLUFactorizationSolver()
        ##self.structure_linear_solver      =   SuperLUSolver()

        #definition of the convergence criteria
        self.conv_criteria = DisplacementCriteria(0.000001,1e-9)
        
        #added for contact
	self.shell_model_part = ModelPart("shell_model_part");
	self.contact_model_part = ModelPart("contact_model_part");
	self.contact_model_part.Properties = self.model_part.Properties

        self.save_structure_model_part_process = SaveStructureModelPartProcess();
        self.save_structure_conditions_process = SaveStructureConditionsProcess();
        self.merge_in_one_model_parts_process = MergeInOneModelPartsProcess();  
        self.domain_size = domain_size
        self.Mesher =TetGenPfemContact()

    #######################################################################
    def Initialize(self):

      
        #definition of time scheme
        self.time_scheme = ResidualBasedPredictorCorrectorBossakScheme(self.damp_factor)
        #definition of the convergence criteria
        #self.conv_criteria = DisplacementCriteria(self.toll,self.absolute_tol)
        #builder_and_solver = ResidualBasedEliminationBuilderAndSolver(self.structure_linear_solver)

        #creating the solution strategy
        CalculateReactionFlag = False
        ReformDofSetAtEachStep = True
        MoveMeshFlag = True
        import strategy_python
        self.solver = strategy_python.SolvingStrategyPython(self.model_part,self.time_scheme,self.structure_linear_solver,self.conv_criteria,CalculateReactionFlag,ReformDofSetAtEachStep,MoveMeshFlag)
       

	#flags necessary for save_structure
	for node in (self.model_part).Nodes:
	    node.SetSolutionStepValue(IS_STRUCTURE, 0 , 1.0) 
	    node.SetSolutionStepValue(IS_FLUID, 0 , 0.0) 

        (self.save_structure_model_part_process).SaveStructureModelPart(self.model_part, self.shell_model_part, self.domain_size); 
	#shell_elements = ElementsArray()
	#(SaveElementBySizeProcess((self.model_part).Elements, shell_elements, 3)).Execute()
                  
	#(self.shell_model_part).Elements = shell_elements      
           #creating the solution strategy
         #self.solver = ResidualBasedNewtonRaphsonStrategy(self.model_part,self.time_scheme,self.structure_linear_solver,self.conv_criteria,30,True,False,True)
#        print "self.echo_level = " , self.echo_level
##        #(self.solver).SetBuilderAndSolver(builder_and_solver)
##        (self.solver).SetEchoLevel(self.echo_level)
##        #(self.solver).SetReformDofSetAtEachStepFlag(True)
##        #(self.solver).SetMoveMeshFlag(True)
##        print "finished initialization of the fluid strategy"

                 
    #######################################################################   
    def Solve(self):
	self.ContactMesh()
	print "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<after contact mesh"
        (self.solver).Solve()

    #######################################################################   
    def SetEchoLevel(self,level):
        (self.solver).SetEchoLevel(level)

    #######################################################################   
    def ContactMesh(self):
	#print (self.model_part).Elements
	#(self.shell_model_part).Elements = (self.model_part).Elements
	shell_elements = ElementsArray()
	(SaveElementBySizeProcess((self.model_part).Elements, shell_elements, 3)).Execute()
                  
	(self.contact_model_part).Elements = shell_elements   

	print ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<inside contact mesh"
	#print self.contact_model_part
	print "After print"
	(self.Mesher).ReGenerateMesh("PfemContactElement3D","Face3D3N",self.contact_model_part)	

	(self.model_part).Elements = ElementsArray()
	#print self.model_part
	#print self.shell_model_part
	#print self.contact_model_part

	print "before merge"
        (self.merge_in_one_model_parts_process).MergeParts(self.model_part, self.shell_model_part)
        (self.merge_in_one_model_parts_process).MergeParts(self.model_part, self.contact_model_part)
	print "after Merge"
	#print self.model_part





