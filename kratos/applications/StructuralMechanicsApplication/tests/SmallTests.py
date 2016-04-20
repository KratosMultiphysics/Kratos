# Import Kratos
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.StructuralMechanicsApplication import *

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

# importing the process factory
import process_factory

def GetFilePath(fileName):
  return os.path.dirname(__file__) + "/" + fileName

class ConstitutiveLawBuild:
  def __init__(self, model_part, domain_size):
    self.model_part = model_part
    self.domain_size = domain_size

  def Initialize(self):
    self.SetConstitutiveLaw()

  def SetConstitutiveLaw(self):
    AssignMaterial(self.model_part.Properties)

def AssignMaterial(Properties):
    prop_id = 1;
    prop = Properties[prop_id]
    mat = LinearElastic3DLaw();
    prop.SetValue(CONSTITUTIVE_LAW, mat.Clone());

def GeneralModelPartBuild(directory_name,file_name):
  # Modelpart for the solid
  model_part = ModelPart("StructuralPart")
  model_part.AddNodalSolutionStepVariable(ALPHA_EAS)
  model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
  model_part.AddNodalSolutionStepVariable(REACTION)

  # Initialize GiD  I/O
  input_file_name = GetFilePath(directory_name + "/" + file_name)

  # Reading the fluid part
  model_part_io_fluid = ModelPartIO(input_file_name)
  model_part_io_fluid.ReadModelPart(model_part)

  # Setting up the buffer size
  model_part.SetBufferSize(3)

  # Adding dofs
  for node in model_part.Nodes:
    node.AddDof(DISPLACEMENT_X, REACTION_X)
    node.AddDof(DISPLACEMENT_Y, REACTION_Y)
    node.AddDof(DISPLACEMENT_Z, REACTION_Z)

  return model_part
      
def GeneralSolverBuild(model_part, domain_size):
  # Set the constitutive law
  constitutive_law = ConstitutiveLawBuild(model_part, domain_size);
  constitutive_law.Initialize();

  max_iters = 30

  linear_solver = SkylineLUFactorizationSolver()

  # Create solver
  builder_and_solver = ResidualBasedBuilderAndSolver(linear_solver)
  mechanical_scheme = ResidualBasedStaticScheme()
  mechanical_convergence_criterion = ResidualCriteria(1e-4, 1e-4)

  mechanical_solver = ResidualBasedNewtonRaphsonStrategy(
            model_part,
            mechanical_scheme,
            linear_solver,
            mechanical_convergence_criterion,
            builder_and_solver,
            max_iters,
            False ,
            True,
            True)

  mechanical_solver.SetEchoLevel(0)
  mechanical_solver.Initialize()
  
  return mechanical_solver

def GeneralSolver(model_part, mechanical_solver, myprocesslist, Dt, final_time):
  
  step = 0
  time = 0.0
  
  while(time <= final_time):

    time += Dt
    step += 1
    
    model_part.CloneTimeStep(time)
    
    for process in myprocesslist:    
        process.ExecuteInitializeSolutionStep()

    mechanical_solver.Solve()
  
  for process in myprocesslist:    
      process.ExecuteFinalize()
  
  return model_part
  
class SprismPatchTests(KratosUnittest.TestCase):

    def setUp(self):
      # Creating the model part
      self.model_part = GeneralModelPartBuild("SPRISM3D6N","patch_test")
      
      # Create the list of processes 
      my_process_list_aux = [{"process_name":"SPRISM_process",
                          "implemented_in_python":True, 
                          "implemented_in_module":"KratosMultiphysics", 
                          "implemented_in_file":"sprism_process", 
                          "parameters":{"model_part_name":"model_part"}}]  
      my_process_list = Parameters(my_process_list_aux)
      Model = { "model_part" : self.model_part}
      process_constructor = process_factory.KratosProcessFactory(Model)
      self.myprocesslist = process_constructor.ConstructListOfProcesses(my_process_list)
      
      # Execute the initial processes
      for process in self.myprocesslist:
          process.ExecuteInitialize()
      
      # Adding the solver
      self.mechanical_solver = GeneralSolverBuild(self.model_part, 3)

    def test_MembranePatch(self):
    
        # Add BC
        for node in self.model_part.Nodes:
            if (node.X >2.40000e-01 -1.0e-5) | (node.Y > 1.20000e-01 -1.0e-5) | (node.X < 1.0e-5) | (node.Y < 1.0e-5):
                node.Fix(DISPLACEMENT_X)
                node.SetSolutionStepValue(DISPLACEMENT_X, 0, 1.0e-7 * (node.X + node.Y / 2))
                node.Fix(DISPLACEMENT_Y)
                node.SetSolutionStepValue(DISPLACEMENT_Y, 0, 1.0e-7 * (node.Y + node.X / 2))

        Dt = 1.0
        final_time = 1.0
        
        self.model_part = GeneralSolver(self.model_part, self.mechanical_solver, self.myprocesslist, Dt, final_time)

        for node in self.model_part.Nodes:
            value = node.GetSolutionStepValue(DISPLACEMENT_X,0)
            self.assertAlmostEqual(value, 1.0e-7 * (node.X + node.Y / 2))
            value = node.GetSolutionStepValue(DISPLACEMENT_Y,0)
            self.assertAlmostEqual(value, 1.0e-7 * (node.Y + node.X / 2))
            
    def test_BendingPatch(self):
        # Add BC
        for node in self.model_part.Nodes:
            if (node.X >2.40000e-01 -1.0e-5) | (node.Y > 1.20000e-01 -1.0e-5) | (node.X < 1.0e-5) | (node.Y < 1.0e-5):
                node.Fix(DISPLACEMENT_X)
                node.SetSolutionStepValue(DISPLACEMENT_X, 0, -1.0e-7 * (node.Z - 0.0005) * (node.X + node.Y / 2))
                node.Fix(DISPLACEMENT_Y)
                node.SetSolutionStepValue(DISPLACEMENT_Y, 0, -1.0e-7 * (node.Z - 0.0005) * (node.Y + node.X / 2))
                node.Fix(DISPLACEMENT_Z)
                node.SetSolutionStepValue(DISPLACEMENT_Z, 0, 0.5 * 1.0e-7 * (node.X ** 2 + node.X * node.Y + node.Y ** 2))

        Dt = 1.0
        final_time = 1.0
        
        self.model_part = GeneralSolver(self.model_part, self.mechanical_solver, self.myprocesslist, Dt, final_time)

        for node in self.model_part.Nodes:
            value = node.GetSolutionStepValue(DISPLACEMENT_X,0)
            self.assertAlmostEqual(value, -1.0e-7 * (node.Z - 0.0005) * (node.X + node.Y / 2))
            value = node.GetSolutionStepValue(DISPLACEMENT_Y,0)
            self.assertAlmostEqual(value,-1.0e-7 * (node.Z - 0.0005) * (node.Y + node.X / 2))
            value = node.GetSolutionStepValue(DISPLACEMENT_Z,0)
            self.assertAlmostEqual(value, 0.5 * 1.0e-7 * (node.X ** 2 + node.X * node.Y + node.Y **2 ))

    def tearDown(self):
        pass