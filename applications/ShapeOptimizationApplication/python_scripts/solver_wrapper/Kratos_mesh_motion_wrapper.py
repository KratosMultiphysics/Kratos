##################################################################################################
# GENERAL MAIN OF THE ANALISYS
##################################################################################################

from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# time control starts
from time import *
print(ctime())

# Import the general variables read from the GiD
import mesh_motion_parameters as parameters

# setting the domain size for the problem to be solved
domain_size = parameters.domain_size

# including kratos path
from KratosMultiphysics import *

# including Applications paths
from KratosMultiphysics.MeshingApplication import *
from KratosMultiphysics.ALEApplication import *

# For GID output
from gid_output import GiDOutput

##################################################################################################

class Toolbox:

    # ********************************************************************************************
    def __init__(self):

        # defining the number of threads:
        num_threads = parameters.NumberofThreads
        parallel = OpenMPUtils()
        print("Num Threads = ", num_threads)
        parallel.SetNumThreads(int(num_threads))

        # problem settings
        problem_name = parameters.problem_name
        domain_size = parameters.domain_size
        input_file_name = parameters.problem_name
        self.time = parameters.start_time

        # defining a model part for the fluid
        self.mesh_model_part = ModelPart("MeshPart")
        self.mesh_model_part.AddNodalSolutionStepVariable(MESH_VELOCITY)

        # Importing mesh motion solver
        MeshSolverType = parameters.MeshSolverType   #Parameters coming from ProjectParameters.mdpa
        if(MeshSolverType == "Laplacian"):
            import mesh_solver_laplacian as mesh_solver_class
            self.mesh_solver_class = mesh_solver_class
            self.mesh_solver_class.AddVariables(self.mesh_model_part)
        elif(MeshSolverType == "StructuralSimilarity"):
            import mesh_solver_structural_similarity as mesh_solver_class
            self.mesh_solver_class = mesh_solver_class
            self.mesh_solver_class.AddVariables(self.mesh_model_part)
        elif(MeshSolverType == "StructuralSimilarityNonlinear"):
            import mesh_solver_structural_similarity_nonlinear as mesh_solver_class
            self.mesh_solver_class = mesh_solver_class
            self.mesh_solver_class.AddVariables(self.mesh_model_part)
        elif(MeshSolverType == "LaplacianComponentwise"):
            import mesh_solver_laplacian_componentwise as mesh_solver_class
            self.mesh_solver_class = mesh_solver_class
            self.mesh_solver_class.AddVariables(self.mesh_model_part)
        else:
            raise NameError("solver type not supported: options are LaplacianComponentwise  - StructuralSimilarity - StructuralSimilarityNonlinear")

        # reading the mesh model part
        model_part_io = ModelPartIO(input_file_name)
        model_part_io.ReadModelPart(self.mesh_model_part)
        self.mesh_model_part.SetBufferSize(parameters.buffer_size)
        self.mesh_model_part.ProcessInfo[DELTA_TIME] = 1

        #Creating the mesh solver
        reform_dofs_at_each_step = parameters.reform_dofs_at_each_step
        Tolerance = parameters.Tolerance
        MaxIter = parameters.MaxIter
        time_order = parameters.time_order

        self.mesh_solver_class.AddDofs(self.mesh_model_part)
        
        if(MeshSolverType == "Laplacian"):
          self.mesh_solver = mesh_solver_class.MeshSolverLaplacian(self.mesh_model_part,domain_size,reform_dofs_at_each_step)
          self.mesh_solver.time_order = time_order
          self.mesh_solver.Initialize()
        elif(MeshSolverType == "StructuralSimilarity"):
          self.mesh_solver = mesh_solver_class.MeshSolverStructuralSimilarity(self.mesh_model_part,reform_dofs_at_each_step)
          self.mesh_solver.time_order = time_order
          self.mesh_solver.Initialize()
        elif(MeshSolverType == "StructuralSimilarityNonlinear"):
          self.mesh_solver = mesh_solver_class.MeshSolverStructuralSimilarityNonlinear(self.mesh_model_part,domain_size,reform_dofs_at_each_step,Tolerance,MaxIter)
          self.mesh_solver.time_order = time_order
          self.mesh_solver.Initialize()
        elif(MeshSolverType == "LaplacianComponentwise"):
          self.mesh_solver = mesh_solver_class.MeshSolverLaplacianComponentwise(self.mesh_model_part,domain_size,reform_dofs_at_each_step)
          self.mesh_solver.time_order = time_order
          self.mesh_solver.Initialize()
        else:
          raise NameError("Selected solver has not been imported !")
        

    # ********************************************************************************************
    def initialize_gid_output(self):

        # For GID output
        self.gid_io = GiDOutput(parameters.problem_name,
                                parameters.VolumeOutput,
                                parameters.GiDPostMode,
                                parameters.GiDMultiFileFlag,
                                parameters.GiDWriteMeshFlag,
                                parameters.GiDWriteConditionsFlag)

        # Initialize design output in GID Format
        self.gid_io.initialize_results(self.mesh_model_part) 

    # ********************************************************************************************
    def finalize_gid_output(self):   
    
        # Finalize design output in GID formad
        self.gid_io.finalize_results()            

    # ********************************************************************************************
    def update_mesh(self,boundary_deformation):

        # Take time
        t0p = clock()

        # Iterate timer of solver
        self.time = self.time + 1
        self.mesh_model_part.CloneTimeStep(self.time)

        print ("Mesh-motion step = ", self.time)

        # Apply given boundary deformation as displacement (Dirichlet) condition
        self.apply_displacement_condition(boundary_deformation)

        # Solve mesh
        self.mesh_solver.Solve()

        # Write design in GID format
        self.gid_io.write_results(self.time, self.mesh_model_part, parameters.nodal_results, []) 

        print("Analysis Finalized ")

        # measure process time
        tfp = clock()

        print(ctime())
        print("Mesh-Motion Completed  [Analysis Time = ", tfp - t0p, "] ")

        # Return new mesh points
        new_mesh_points = {}
        for node in self.mesh_model_part.Nodes:
            new_mesh_points[node.Id] = [node.X,node.Y,node.Z]
        return new_mesh_points

    # ********************************************************************************************
    def update_reference_to_deformed_mesh(self):

        for node in self.mesh_model_part.Nodes:
            node.X0 = node.X
            node.Y0 = node.Y
            node.Z0 = node.Z

    # ********************************************************************************************
    def apply_displacement_condition(self,boundary_deformation):

        # Set boundary deformation update
        for node in self.mesh_model_part.GetNodes(3):
            disp = Vector(3)
            disp[0] = boundary_deformation[node.Id][0]
            disp[1] = boundary_deformation[node.Id][1]
            disp[2] = boundary_deformation[node.Id][2]
            node.SetSolutionStepValue(DISPLACEMENT,1,disp)
            node.Fix(DISPLACEMENT_X)
            node.Fix(DISPLACEMENT_Y)
            node.Fix(DISPLACEMENT_Z)

            if(node.Id == 1):
                print("node.Z = ",node.Z)

        # Treatment of fixed boundaries

        # farfield patch
        for node in self.mesh_model_part.GetNodes(1):
            node.Fix(DISPLACEMENT_X)
            node.Fix(DISPLACEMENT_Y)
            node.Fix(DISPLACEMENT_Z)
            disp = Vector(3)
            disp[0] = 0.0
            disp[1] = 0.0
            disp[2] = 0.0
            node.SetSolutionStepValue(DISPLACEMENT,1,disp)

        # Symmetry patch
        for node in self.mesh_model_part.GetNodes(2):
            node.SetSolutionStepValue(DISPLACEMENT_Y,1,0.0)
            node.Fix(DISPLACEMENT_Y)

    # ********************************************************************************************
