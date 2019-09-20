domain_size = 3

# Fluid solver configuration
class StructuralSolverConfiguration:
    domain_size = 3
    solver_type =  "structural_solver_static"

    #direct solver
    class linear_solver_config:
        solver_type = "Super LU"

    #decomment this to use an iterative solver
    #class linear_solver_config:
        #solver_type = "AMGCL"
        #tolerance = 1E-4
        #max_iteration = 400
        ##preconditioner = "Diagonal"
        #scaling = False
        #smoother_type = "ILU0"
        #krylov_type = "BICGSTAB2"
        #verbosity = 1
        
#general problem settings
Dt = 0.1
Start_time = 0.0
max_time = 20.0
nsteps = 1000

                   
#output settings
output_time = 0.001
output_step = 1
VolumeOutput = True
nodal_results=["VELOCITY","TEMPERATURE"]
gauss_points_results=[]
GiDPostMode = "Binary"
GiDWriteMeshFlag = True
GiDWriteConditionsFlag = True
GiDWriteParticlesFlag = True
GiDMultiFileFlag = "Single"

problem_name="hexa_laplacian"
problem_path="/home/riccardo/examples/cyl_3d.gid"
kratos_path="/home/riccardo/Kratos"
