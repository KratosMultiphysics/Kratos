#Problem Data
domain_size = 2
NumberofThreads = 6
buffer_size = 3
start_time = 0
MeshSolverType = "StructuralSimilarity" #StructuralSimilarity   #StructuralSimilarityNonlinear  #Laplacian
problem_name="mesh_ONERAM6_inv_ffd"

# Mesh solver
reform_dofs_at_each_step = False
Tolerance = 10e-5
MaxIter = 40
time_order = 2

# Post-processing
nodal_results=["MESH_VELOCITY","DISPLACEMENT"]
VolumeOutput = True
GiDPostMode = "Binary"
GiDWriteMeshFlag = True
GiDWriteConditionsFlag = True
GiDWriteParticlesFlag = False
GiDMultiFileFlag = "Single"
