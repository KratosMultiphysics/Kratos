domain_size = 2

#SolverType = "FractionalStep"
NumberofThreads = 8
MeshSolverType = "StructuralSimilarity" #StructuralSimilarity  #Laplacian


#general problem settings
AutomaticDeltaTime = "Fixed"
Dt = 0.01
Start_time = 0.0
max_time = 1
nsteps = 1000


groups_dictionary = {
        "Mesh" : 1,
                   }
#output settings
output_time = 0.01
output_step = 1
VolumeOutput = True
nodal_results=["DISPLACEMENT", "MESH_VELOCITY"]
gauss_points_results=[]
GiDPostMode = "Binary"
GiDWriteMeshFlag = True
GiDWriteConditionsFlag = True
GiDWriteParticlesFlag = False
GiDMultiFileFlag = "Single"

problem_name="03_MeshTest_3D8N"
kratos_path="D:\Kratos"
