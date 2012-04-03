Dt = 0.0005

##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size =3


import time 
from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.MeshingApplication import *

#defining a model part
model_part = ModelPart("FluidPart");  

#providing the variable list to the model part
import fractional_step_solver
fractional_step_solver.AddVariables(model_part)
model_part.AddNodalSolutionStepVariable(REACTION)

#introducing input file name
input_file_name = "cavity3D"

#reading the fluid part
gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
gid_io = GidIO(input_file_name,gid_mode,multifile,deformed_mesh_flag, write_conditions)
model_part_io_fluid = ModelPartIO(input_file_name)
model_part_io_fluid.ReadModelPart(model_part)

print model_part

#the buffer size should be set up here after the mesh is read for the first time
model_part.SetBufferSize(3)

#adding dofs
fractional_step_solver.AddDofs(model_part)



#creating a fluid solver object
fluid_solver = fractional_step_solver.IncompressibleFluidSolver(model_part,domain_size)
#fluid_solver = flexible_incompressible_fluid_solver.IncompressibleFluidSolver(model_part,domain_size)
fluid_solver.laplacian_form =1;
fluid_solver.predictor_corrector = True;
fluid_solver.vel_toll = 1e-3
fluid_solver.time_order = 2
fluid_solver.ReformDofAtEachIteration = False
fluid_solver.echo_level = 0
fluid_solver.compute_reactions = True


##pILUPrecond = ILU0Preconditioner() 
##pDiagPrecond = DiagonalPreconditioner() 
##fluid_solver.pressure_linear_solver =  BICGSTABSolver(1e-3, 5000,pDiagPrecond)
##fluid_solver.pressure_linear_solver = SkylineLUFactorizationSolver();

##for node in model_part.Nodes:
##    if(node.X < 0.001 and node.Y<0.001):
##        node.Fix(PRESSURE)

fluid_solver.Initialize()

def fu(x):
    return (x**2) * ( (1.0-x)**2 )

def df(x):
    return 2.0*x*( (1.0-x)**2 ) - 2.0*(x**2)*(1.0-x)

def d2f(x):
    return 2.0*( (1.0-x)**2 ) - 8.0*x*(1.0-x) + 2.0*(x**2)

def d3f(x):
    return -12.0*(1.0-x) + 12.0*x
    


def MoveNodes(nodes,Dt,substeps,h):
  dt_small = Dt / float(substeps)
  cfl = 0.0
  
  for node in nodes:
    if(node.IsFixed(VELOCITY_X) == False):
      for step in range(0,substeps):
        x = node.X
        y = node.Y
        ht = 1.0 #h(t)
        vx = 100.0 * ht * fu(x) * df(y)
        vy = -100.0 * ht * fu(y) * df(x)
        node.X = x + vx*dt_small
        node.Y = y + vy*dt_small
		
	    
      
  
outfile = open("reactions.out",'w')

def WriteReactions(time,nodes):
    outfile.write("time = " + str(time) + "\n" )
    for node in model_part.Nodes:
        vel = node.GetSolutionStepValue(REACTION)
        a = str(node.Id) + " " + str(vel[0]) + " " +  str(vel[1]) + " " + str(vel[2]) + "\n"
        outfile.write(a)
    outfile.write("\n")
    
#settings to be changed
nsteps = 100



##zero = Vector(3);
##zero[0] = 0.0;
##zero[1] = 0.0;
##zero[2] = 0.0;
for step in range(0,50):
	print "line49"

	timeX = Dt*step
	model_part.CloneTimeStep(timeX)
	
	reconnector = TetrahedraReconnectUtility(model_part)
	
	print timeX
    #print model_part.ProcessInfo()[TIME]

    #solving the fluid problem
	substeps = 10
	MoveNodes(model_part.Nodes,Dt,substeps,1)

    #do mesh improvement
	start = time.time()
	print "trying to obtain quality"
	
    #Parameters ModelPart& r_model_part, 
    #						int iterations , 
    #                        bool processByNode, 
    #						bool ProcessByFace, 
    #						bool ProcessByEdge, 
    #						 bool saveToFile, 
    #						bool removeFreeVertexes
    #						bool doInParallel
	#						bool reInsertNodes
	#						bool debugMode
	ProcessByFace=1
	ProcessByEdge=1
	ProcessByNode=0
	reInsertNodes = 0
	numOptIterations = 3
	debugMode = 0
	saveToFile = 0
	doInParallel = 0
	reInsertNodes = 0
	removeFreeVertexes = 0
	start2 = time.time()
	reconnector.setMaxNumThreads(16)
	reconnector.setBlockSize(16)
	reconnector.OptimizeQuality(model_part , numOptIterations ,ProcessByNode,ProcessByFace,ProcessByEdge, saveToFile, removeFreeVertexes,doInParallel,reInsertNodes,debugMode)
	elapsed2 = (time.time() - start2)
	meshIsValid =reconnector.EvaluateQuality()	 
	start3 = time.time()
        reconnector.FinalizeOptimization(removeFreeVertexes)
	elapsed3 = (time.time() - start3)
	
	
	
	elapsed = (time.time() - start)
	print "*************************************************"
	print "Elapsed time.Total:"+ '{0:2f} '.format(elapsed)+ " Reconnect"+'{0:2f} '.format(elapsed2)+ " Del:"+'{0:2f} '.format(elapsed3)
	
	print "*************************************************"
	mesh_name = timeX #if we want the mesh to change at each time step then ****mesh_name = time****
	gid_io.InitializeMesh( mesh_name)
	gid_io.WriteMesh( model_part.GetMesh() )
	gid_io.FinalizeMesh()
	gid_io.InitializeResults(mesh_name , model_part.GetMesh())


	#print the results
	gid_io.WriteNodalResults(PRESSURE,model_part.Nodes,timeX,0)
	gid_io.WriteNodalResults(VELOCITY,model_part.Nodes,timeX,0)
	gid_io.WriteNodalResults(REACTION,model_part.Nodes,timeX,0)

	gid_io.FinalizeResults()
	
	if not meshIsValid:
		print 'Mesh has negative elements. Reached iteration' , step
		break
	


          
        

