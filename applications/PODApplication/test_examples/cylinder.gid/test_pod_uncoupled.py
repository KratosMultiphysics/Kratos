# -*- coding: utf-8 -*-
import string
import scipy
from scipy import array

import fluid_only_var

##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = fluid_only_var.domain_size

##################################################################
##################################################################
## ATTENTION: here the order is important

#including kratos path
kratos_libs_path            = '../../../../libs' ##kratos_root/libs
kratos_applications_path    = '../../../../applications' ##kratos_root/applications
import sys
sys.path.append(kratos_libs_path)
sys.path.append(kratos_applications_path)

#importing Kratos main library
from Kratos import *
kernel = Kernel()   #defining kernel

#importing applications
import applications_interface
applications_interface.Import_IncompressibleFluidApplication = True
applications_interface.Import_ExternalSolversApplication = True
applications_interface.Import_PodApplication = True
applications_interface.ImportApplications(kernel, kratos_applications_path)

## from now on the order is not anymore crucial
##################################################################
##################################################################
from KratosIncompressibleFluidApplication import *
from KratosExternalSolversApplication import *
from KratosPodApplication import *

#print kernel

#defining a model part for the fluid and one for the structure
fluid_model_part = ModelPart("FluidPart");  

#############################################
##importing the solvers needed
SolverType = fluid_only_var.SolverType
if(SolverType == "FractionalStep"):
    import incompressible_fluid_solver
    incompressible_fluid_solver.AddVariables(fluid_model_part)
elif(SolverType == "monolithic_solver_eulerian"):
    import monolithic_solver_eulerian_standard_biggerlocalM
    monolithic_solver_eulerian_standard_biggerlocalM.AddVariables(fluid_model_part)
elif(SolverType == "monolithic_solver_eulerian_compressible"):
    import monolithic_solver_eulerian_compressible
    monolithic_solver_eulerian_compressible.AddVariables(fluid_model_part)
else:
    raise "solver type not supported: options are FractionalStep - Monolithic"

#introducing input file name
input_file_name = fluid_only_var.problem_name

#reading the fluid part
gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
gid_io = GidIO(input_file_name,gid_mode,multifile,deformed_mesh_flag, write_conditions)
model_part_io_fluid = ModelPartIO(input_file_name)
model_part_io_fluid.ReadModelPart(fluid_model_part)

#setting up the buffer size: SHOULD BE DONE AFTER READING!!!
fluid_model_part.SetBufferSize(3)

##adding dofs
if(SolverType == "FractionalStep"):
    incompressible_fluid_solver.AddDofs(fluid_model_part)
elif(SolverType == "monolithic_solver_eulerian"):
    monolithic_solver_eulerian_standard_biggerlocalM.AddDofs(fluid_model_part)
elif(SolverType == "monolithic_solver_eulerian_compressible"):
    monolithic_solver_eulerian_compressible.AddDofs(fluid_model_part)

#########select here the laplacian form!!!!!!!!!!!!!!!!!
laplacian_form = fluid_only_var.laplacian_form 
if(laplacian_form >= 2):
    for node in fluid_model_part.Nodes:
        node.Free(PRESSURE)

##check to ensure that no node has zero density or pressure
for node in fluid_model_part.Nodes:
    if(node.GetSolutionStepValue(DENSITY) == 0.0):
        print "node ",node.Id," has zero density!"
        raise 'node with zero density found'
    if(node.GetSolutionStepValue(VISCOSITY) == 0.0):
        print "node ",node.Id," has zero viscosity!"
        raise 'node with zero VISCOSITY found'

#reading the intial values

def Intial_Values_Read(infile,Number_Nodes):
      
      File = open(infile,'r')
      #IN=[]
      #in_=Vector(Number_Pods)
      #in_boundary=Vector(Number_Pods)
      t_=1
      print "begin reading Intial VX values"  
      #i=0
      while (t_<=Number_Nodes):
	   searchname_=str('Result "Velocity_avrged" "Kratos" 1 Vector OnNodes')
           line = File.readline()
           #print "searching for:   " +searchname_        
           if not line:
             #print "aaaaaaaaaaaaaaaaaaaaaA"
             break
           else:   
	     if(line.find(searchname_) >= 0):
	        print "begin reading "+searchname_ 
                line = File.readline()
                line = File.readline()
                j=0
                if(line.find(str('Values')) >= 0):
		   go_on = True 
                   while go_on:
                     line = File.readline()
                     if not line:
                        break
                   
                     if (line.find("End values") >= 0):
                        print 'END Reading elements'
                        go_on = False;
                        
                        #searchname_=str('Result "')+Quantity_+str ('" "Kratos" ')+str(t_)+str(' ')+Dimension_+str(' OnNodes')
                        #print searchname_
                     else:
		        aaa = line.split();
		        #print aaa
		        #in_[j]=float(aaa[0])
		        j=j+1
		        fluid_model_part.Nodes[j].SetValue(VELOCITY_X,float(aaa[1]))
		        fluid_model_part.Nodes[j].SetValue(VELOCITY_Y,float(aaa[2]))
#		        fluid_model_part.Nodes[j].SetSolutionStepValue(VELOCITY_X,0,float(aaa[1]))
		        #fluid_model_part.Nodes[j].SetSolutionStepValue(VELOCITY_Y,0,float(aaa[2]))
		       
		   t_=t_+1
      File.close()	   
      #File = open(infile,'r')
      #t_=1
      #print "begin reading Intial VY values"  
      ##i=0
      #while (t_<=Number_Nodes):
	   #searchname_=str('Result "VelocityY_average" "Kratos" 1 Scalar OnNodes')
           #line = File.readline()
           ##print "searching for:   " +searchname_        
           #if not line:
             #break
           #else:   
	     #if(line.find(searchname_) >= 0):
	        #print "begin reading "+searchname_ 
                #line = File.readline()
                #line = File.readline()
                #j=0
                #if(line.find(str('Values')) >= 0):
		   #go_on = True 
                   #while go_on:
                     #line = File.readline()
                     #if not line:
                        #break
                   
                     #if (line.find("End values") >= 0):
                        #print 'END Reading elements'
                        #go_on = False;
                        
                        ##searchname_=str('Result "')+Quantity_+str ('" "Kratos" ')+str(t_)+str(' ')+Dimension_+str(' OnNodes')
                        ##print searchname_
                     #else:
		        
		        #aaa = line.split();
		        ##in_[j]=float(aaa[0])
		        #j=j+1
		        ##fluid_model_part.Nodes[j].SetSolutionStepValue(VELOCITY_X,0,float(aaa[1]))
###		        fluid_model_part.Nodes[j].SetSolutionStepValue(VELOCITY_Y,0,float(aaa[1]))
		        #fluid_model_part.Nodes[j].SetValue(VELOCITY_Y,float(aaa[1]))
		       
		   #t_=t_+1               
      #File.close()
      File = open(infile,'r')
      t_=1
      print "begin reading Intial PRESSURE values"  
      #i=0
      while (t_<=Number_Nodes):
	   searchname_=str('Result "PRESSURE_avrged" "Kratos" 1 Scalar OnNodes')
           line = File.readline()
           #print "searching for:   " +searchname_        
           if not line:
             break
           else:   
	     if(line.find(searchname_) >= 0):
	        print "begin reading "+searchname_ 
                line = File.readline()
                line = File.readline()
                j=0
                if(line.find(str('Values')) >= 0):
		   go_on = True 
                   while go_on:
                     line = File.readline()
                     if not line:
                        break
                   
                     if (line.find("End values") >= 0):
                        print 'END Reading elements'
                        go_on = False;
                        
                        #searchname_=str('Result "')+Quantity_+str ('" "Kratos" ')+str(t_)+str(' ')+Dimension_+str(' OnNodes')
                        #print searchname_
                     else:
		        
		        aaa = line.split();
		        #in_[j]=float(aaa[0])
		        j=j+1
		        #fluid_model_part.Nodes[j].SetSolutionStepValue(VELOCITY_X,0,float(aaa[1]))
##		        fluid_model_part.Nodes[j].SetSolutionStepValue(VELOCITY_Y,0,float(aaa[1]))
		        fluid_model_part.Nodes[j].SetValue(PRESSURE,float(aaa[1]))
		       
		   t_=t_+1               
      File.close()


# fill the PODS into the object node (just a lo0op) check pod application CPP
a = Vector(3)  #depending on the number of MOdes i take
a[0]=1
a[1] = 2
a[2] = 3
def Reduced_Pods_Read(infile,Variable,Number_Nodes,Number_Pods,BC_name):
      
      File = open(infile,'r')
      #IN=[]
      in_=Vector(Number_Pods)
      in_boundary=Vector(Number_Pods)
      t_=1
      print "begin reading Reduced PODS File"  
      #i=0
      while (t_<=Number_Nodes):
	   searchname_=str('Reduced POD values for    ')+str(t_)
           line = File.readline()
           #print "searching for:   " +searchname_        
           if not line:
             break
           else:   
	     if(line.find(searchname_) >= 0):
	        print "begin reading "+searchname_ 
                line = File.readline()
                
                j=0
                if(line.find(str('Values')) >= 0):
		   go_on = True 
                   while go_on:
                     line = File.readline()
                     if not line:
                        break
                   
                     if (line.find("End Values") >= 0):
                        print 'END Reading elements'
                        go_on = False;
                        
                        #searchname_=str('Result "')+Quantity_+str ('" "Kratos" ')+str(t_)+str(' ')+Dimension_+str(' OnNodes')
                        #print searchname_
                     else:
		        
		        aaa = line.split();
		        in_[j]=float(aaa[0])
		        j=j+1
		   
		   if(BC_name==False):
		     fluid_model_part.Nodes[t_].SetValue(Variable,in_[0:j])
		   #else:
	           
	             #if(fluid_model_part.Nodes[t_].IsFixed(BC_name) == False):
		        #fluid_model_part.Nodes[t_].SetValue(Variable,in_[0:j])
		     #else:
		        #fluid_model_part.Nodes[t_].SetValue(Variable,in_boundary[0:j])
		   
		   t_=t_+1
	               
      File.close()
      
     
 #POD_VELOCITY_X    
npods = 20
nnodes = len(fluid_model_part.Nodes)


dirname = './PODS_verticies_combined/'
#dirname = './20sec_200shots_4podmodes_average.intial/'

Intial_Values_Read(dirname+'averaged_snapshots.res',nnodes)
#Intial_Values_Read('cylinder_0.post.res',nnodes)

Reduced_Pods_Read(dirname+'VH_Pressure_reduced.txt',POD_PRESSURE,nnodes,npods,False)
Reduced_Pods_Read(dirname+'VH_VY_reduced.txt',POD_VELOCITY_Y,nnodes,npods,False)
Reduced_Pods_Read(dirname+'VH_VX_reduced.txt',POD_VELOCITY_X,nnodes,npods,False)


#Reduced_Pods_Read('VH_VZ_reduced.txt',POD_VELOCITY_Z,nnodes,npods)

#Reduced_Pods_Read('VH_V_all_reduced.txt',POD_VELOCITY_X,967,5)
#Reduced_Pods_Read('VH_V_all_reduced.txt',POD_VELOCITY_Y,967,5)
#Reduced_Pods_Read('VH_V_all_reduced.txt',POD_VELOCITY_Z,967,5)


#print fluid_model_part.Nodes[30].GetValue(POD_PRESSURE)
#print fluid_model_part.Nodes[30].GetValue(POD_VELOCITY_X)
#print fluid_model_part.Nodes[30].GetValue(POD_VELOCITY_Y)
#print fluid_model_part.Nodes[30].GetValue(POD_VELOCITY_Z)

#print "sdfgjklsdfgljdsfgkdklñsfgdsfjjsdflkñgjklñdfjsglñkjsdlñkgjklñdfsjgldfjslgjsdflñjglñfdjslgñjdslñfjglñsdfjglñfdsjlñgjñsdlfg"

#creating the solvers
#fluid solver
if(SolverType == "FractionalStep"):
    fluid_solver = incompressible_fluid_solver.IncompressibleFluidSolver(fluid_model_part,domain_size)
    fluid_solver.laplacian_form = laplacian_form; #standard laplacian form
    fluid_solver.predictor_corrector = fluid_only_var.predictor_corrector
    fluid_solver.max_press_its = fluid_only_var.max_press_its
    fluid_solver.Initialize()
elif(SolverType == "monolithic_solver_eulerian"): 
    fluid_solver = monolithic_solver_eulerian_standard_biggerlocalM.MonolithicSolver(fluid_model_part,domain_size)
    oss_swith = fluid_only_var.use_oss
    dynamic_tau = fluid_only_var.dynamic_tau
    fluid_solver.max_iter = 30
    fluid_model_part.ProcessInfo.SetValue(OSS_SWITCH, oss_swith);				
    fluid_model_part.ProcessInfo.SetValue(DYNAMIC_TAU, dynamic_tau);
    fluid_solver.Initialize()
elif(SolverType == "monolithic_solver_eulerian_compressible"): 
    fluid_solver = monolithic_solver_eulerian_compressible.MonolithicSolver(fluid_model_part,domain_size)
    oss_swith = fluid_only_var.use_oss
    dynamic_tau = fluid_only_var.dynamic_tau
    fluid_model_part.ProcessInfo.SetValue(OSS_SWITCH, oss_swith);				
    fluid_model_part.ProcessInfo.SetValue(DYNAMIC_TAU, dynamic_tau);
    fluid_solver.Initialize()


#print "fluid solver created"

#settings to be changed
Dt = fluid_only_var.Dt 
full_Dt = Dt 
initial_Dt = 1.0 * full_Dt #0.05 #0.01
final_time = fluid_only_var.max_time
output_step = fluid_only_var.output_step

fluid_only_var.output_step

out = 0

#x = raw_input("stopped to allow debug: set breakpoints and press enter to continue");

#mesh to be printed
mesh_name = 0.0
gid_io.InitializeMesh( mesh_name)
gid_io.WriteMesh( fluid_model_part.GetMesh() )
gid_io.FinalizeMesh()



gid_io.InitializeResults(mesh_name,(fluid_model_part).GetMesh())

for node in fluid_model_part.Nodes:
  vel = node.GetValue(VELOCITY)
  node.SetSolutionStepValue(VELOCITY,0,vel)
  
  press = node.GetValue(PRESSURE)
  node.SetSolutionStepValue(PRESSURE,0,press)

  print node.GetValue(POD_VELOCITY_X)
  print node.GetValue(POD_VELOCITY_Y)
  print node.GetValue(POD_PRESSURE)


time = 0.0
step = 0
while(time < final_time):
  
    
    #fluid_model_part.Nodes[587].GetValue(VELOCITY)
    if(step < 5):
        Dt = initial_Dt
    else:
        Dt = full_Dt
        
    time = time + Dt
    fluid_model_part.CloneTimeStep(time)

    if(step >= 3):
        #if(step == 3):
            #print "*****************************"
            #for node in fluid_model_part.Nodes:
                #vel = node.GetValue(VELOCITY)
                #node.SetSolutionStepValue(VELOCITY,0,vel)
	#for node in fluid_model_part.Nodes:
	  #vel = node.GetValue(VELOCITY)
	  #node.SetSolutionStepValue(VELOCITY,0,vel)
        
	#print VELOCITY
        fluid_solver.Solve()
        
        #if(time < 20.0):
	  #for node in fluid_model_part.Nodes:
	    #node.SetSolutionStepValue(BODY_FORCE_Y,0,-1.0)
	#else:
	  #for node in fluid_model_part.Nodes:
	    #node.SetSolutionStepValue(BODY_FORCE_Y,0,0.0)
        
    

    if(out == output_step):
        if(SolverType == "FractionalStep"):
            gid_io.WriteNodalResults(PRESSURE,fluid_model_part.Nodes,time,0)
            gid_io.WriteNodalResults(VELOCITY,fluid_model_part.Nodes,time,0)
        else:
            gid_io.WriteNodalResults(PRESSURE,fluid_model_part.Nodes,time,0)
            gid_io.WriteNodalResults(VELOCITY,fluid_model_part.Nodes,time,0)
            gid_io.WriteNodalResults(REACTION,fluid_model_part.Nodes,time,0)
 
        out = 0

    out = out + 1
    step = step + 1
      
gid_io.FinalizeResults()
          
 


