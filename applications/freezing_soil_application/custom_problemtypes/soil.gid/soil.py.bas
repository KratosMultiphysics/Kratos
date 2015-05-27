##################################################################
##### ekate - Enhanced KRATOS for Advanced Tunnel Enineering #####
##### copyright by CIMNE, Barcelona, Spain                   #####
#####          and Janosch Stascheit for TUNCONSTRUCT        #####
##### all rights reserved                                    #####
##################################################################
#setting the domain size for the problem to be solved
domain_size = 3
##################################################################
##################################################################
## ATTENTION: here the order is important                    #####
##################################################################
## including kratos path                                     #####
## ATTENTION: the following lines have to be adapted to      #####
##            match your acrtual configuration               #####
##################################################################
import sys
import os
kratos_root_path=os.environ['KRATOS_ROOT_PATH'] 
##setting up paths
kratos_libs_path = kratos_root_path+'libs' ##kratos_root/libs
kratos_applications_path = kratos_root_path+'applications' ##kratos_root/applications
##################################################################
##################################################################
sys.path.append(kratos_libs_path)
sys.path.append(kratos_applications_path)

##################################################################
##################################################################
sys.path.append('./rEpLaCeMeNtStRiNg.gid')
import rEpLaCeMeNtStRiNg_include
from rEpLaCeMeNtStRiNg_include import **
# calculate insitu-stress for geology_virgin.gid
model = rEpLaCeMeNtStRiNg_include.Model('rEpLaCeMeNtStRiNg','./')
model.InitializeModel()

##################################################################
###  SIMULATION  #################################################
##################################################################    
if(*GenData(Enable_Kratos_Watch,int)==1):
      watchElemNum= *GenData(Watch_Element_Number,int)
      model.model_part.Elements[watchElemNum].SetValue(KRATOS_WATCH_FLAG, 1)

*set var plasticflag(int)=GenData(Plastic_Analysis,int)
*set var prestressflag(int)=GenData(Assign_Prestress,int)
*loop elems
model.model_part.Elements[*ElemsNum].SetValue(PLASTIC_FLAG, *plasticflag)  
model.model_part.Elements[*ElemsNum].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
*end elems 
FlagPrestress = *prestressflag

*if(strcmp(GenData(Simulation_Script),"standard")==0)   
print ("++++++++++ << Start Standard Simulation: "+'rEpLaCeMeNtStRiNg.gid >> ++++++++++') 
Dt= *GenData(time_step_length,real) 
n= *GenData(time_steps,int) 
small_n= *GenData(small_time_steps,int) 
small_Dt= Dt/float(small_n)
time = small_Dt 
for small_step in range(1,small_n):	
	model.Solve(time)
	if(*GenData(show_small_time,int)==1):
		model.WriteOutput(time)
	print("~~~~~~~~~~~~~~ ( time= "+str(time)+": Dt= "+str(small_Dt)+" / smallStep= "+str(small_step)+" ) ~~~~~~~~~~~~~~") 
	time = time+small_Dt
for step in range(1,n+1):	
	model.Solve(time)
	model.WriteOutput(time)
	print("~~~~~~~~~~~~~~ ( time= "+str(time)+": Dt= "+str(Dt)+" / step= "+str(step)+" ) ~~~~~~~~~~~~~~") 
	time = time+Dt
if(FlagPrestress==1):
	print ("+++++++ Assign Prestress +++++++") 
*loop elems
	model.model_part.Elements[*ElemsNum].SetValue(ASSIGN_PRESTRESS_FLAG, FlagPrestress)
*end elems 
	stable_Dt= *GenData(stable_time_step_length,real)
	stable_n= *GenData(stable_time_steps,int) 
	stable_small_n= *GenData(stable_small_time_steps,int) 
	stable_small_Dt= stable_Dt/float(stable_small_n)
	time = time+stable_small_Dt
	for stable_small_step in range(1,stable_small_n):	
		model.Solve(time)
		if(*GenData(show_stable_small_time,int)==1):
			model.WriteOutput(time)
		print("~~~~~~~~~~~~~~ ( time= "+str(time)+": stable_Dt= "+str(stable_small_Dt)+" / stable_step= "+str(stable_small_step)+" ) ~~~~~~~~~~~~~~") 
		time = time+stable_small_Dt
	for stable_step in range(1,stable_n+1):	
		model.Solve(time)
		model.WriteOutput(time)
		print("~~~~~~~~~~~~~~ ( time= "+str(time)+": stable_Dt= "+str(stable_Dt)+" / stable_step= "+str(stable_step)+" ) ~~~~~~~~~~~~~~") 
		time = time+stable_Dt
print "++++++++++ << Standard Simulation Done! >> ++++++++++" 
sys.exit(0)
*else 
print ("++++++++++ << Start Custom Simulation: "+'rEpLaCeMeNtStRiNg.gid >> ++++++++++') 
*# user-defined script is used (will be appended automatically)
# =====================
# | USER SCRIPT FOR CALCULATION OF AUSBLAS.GID |
# vvvvvvvvvvvvvvvvvvvvv  
*endif
