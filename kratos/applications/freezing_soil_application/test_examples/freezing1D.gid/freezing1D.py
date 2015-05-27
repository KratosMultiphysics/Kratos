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
sys.path.append('./freezing1D.gid')
import freezing1D_include
from freezing1D_include import *
# calculate insitu-stress for geology_virgin.gid
model = freezing1D_include.Model('freezing1D','./')
model.InitializeModel()

##################################################################
###  SIMULATION  #################################################
##################################################################    
if(1==1):
      watchElemNum= 100
      model.model_part.Elements[watchElemNum].SetValue(KRATOS_WATCH_FLAG, 1)

model.model_part.Elements[1].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[1].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[2].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[2].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[3].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[3].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[4].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[4].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[5].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[5].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[6].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[6].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[7].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[7].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[8].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[8].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[9].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[9].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[10].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[10].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[11].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[11].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[12].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[12].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[13].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[13].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[14].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[14].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[15].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[15].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[16].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[16].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[17].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[17].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[18].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[18].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[19].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[19].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[20].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[20].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[21].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[21].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[22].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[22].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[23].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[23].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[24].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[24].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[25].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[25].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[26].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[26].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[27].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[27].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[28].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[28].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[29].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[29].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[30].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[30].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[31].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[31].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[32].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[32].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[33].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[33].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[34].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[34].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[35].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[35].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[36].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[36].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[37].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[37].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[38].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[38].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[39].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[39].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[40].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[40].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[41].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[41].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[42].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[42].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[43].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[43].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[44].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[44].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[45].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[45].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[46].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[46].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[47].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[47].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[48].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[48].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[49].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[49].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[50].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[50].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[51].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[51].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[52].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[52].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[53].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[53].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[54].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[54].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[55].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[55].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[56].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[56].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[57].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[57].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[58].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[58].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[59].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[59].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[60].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[60].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[61].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[61].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[62].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[62].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[63].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[63].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[64].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[64].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[65].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[65].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[66].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[66].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[67].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[67].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[68].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[68].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[69].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[69].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[70].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[70].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[71].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[71].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[72].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[72].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[73].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[73].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[74].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[74].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[75].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[75].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[76].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[76].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[77].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[77].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[78].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[78].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[79].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[79].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[80].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[80].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[81].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[81].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[82].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[82].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[83].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[83].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[84].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[84].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[85].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[85].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[86].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[86].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[87].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[87].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[88].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[88].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[89].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[89].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[90].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[90].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[91].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[91].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[92].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[92].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[93].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[93].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[94].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[94].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[95].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[95].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[96].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[96].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[97].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[97].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[98].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[98].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[99].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[99].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
model.model_part.Elements[100].SetValue(PLASTIC_FLAG, 1)  
model.model_part.Elements[100].SetValue(ASSIGN_PRESTRESS_FLAG, 0)
FlagPrestress = 0

print ("++++++++++ << Start Custom Simulation: "+'freezing1D.gid >> ++++++++++') 
# =====================
# | USER SCRIPT FOR CALCULATION OF AUSBLAS.GID |
# vvvvvvvvvvvvvvvvvvvvv  

time= 0.0 

Ti= 2.0
Tx= -8.0  

print ("############## Initial condition: Ti= "+str(Ti)+" ############## ")
for i in range(1, len(model.model_part.Nodes)+1):     
	#model.model_part.Nodes[i].SetSolutionStepValue(DISPLACEMENT_X, 0.0)   
	#model.model_part.Nodes[i].SetSolutionStepValue(DISPLACEMENT_Y, 0.0)   
	#model.model_part.Nodes[i].SetSolutionStepValue(DISPLACEMENT_Z, 0.0)   
	#model.model_part.Nodes[i].Fix(DISPLACEMENT_X)   
	#model.model_part.Nodes[i].Fix(DISPLACEMENT_Y)   
	#model.model_part.Nodes[i].Fix(DISPLACEMENT_Z)     
	model.model_part.Nodes[i].SetSolutionStepValue(WATER_PRESSURE, 0)    
	#model.model_part.Nodes[i].Fix(WATER_PRESSURE)    
	model.model_part.Nodes[i].SetSolutionStepValue(TEMPERATURE, Ti)   

for i in range(0, len(model.layer_nodes_sets['topSurface'])):
	model.model_part.Nodes[model.layer_nodes_sets['topSurface'][i]].SetSolutionStepValue(CONVECTION_COEFFICIENT, 100)  
	model.model_part.Nodes[model.layer_nodes_sets['topSurface'][i]].Fix(CONVECTION_COEFFICIENT)    
	model.model_part.Nodes[model.layer_nodes_sets['topSurface'][i]].SetSolutionStepValue(AMBIENT_TEMPERATURE, Tx)  
	model.model_part.Nodes[model.layer_nodes_sets['topSurface'][i]].Fix(AMBIENT_TEMPERATURE)   
	#model.model_part.Nodes[model.layer_nodes_sets['topSurface'][i]].SetSolutionStepValue(TEMPERATURE, Tx)  
	#model.model_part.Nodes[model.layer_nodes_sets['topSurface'][i]].Fix(TEMPERATURE)      
for i in range(0, len(model.layer_nodes_sets['bottomSurface'])): 
	model.model_part.Nodes[model.layer_nodes_sets['bottomSurface'][i]].Fix(WATER_PRESSURE)   
	model.model_part.Nodes[model.layer_nodes_sets['bottomSurface'][i]].Fix(TEMPERATURE)   	
	model.model_part.Nodes[model.layer_nodes_sets['bottomSurface'][i]].SetSolutionStepValue(DISPLACEMENT_X, 0.0)   
	model.model_part.Nodes[model.layer_nodes_sets['bottomSurface'][i]].SetSolutionStepValue(DISPLACEMENT_Y, 0.0)   
	model.model_part.Nodes[model.layer_nodes_sets['bottomSurface'][i]].SetSolutionStepValue(DISPLACEMENT_Z, 0.0)   
	model.model_part.Nodes[model.layer_nodes_sets['bottomSurface'][i]].Fix(DISPLACEMENT_X)   
	model.model_part.Nodes[model.layer_nodes_sets['bottomSurface'][i]].Fix(DISPLACEMENT_Y)   
	model.model_part.Nodes[model.layer_nodes_sets['bottomSurface'][i]].Fix(DISPLACEMENT_Z)     
 

print ("############## Transient Phase: Set topSurface Tx= "+str(Tx)+" ############## ")
listSteps= [30,6,2, 360,3,3, 1200,3,5, 3600,3,6]  # 1h + 5h +  18h = 1d  
print "############## >> 1st Hour:" 
## [Dt, n, inner_n];  
listSteps= [1200,3,6, 3600,7,6]  # 1h + 5h +  18h = 1d  
listSteps= [10,20,2,20,10,1]  # 1h   
listSteps= [10,20,2, 50,10,1]  
listSteps= [10,10,1, 50,10,1, 100,10,1]  
listSteps= [10,10,2, 50,10,1, 100,10,1]  
listSteps= [30,6,3, 60,6,3, 120,6,3, 360,6,3, 1200,3,5, 3600,7,6 ]  # 1h + 5h +  18h = 2d  
 
listSteps= [5,10,1, 30,6,3]  # 1h + 5h +  18h = 2d  
listSteps= [20,20,1]  # 1h + 5h +  18h = 2d  
listSteps= [30,6,3, 60,6,3, 120,6,3, 600,6,3, 1200,6,3, 3600,6,3]   
listSteps= [30,6,2, 60,6,3, 360,3,2, 1200,3,5, 3600,7,6]  # 1h + 5h +  18h = 2d  

for step in range(0, len(listSteps)/3):
	Dt= listSteps[step*3]
	n= listSteps[step*3+1]
	inner_n= listSteps[step*3+2]
	for step in range(1,n+1):
		for inner_step in range(1,inner_n+1):
			time = time + Dt
			model.Solve( time )
		model.WriteOutput(time) 
		print("~~~~~~~~~~~~~~ ( Dt= "+str(Dt)+", step= "+str(step)+", time= "+str(time)+" s = "+str(time/60.0)+" m ) ~~~~~~~~~~~~~~") 
 
for i in range(0, len(model.layer_nodes_sets['topSurface'])):
	model.model_part.Nodes[model.layer_nodes_sets['topSurface'][i]].SetSolutionStepValue(TEMPERATURE, Tx)  
	model.model_part.Nodes[model.layer_nodes_sets['topSurface'][i]].SetSolutionStepValue(TEMPERATURE_NULL, Tx)  
	model.model_part.Nodes[model.layer_nodes_sets['topSurface'][i]].SetSolutionStepValue(TEMPERATURE_EINS, Tx)  
	model.model_part.Nodes[model.layer_nodes_sets['topSurface'][i]].Fix(TEMPERATURE)  
  
#print "############## >> 2st Day: " 
listSteps= [3600,8,6, 7200,12,6, 14400,10,6, 43200,10,6, 86400,30,6 ]

for step in range(0, len(listSteps)/3):
	Dt= listSteps[step*3]
	n= listSteps[step*3+1]
	inner_n= listSteps[step*3+2]
	for step in range(1,n+1):
		for inner_step in range(1,inner_n+1):
			time = time + Dt
			model.Solve( time )
		model.WriteOutput(time)  
		print("~~~~~~~~~~~~~~ ( Dt= "+str(Dt)+", step= "+str(step)+", time= "+str(time)+" s = "+str(time/3600.0/24)+" d ) ~~~~~~~~~~~~~~") 
sys.exit(0) 
  