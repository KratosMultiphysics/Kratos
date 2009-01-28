##################################################################################
##################################################################################
###########  PYTHON SKRIPT ZUM ERSTELLEN EINES COUPLING PROBLEM TYPES  ###########
##################################################################################
##################################################################################

import pt_generation;
import filegeneration;

##################################################################################
#############       Projectname            #######################################
##################################################################################

projectname = 'structural_application'
projectpath=''

##only Fluid
#mode = 1
#
#only Structure
mode = 2
#
#Fluid-Structure-Interaction
#mode=3

##################################################################################
#####        Filegeneration        ###############################################
##################################################################################

filegeneration.GenerateAllFiles(projectname,projectpath,mode)

##################################################################################
####     Boundary - Conditions    ################################################
##################################################################################
###p=point,l=line, s=surface, v=volume ... combinations plsv ####

### FLUID ####

### STRUCTURE ###
pt_generation.AddScalarBC('PRESSURE','plvs','nodes',projectname,'Structure') 
pt_generation.AddScalarBC('NEGATIVE_FACE_PRESSURE','plsv','nodes',projectname,'Structure')
pt_generation.AddScalarBC('POSITIVE_FACE_PRESSURE','plsv','nodes',projectname,'Structure')
pt_generation.AddScalarBC('IS_INTERFACE','plsv','nodes',projectname,'Structure')
pt_generation.AddScalarBC('FLAG_VARIABLE','plsv','nodes',projectname,'Structure')

pt_generation.AddVecBC('DISPLACEMENT','plsv','nodes',projectname,'Structure')
pt_generation.AddVecBC('ROTATION','plsv','nodes',projectname,'Structure')
pt_generation.AddVecBC('VELOCITY','plsv','nodes',projectname,'Structure')
pt_generation.AddVecBC('DAMAGE','plsv','nodes',projectname,'Structure')

##pt_generation.AddScalarBC('NEGATIVE_FACE_PRESSURE','ls','nodes',projectname,'Structure')
##pt_generation.AddVecBC('DISPLACEMENT','plsv','nodes',projectname,'Structure')

#######################################################################################
####    Conditions         ############################################################
#######################################################################################

### FLUID ####

##
##### STRUCTURE ###
pt_generation.AddCondition('Face2D','l',projectname,'Structure')
pt_generation.AddCondition('Face3D3N','s',projectname,'Structure')
pt_generation.AddCondition('Face3D4N','s',projectname,'Structure')
pt_generation.AddCondition('Face3D9N','s',projectname,'Structure')
pt_generation.AddCondition('Condition2D2N','l',projectname,'Structure')
pt_generation.AddCondition('Condition3D3N','s',projectname,'Structure')


#######################################################################################
####   Elements   #####################################################################
#######################################################################################

### FLUID ####

### STRUCTURE ###
pt_generation.AddElement('TotalLagrangian2D3N','s',projectname,'Structure')
pt_generation.AddElement('TotalLagrangian2D6N','s',projectname,'Structure')
pt_generation.AddElement('TotalLagrangian2D4N','s',projectname,'Structure')
pt_generation.AddElement('TotalLagrangian2D8N','s',projectname,'Structure')
pt_generation.AddElement('TotalLagrangian3D4N','v',projectname,'Structure')
pt_generation.AddElement('TotalLagrangian3D8N','v',projectname,'Structure')
pt_generation.AddElement('TotalLagrangian3D27N','v',projectname,'Structure')
pt_generation.AddElement('MembraneElement','s',projectname,'Structure')
pt_generation.AddElement('IsoShellElement','s',projectname,'Structure')
pt_generation.AddElement('AnisoShellElement','s',projectname,'Structure')
pt_generation.AddElement('AnisoLinearShellElement','s',projectname,'Structure')


########################################################################################
########    Materials    ###############################################################
########################################################################################


#pt_generation.AddFluidMaterial('Air','1.2','0.000017','0',projectname)
#pt_generation.AddFluidMaterial('Fluid','1','1','0',projectname)
pt_generation.AddStructureMaterial('Aluminium','2700','70000','0.3','1','1',projectname)
pt_generation.AddStructureMaterial('StrucMat','2700','70000','0.3','1','1',projectname)
