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

projectname = 'incompressible_fluid_application'
##projectname = 'ALEapplication'
projectpath=''
#projectpath = '..\\neu\\'

##only Fluid
mode = 1
#
#only Structure
#mode = 2
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
pt_generation.AddScalarBC('PRESSURE','plvs','nodes',projectname,'Fluid') 
pt_generation.AddScalarBC('VISCOSITY','plsv','nodes',projectname,'Fluid')
pt_generation.AddScalarBC('DENSITY','plsv','nodes',projectname,'Fluid')
pt_generation.AddScalarBC('TEMPERATURE','plsv','nodes',projectname,'Fluid')
pt_generation.AddScalarBC('CONDUCTIVITY','plsv','nodes',projectname,'Fluid')
pt_generation.AddScalarBC('SPECIFIC_HEAT','plsv','nodes',projectname,'Fluid')
pt_generation.AddScalarBC('HEAT_FLUX','plsv','nodes',projectname,'Fluid')
pt_generation.AddScalarBC('IS_POROUS','plsv','nodes',projectname,'Fluid')
pt_generation.AddScalarBC('POROSITY','plsv','nodes',projectname,'Fluid')
pt_generation.AddScalarBC('EXTERNAL_PRESSURE','plsv','nodes',projectname,'Fluid')

pt_generation.AddScalarBC('IS_STRUCTURE','plsv','nodes',projectname,'Fluid')
pt_generation.AddScalarBC('IS_BOUNDARY','plsv','nodes',projectname,'Fluid')
pt_generation.AddScalarBC('IS_INTERFACE','plsv','nodes',projectname,'Fluid')
pt_generation.AddScalarBC('FLAG_VARIABLE','plsv','nodes',projectname,'Fluid')

pt_generation.AddVecBC('DISPLACEMENT','plsv','nodes',projectname,'Fluid')
pt_generation.AddVecBC('VELOCITY','plsv','nodes',projectname,'Fluid')
pt_generation.AddVecBC('BODY_FORCE','plsv','nodes',projectname,'Fluid')


### STRUCTURE ###
##pt_generation.AddScalarBC('NEGATIVE_FACE_PRESSURE','ls','nodes',projectname,'Structure')
##pt_generation.AddVecBC('DISPLACEMENT','plsv','nodes',projectname,'Structure')

#######################################################################################
####    Conditions         ############################################################
#######################################################################################

### FLUID ####
pt_generation.AddCondition('Condition2D','l',projectname,'Fluid')
pt_generation.AddCondition('Condition3D','s',projectname,'Fluid')
pt_generation.AddCondition('Fluid2DNeumann','l',projectname,'Fluid')
pt_generation.AddCondition('Fluid3DNeumann','s',projectname,'Fluid')

##
##### STRUCTURE ###
##pt_generation.AddCondition('Face2D','l',projectname,'Structure')
##pt_generation.AddCondition('Face3D','s',projectname,'Structure')


#######################################################################################
####   Elements   #####################################################################
#######################################################################################

### FLUID ####
pt_generation.AddElement('Fluid2D','s',projectname,'Fluid')
pt_generation.AddElement('NDFluid2D','s',projectname,'Fluid')
pt_generation.AddElement('NDFluid2DCrankNicolson','s',projectname,'Fluid')
pt_generation.AddElement('Fluid3D','v',projectname,'Fluid')
pt_generation.AddElement('Fluid2DCoupled','s',projectname,'Fluid')
pt_generation.AddElement('Fluid3DCoupled','v',projectname,'Fluid')
pt_generation.AddElement('ConvDiff2D','s',projectname,'Fluid')
pt_generation.AddElement('ConvDiff3D','v',projectname,'Fluid')

### STRUCTURE ###
##pt_generation.AddElement('TotalLagrangian','sv',projectname,'Structure')


########################################################################################
########    Materials    ###############################################################
########################################################################################


pt_generation.AddFluidMaterial('Air','1.2','0.000017','0',projectname)
#pt_generation.AddFluidMaterial('Fluid','1','1','0',projectname)
#pt_generation.AddStructureMaterial('Aluminium','2700','70000','0.3','1','1',projectname)
#pt_generation.AddStructureMaterial('StrucMat','2700','70000','0.3','1','1',projectname)
