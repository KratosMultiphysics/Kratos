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
##################################################################
import sys
##################################################################
sys.path.append('./beam_3d.gid')
import beam_3d_include
from beam_3d_include import *
# calculate insitu-stress for geology_virgin.gid
model = beam_3d_include.Model('beam_3d','./')
linear_elastic = False
model.InitializeModel( linear_elastic )

##################################################################
###  SIMULATION  #################################################
##################################################################
time = 0.0
model.SetCalculateInSituStress( False )
model.Solve( time, 0, 0, 0, 0 )
model.WriteOutput( time )
print "Calculation done"
sys.exit(0)
