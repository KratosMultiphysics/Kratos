from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
domain_size = 3

Dt = 0.1
max_time = 0.1
nsteps = 10
DestinationMesh = r'/home/jcotela/Desktop/destination2.gid/destination2.mdpa'
DestinationMesh = DestinationMesh[0:DestinationMesh.rindex('.mdpa')]
DestinationMesh = DestinationMesh.replace('/', '//')
# Declare Python Variables

problem_name = "meshtest"
problem_path = "/home/jcotela/Desktop/meshtest.gid"
kratos_path = "/home/jcotela/kratos"
