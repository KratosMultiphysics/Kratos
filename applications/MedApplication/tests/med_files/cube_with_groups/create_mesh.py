#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.6.0 with dump python functionality
###

import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()

###
### GEOM component
###

import GEOM
from salome.geom import geomBuilder
import math
import SALOMEDS


geompy = geomBuilder.New()

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)
Box_1 = geompy.MakeBoxDXDYDZ(200, 200, 200)
[Face_1] = geompy.SubShapes(Box_1, [13])
[Face_1] = geompy.GetExistingSubObjects(Box_1, False)
[Face_1] = geompy.GetExistingSubObjects(Box_1, False)
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Box_1, 'Box_1' )
geompy.addToStudyInFather( Box_1, Face_1, 'Face_1' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

#hyp_0.SetLength( 100 ) ### not created Object
NETGEN_3D_Parameters_1 = smesh.CreateHypothesisByAverageLength( 'NETGEN_Parameters', 'NETGENEngine', 100, 0 )
NETGEN_1D_2D_3D = smesh.CreateHypothesis('NETGEN_2D3D', 'NETGENEngine')

#hyp_3.SetLength( 50 ) ### not created Object
NETGEN_3D_Parameters_2 = smesh.CreateHypothesisByAverageLength( 'NETGEN_Parameters', 'NETGENEngine', 50, 0 )
cube = smesh.Mesh(Box_1)
status = cube.AddHypothesis( Box_1, NETGEN_3D_Parameters_2 )
status = cube.AddHypothesis(NETGEN_1D_2D_3D)
interface = cube.GroupOnGeom(Face_1,'Face_1',SMESH.FACE)
isDone = cube.Compute()
[ interface ] = cube.GetGroups()
interface_nodes = cube.GroupOnGeom(Face_1,'Group_1',SMESH.NODE)
[ interface, interface_nodes ] = cube.GetGroups()
interface_nodes.SetName( 'interface_nodes' )
interface.SetName( 'interface' )
smesh.SetName(cube, 'cube')

## Set names of Mesh objects
smesh.SetName(NETGEN_1D_2D_3D, 'NETGEN 1D-2D-3D')
smesh.SetName(interface_nodes, 'interface_nodes')
smesh.SetName(NETGEN_3D_Parameters_2, 'NETGEN 3D Parameters_2')
smesh.SetName(NETGEN_3D_Parameters_1, 'NETGEN 3D Parameters_1')
smesh.SetName(cube.GetMesh(), 'cube')
smesh.SetName(interface, 'interface')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
