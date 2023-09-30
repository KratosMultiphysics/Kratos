#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.6.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, r'/home/philipp/software/Kratos_master/applications/MedApplication/tests/med_files/cube_with_adjacent_groups')

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
[Face_1,Face_2] = geompy.SubShapes(Box_1, [13, 23])
[Edge_1] = geompy.SubShapes(Box_1, [15])
[Face_1, Face_2, Edge_1] = geompy.GetExistingSubObjects(Box_1, False)
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Box_1, 'Box_1' )
geompy.addToStudyInFather( Box_1, Face_1, 'Face_1' )
geompy.addToStudyInFather( Box_1, Face_2, 'Face_2' )
geompy.addToStudyInFather( Box_1, Edge_1, 'Edge_1' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

#hyp_3.SetLength( 75 ) ### not created Object
NETGEN_3D_Parameters_1 = smesh.CreateHypothesisByAverageLength( 'NETGEN_Parameters', 'NETGENEngine', 75, 0 )
cube = smesh.Mesh(Box_1)
status = cube.AddHypothesis( Box_1, NETGEN_3D_Parameters_1 )
NETGEN_1D_2D_3D = cube.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D)
Face_1_1 = cube.GroupOnGeom(Face_1,'Face_1',SMESH.FACE)
Face_2_1 = cube.GroupOnGeom(Face_2,'Face_2',SMESH.FACE)
Edge_1_1 = cube.GroupOnGeom(Edge_1,'Edge_1',SMESH.EDGE)
isDone = cube.Compute()
[ Face_1_1, Face_2_1, Edge_1_1 ] = cube.GetGroups()
smesh.SetName(cube, 'cube')
try:
  cube.ExportMED(r'/home/philipp/software/Kratos_master/applications/MedApplication/tests/med_files/cube_with_adjacent_groups/cube.med',auto_groups=0,version=40,overwrite=1,meshPart=None,autoDimension=1)
  pass
except:
  print('ExportMED() failed. Invalid file name?')


## Set names of Mesh objects
smesh.SetName(NETGEN_1D_2D_3D.GetAlgorithm(), 'NETGEN 1D-2D-3D')
smesh.SetName(NETGEN_3D_Parameters_1, 'NETGEN 3D Parameters_1')
smesh.SetName(Face_1_1, 'Face_1')
smesh.SetName(Face_2_1, 'Face_2')
smesh.SetName(cube.GetMesh(), 'cube')
smesh.SetName(Edge_1_1, 'Edge_1')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
