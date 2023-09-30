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
[Face_1,Face_2] = geompy.SubShapes(Box_1, [13, 23])
[Edge_1] = geompy.SubShapes(Box_1, [15])
[Face_3,Face_4] = geompy.SubShapes(Box_1, [33, 31])
[Face_1, Face_2, Edge_1] = geompy.GetExistingSubObjects(Box_1, False)
Auto_group_for_Group_1 = geompy.CreateGroup(Box_1, geompy.ShapeType["FACE"])
geompy.UnionList(Auto_group_for_Group_1, [Face_3, Face_4])
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Box_1, 'Box_1' )
geompy.addToStudyInFather( Box_1, Face_1, 'Face_1' )
geompy.addToStudyInFather( Box_1, Face_2, 'Face_2' )
geompy.addToStudyInFather( Box_1, Edge_1, 'Edge_1' )
geompy.addToStudyInFather( Box_1, Face_3, 'Face_3' )
geompy.addToStudyInFather( Box_1, Face_4, 'Face_4' )
geompy.addToStudyInFather( Box_1, Auto_group_for_Group_1, 'Auto_group_for_Group_1' )

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
[ Face_1_1, Face_2_1, Edge_1_1 ] = cube.GetGroups()
isDone = cube.Compute()
[ Face_1_1, Face_2_1, Edge_1_1 ] = cube.GetGroups()
isDone = cube.Compute()
smesh.SetName(cube, 'cube')
[ Face_1_1, Face_2_1, Edge_1_1 ] = cube.GetGroups()
boundary = cube.GroupOnGeom(Auto_group_for_Group_1,'Group_1',SMESH.FACE)
[ Face_1_1, Face_2_1, Edge_1_1, boundary ] = cube.GetGroups()
boundary.SetName( 'boundary' )
[ Face_1_1, Face_2_1, Edge_1_1, boundary ] = cube.GetGroups()
boundary_top = cube.GroupOnGeom(Face_3,'Group_1',SMESH.FACE)
[ Face_1_1, Face_2_1, Edge_1_1, boundary, boundary_top ] = cube.GetGroups()
boundary_top.SetName( 'boundary.top' )
[ Face_1_1, Face_2_1, Edge_1_1, boundary, boundary_top ] = cube.GetGroups()
boundary_bottom = cube.GroupOnGeom(Face_4,'Group_1',SMESH.FACE)
[ Face_1_1, Face_2_1, Edge_1_1, boundary, boundary_top, boundary_bottom ] = cube.GetGroups()
boundary_bottom.SetName( 'boundary.bottom' )
isDone = cube.Compute()
[ Face_1_1, Face_2_1, Edge_1_1, boundary, boundary_top, boundary_bottom ] = cube.GetGroups()
smesh.SetName(cube, 'cube')

## Set names of Mesh objects
smesh.SetName(NETGEN_1D_2D_3D.GetAlgorithm(), 'NETGEN 1D-2D-3D')
smesh.SetName(NETGEN_3D_Parameters_1, 'NETGEN 3D Parameters_1')
smesh.SetName(Face_1_1, 'Face_1')
smesh.SetName(Face_2_1, 'Face_2')
smesh.SetName(boundary, 'boundary')
smesh.SetName(boundary_top, 'boundary.top')
smesh.SetName(boundary_bottom, 'boundary.bottom')
smesh.SetName(cube.GetMesh(), 'cube')
smesh.SetName(Edge_1_1, 'Edge_1')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
