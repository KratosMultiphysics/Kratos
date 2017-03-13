from KratosMultiphysics import *
from KratosMultiphysics.NURBSBRepApplication import *
import KratosMultiphysics.ExternalSolversApplication 

import json

#cad_geometry_input_filename = 'geometry.json'
#with open(cad_geometry_input_filename) as cad_data1:
#    cad_geometry = json.load(cad_data1)

#import define_output
cad_geometry_file = open("geometry.json",'r')
cad_geometry = Parameters( cad_geometry_file.read())

model_part = ModelPart("NurbsCADGeometry")

geometry_reader = BrepModelGeometryReader(cad_geometry)

#brep_model_vector = []
#brep_model_vector = geometry_reader.ReadGeometry(model_part)
modeler = NurbsBrepModeler(geometry_reader, model_part)
#modeler.SetUp()