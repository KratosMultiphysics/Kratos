from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.IgaApplication
import matplotlib.pyplot as plt 


iga_model = KratosMultiphysics.Model()
iga_model_part = iga_model.CreateModelPart("IgaModelPart")



with open("geometry.json",'r') as geometry_file:
    iga_geometry_parameters = KratosMultiphysics.Parameters( geometry_file.read())

iga_geometry_reader = KratosMultiphysics.IgaApplication.BrepJsonIO()

skin_model = KratosMultiphysics.Model()
skin_model_part = skin_model.CreateModelPart("SkinModelPart")


embedded_iga_modeler = KratosMultiphysics.IgaApplication.EmbeddedIgaModeler(iga_model_part)
embedded_iga_modeler.ImportGeometry(iga_geometry_reader, iga_geometry_parameters)
embedded_iga_modeler.CreateElements2D(skin_model_part)

coords = embedded_iga_modeler.PrintCurveTessellationPoints()
plt.title("Quad Tessellation")

for i in range(len(coords)):
    plt.plot(coords[i][0], coords[i][1], 'b+')
    # plt.annotate(str(i), xy=(coords[i][0], coords[i][1]))
plt.show()







