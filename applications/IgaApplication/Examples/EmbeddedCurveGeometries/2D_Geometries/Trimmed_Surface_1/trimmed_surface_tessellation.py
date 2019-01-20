from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.IgaApplication
import matplotlib.pyplot as plt 
import numpy as np


iga_model = KratosMultiphysics.Model()
iga_model_part = iga_model.CreateModelPart("IgaModelPart")



with open("geometry.json",'r') as geometry_file:
    iga_geometry_parameters = KratosMultiphysics.Parameters( geometry_file.read())

iga_geometry_reader = KratosMultiphysics.IgaApplication.BrepJsonIO()

embedded_iga_modeler = KratosMultiphysics.IgaApplication.EmbeddedIgaModeler(iga_model_part)
embedded_iga_modeler.ImportGeometry(iga_geometry_reader, iga_geometry_parameters)

##################################################################
## Print tessellated Points and Triangles
##################################################################

tes_coords = embedded_iga_modeler.PrintParameterCurveTessellationPoints()

plt.title("Tessellation and Triangulation of a trimmed surface")

for i in range(len(tes_coords)): 
    plt.plot(tes_coords[i][0],tes_coords[i][1], 'b+')
    plt.annotate(str(i),xy=(tes_coords[i][0],tes_coords[i][1]))


tri_coords = embedded_iga_modeler.PrintTriangulationPoints()

for i in range (0,len(tri_coords), 3):
    plt.plot([tri_coords[i][0], tri_coords[i+1][0], tri_coords[i+2][0], tri_coords[i][0]], [tri_coords[i][1], tri_coords[i+1][1], tri_coords[i+2][1], tri_coords[i][1]])


##################################################################
## Delaunay Triangulation of tessellated patch
##################################################################

# points1 = []
# for i in range(len(tesX)):
#     points1.append([tesX[i], tesY[i]])


# points = np.array(points1)

# from scipy.spatial import Delaunay
# tri = Delaunay(points)

# plt.triplot(points[:,0], points[:,1], tri.simplices.copy())
# plt.plot(points[:,0], points[:,1], 'o')







plt.show()








