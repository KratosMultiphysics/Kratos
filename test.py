import KratosMultiphysics as Kratos
import numpy as np

model = Kratos.Model()
model_part = model.CreateModelPart("test")
model_part.AddNodalSolutionStepVariable(Kratos.ACCELERATION)

_ = model_part.CreateNewNode(1, 0, 0, 0)
_ = model_part.CreateNewNode(2, 1, 0, 0)
_ = model_part.CreateNewNode(3, 1, 1, 0)
_ = model_part.CreateNewNode(4, 0, 1, 0)

properties = model_part.CreateNewProperties(1)
_ = model_part.CreateNewElement("Element2D3N", 1, [1,2,3], properties)
_ = model_part.CreateNewElement("Element2D3N", 2, [1,3,4], properties)

for node in model_part.Nodes:
    node.SetSolutionStepValue(Kratos.ACCELERATION, Kratos.Array3([node.Id, node.Id + 1, node.Id + 2]))
    node.SetValue(Kratos.VELOCITY, Kratos.Array3([node.Id, node.Id + 1, node.Id + 2]))

for element in model_part.Elements:
    element.SetValue(Kratos.VELOCITY, Kratos.Array3([element.Id, element.Id + 1, element.Id + 2]))

# reading from nodal historical
nodal_hist_exp = Kratos.Expression.NodalExpression(model_part)
Kratos.Expression.VariableExpressionIO.Read(nodal_hist_exp, Kratos.ACCELERATION, True)

# reading from nodal non-historical
nodal_nonhist_exp = Kratos.Expression.NodalExpression(model_part)
Kratos.Expression.VariableExpressionIO.Read(nodal_nonhist_exp, Kratos.VELOCITY, False)

# reading from elemental values in Kratos
element_exp = Kratos.Expression.ElementExpression(model_part)
Kratos.Expression.VariableExpressionIO.Read(element_exp, Kratos.VELOCITY)

# reading from integration points. Following will give an empty array since
# there aren't any gauss point values defined in Element2D3N.
elem_gauss_exp = Kratos.Expression.ElementExpression(model_part)
Kratos.Expression.IntegrationPointExpressionIO.Read(elem_gauss_exp, Kratos.PRESSURE)

# all the expressions are OpenMP parallelized and MPI parallelized
# you can do following arithmatic operations as well which has O(1) cost (irrespective of
# the number of nodes/conditions and elements in model part)
new_exp = element_exp * (-1.0) + 3

# convert any expression (NodalExpression, ElementExpression, ConditionExpression)
# to numpy array. The shape of the numpy array will be [n, m]. "n" is number
# of elements. "m" is 3 (since new_exp, element_exp was reading Kratos.VELOCITY)
numpy_array = new_exp.Evaluate()
numpy_array *= 2.0

# convert the numpy array back to an expression
numpy_expression = Kratos.Expression.ElementExpression(model_part)
Kratos.Expression.CArrayExpressionIO.Read(numpy_expression, numpy_array)

# writing expression data back to model part elements
Kratos.Expression.VariableExpressionIO.Write(numpy_expression, Kratos.ACCELERATION)

# writing expression data back to model part gauss points
Kratos.Expression.IntegrationPointExpressionIO.Write(elem_gauss_exp, Kratos.DISPLACEMENT_X)

# following part is for the hdf5 input output of the expressions

import KratosMultiphysics.HDF5Application as KratosHDF5
from KratosMultiphysics.HDF5Application.core.file_io import OpenHDF5File

# first you need to create an HDF5 file
hdf5_write_params = Kratos.Parameters("""{
    "file_name": "test.h5",
    "file_access_mode": "truncate"
}""")
with OpenHDF5File(hdf5_write_params, model_part) as h5_file:
    prefix_settings = Kratos.Parameters("""{
        "prefix": "/WHERE_YOU_WANT_TO_PUT_DATA_IN_H5/"
    }""")
    expio = KratosHDF5.ExpressionIO(prefix_settings, h5_file)

    expio.Write("elem_exp", element_exp)
    expio.Write("nodal_hist_exp", nodal_hist_exp)
    expio.Write("nodal_nonhist_exp", nodal_nonhist_exp)
    expio.Write("numpy_expression", numpy_expression)

# now the HDF5 reading part
hdf5_read_params = Kratos.Parameters("""{
    "file_name": "test.h5",
    "file_access_mode": "read_only"
}""")
with OpenHDF5File(hdf5_read_params, model_part) as h5_file:
    prefix_settings = Kratos.Parameters("""{
        "prefix": "/WHERE_YOU_WANT_TO_PUT_DATA_IN_H5/"
    }""")
    expio = KratosHDF5.ExpressionIO(prefix_settings, h5_file)

    new_elem_exp = Kratos.Expression.ElementExpression(model_part)
    new_nodal_hist_exp = Kratos.Expression.NodalExpression(model_part)
    new_nodal_nonhist_exp = Kratos.Expression.NodalExpression(model_part)
    new_numpy_expression = Kratos.Expression.ElementExpression(model_part)

    expio.Read("elem_exp", new_elem_exp)
    expio.Read("nodal_hist_exp", new_nodal_hist_exp)
    expio.Read("nodal_nonhist_exp", new_nodal_nonhist_exp)
    expio.Read("numpy_expression", new_numpy_expression)

print("Errors: ")
print(np.linalg.norm((new_elem_exp - element_exp).Evaluate()))
print(np.linalg.norm((new_nodal_hist_exp - nodal_hist_exp).Evaluate()))
print(np.linalg.norm((new_nodal_nonhist_exp - nodal_nonhist_exp).Evaluate()))
print(np.linalg.norm((new_numpy_expression - numpy_expression).Evaluate()))


