from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from beam_sections_python_utility import SetProperties


def AssignMaterial(Properties):
# GUI property identifier: Property1
# GUI material identifier: Steel_AISI1059
    prop_id = 1
    prop = Properties[prop_id]
    mat = HyperElasticUPPlaneStrain2DLaw()
    prop.SetValue(CONSTITUTIVE_LAW, mat.Clone())
