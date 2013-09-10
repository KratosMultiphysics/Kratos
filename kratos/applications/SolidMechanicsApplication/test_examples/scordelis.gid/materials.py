
# Importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from apply_sections import SetProperties
def AssignMaterial(Properties):
# GUI property identifier: Property1
# GUI material identifier: Steel_AISI1059
    prop_id = 1;
    prop = Properties[prop_id]
    mat = Isotropic2D();
    prop.SetValue(CONSTITUTIVE_LAW, mat.Clone());
