
# Importing the Kratos Library
from KratosMultiphysics import **
from KratosMultiphysics.SolidMechanicsApplication import **
from beam_sections_python_utility import SetProperties

*loop materials
def AssignMaterial(Properties):
# GUI property identifier: Property1
# GUI material identifier: Steel_AISI1059
*format "%i"
    prop_id = *MatNum;
    prop = Properties[prop_id]
    mat = *MatProp(CONSTITUTIVE_LAW_NAME)();
    prop.SetValue(CONSTITUTIVE_LAW, mat.Clone());
*end materials