
# Importing the Kratos Library
from KratosMultiphysics import **
from KratosMultiphysics.SolidMechanicsApplication import **
from KratosMultiphysics.PfemSolidMechanicsApplication import **
from beam_sections_python_utility import SetProperties

def AssignMaterial(Properties):
*loop materials
# GUI property identifier: Property1
# GUI material identifier: Steel_AISI1059
*format "%i"
    prop_id = *MatNum;
    prop = Properties[prop_id]
    mat = *MatProp(CONSTITUTIVE_LAW_NAME)();
    prop.SetValue(CONSTITUTIVE_LAW, mat.Clone());
*end materials