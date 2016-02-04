
# Importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
#from beam_sections_python_utility import SetProperties
def AssignMaterial(Properties):
# GUI property identifier: Property1
# GUI material identifier: Concrete
    prop_id = 1;
    prop = Properties[prop_id]
    mat = IsotropicDamageSimoJuPlaneStress2DLaw();
    prop.SetValue(CONSTITUTIVE_LAW, mat.Clone());
