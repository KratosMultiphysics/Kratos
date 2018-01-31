import KratosMultiphysics
import KratosMultiphysics.ConstitutiveModelsApplication
import MainMaterial

MainMaterial.Solution("shear_traction_parameters.json","ogden_materials.json").Run()

