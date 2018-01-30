import KratosMultiphysics
import KratosMultiphysics.ConstitutiveModelsApplication
import MainMaterial

MainMaterial.Solution("shear_parameters.json","neohookean_materials.json").Run()

