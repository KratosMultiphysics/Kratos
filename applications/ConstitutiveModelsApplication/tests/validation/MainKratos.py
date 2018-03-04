import KratosMultiphysics
import KratosMultiphysics.ConstitutiveModelsApplication
import MainMaterial

#MainMaterial.Solution("shear_traction_parameters.json","isochoric_ogden_materials.json").Run()
#MainMaterial.Solution("shear_traction_parameters.json","ogden_materials.json").Run()
#MainMaterial.Solution("shear_traction_parameters.json","neohookean_materials.json").Run()

#MainMaterial.Solution("shear_parameters.json","isochoric_ogden_materials.json").Run()
#MainMaterial.Solution("shear_parameters.json","ogden_materials.json").Run()
#MainMaterial.Solution("shear_parameters.json","neohookean_materials.json").Run()

#MainMaterial.Solution("traction_parameters.json","isochoric_ogden_materials.json").Run()
#MainMaterial.Solution("traction_parameters.json","ogden_materials.json").Run()
MainMaterial.Solution("traction_parameters.json","neohookean_materials.json").Run()
