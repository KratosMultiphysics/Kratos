import KratosMultiphysics
# import KratosMultiphysics.SolidMechanicsApplication
import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.FemToDemApplication
import MainFemDem

model = KratosMultiphysics.Model()
MainFemDem.FEM_Solution(model).Run()