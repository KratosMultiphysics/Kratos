import KratosMultiphysics
import KratosMultiphysics.FemToDemApplication
import MainFemDem

model = KratosMultiphysics.Model()
MainFemDem.FEM_Solution(model).Run()

