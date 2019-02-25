import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication
import KratosMultiphysics.FemToDemApplication
import CouplingFemDem3D

model = KratosMultiphysics.Model()
CouplingFemDem3D.FEMDEM3D_Solution(model).Run()