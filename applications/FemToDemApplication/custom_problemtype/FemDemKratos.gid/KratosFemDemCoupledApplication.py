import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication
import KratosMultiphysics.FemToDemApplication
import CouplingFemDem

model = KratosMultiphysics.Model()
CouplingFemDem.FEMDEM_Solution(model).Run()