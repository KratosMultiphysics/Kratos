import KratosMultiphysics
import KratosMultiphysics.ExternalSolversApplication as KratosSolvers
import KratosMultiphysics.ConvectionDiffusionApplication as KratosConvDiff
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid
import KratosMultiphysics.PoromechanicsApplication as KratosPoro
import KratosMultiphysics.DamApplication as KratosDam

import dam_main as Main

solution = Main.Solution()
solution.Run()
