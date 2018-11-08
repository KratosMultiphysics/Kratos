from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest

from math import sqrt, sin, cos, pi, exp, atan

class TestComputeMomentOfInertia(KratosUnittest.TestCase):
    # muting the output
    #KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING) TODO reenable this

    def _apply_beam_material_properties(self,mp,dim):
        #define properties
        mp.GetProperties()[0].SetValue(KratosMultiphysics.YOUNG_MODULUS,210e9)
        mp.GetProperties()[0].SetValue(KratosMultiphysics.DENSITY,7850)
        mp.GetProperties()[0].SetValue(StructuralMechanicsApplication.CROSS_AREA,0.01)
        mp.GetProperties()[0].SetValue(KratosMultiphysics.POISSON_RATIO,0.30)
        mp.GetProperties()[0].SetValue(StructuralMechanicsApplication.TORSIONAL_INERTIA,0.00001)
        mp.GetProperties()[0].SetValue(StructuralMechanicsApplication.I22,0.00001)
        mp.GetProperties()[0].SetValue(StructuralMechanicsApplication.I33,0.00001)

        cl = StructuralMechanicsApplication.LinearElastic3DLaw()
        mp.GetProperties()[0].SetValue(KratosMultiphysics.CONSTITUTIVE_LAW,cl)

    def _apply_shell_material_properties(self,mp):
        #define properties
        mp.GetProperties()[1].SetValue(KratosMultiphysics.YOUNG_MODULUS,100e3)
        mp.GetProperties()[1].SetValue(KratosMultiphysics.POISSON_RATIO,0.3)
        mp.GetProperties()[1].SetValue(KratosMultiphysics.THICKNESS,1.0)
        mp.GetProperties()[1].SetValue(KratosMultiphysics.DENSITY,1.0)

        cl = StructuralMechanicsApplication.LinearElasticPlaneStress2DLaw()

        mp.GetProperties()[1].SetValue(KratosMultiphysics.CONSTITUTIVE_LAW,cl)

    def _apply_orthotropic_shell_material_properties(self,mp):
        #define properties
        # we specify only the properties we need (others are youngs modulus etc)
        num_plies = 3
        orthotropic_props = KratosMultiphysics.Matrix(num_plies,16)
        for row in range(num_plies):
            for col in range(16):
                orthotropic_props[row,col] = 0.0

        # Orthotropic mechanical moduli
        orthotropic_props[0,0] = 0.005 # lamina thickness
        orthotropic_props[0,2] = 2200  # density
        orthotropic_props[1,0] = 0.01  # lamina thickness
        orthotropic_props[1,2] = 1475  # density
        orthotropic_props[2,0] = 0.015 # lamina thickness
        orthotropic_props[2,2] = 520   # density

        mp.GetProperties()[1].SetValue(StructuralMechanicsApplication.SHELL_ORTHOTROPIC_LAYERS,orthotropic_props)

        cl = StructuralMechanicsApplication.LinearElasticOrthotropic2DLaw()

        mp.GetProperties()[1].SetValue(KratosMultiphysics.CONSTITUTIVE_LAW,cl)

    def _apply_solid_material_properties(self,mp):
        #define properties
        mp.GetProperties()[1].SetValue(KratosMultiphysics.YOUNG_MODULUS,100e3)
        mp.GetProperties()[1].SetValue(KratosMultiphysics.POISSON_RATIO,0.3)
        mp.GetProperties()[1].SetValue(KratosMultiphysics.THICKNESS,1.0)
        mp.GetProperties()[1].SetValue(KratosMultiphysics.DENSITY,1.0)

        cl = StructuralMechanicsApplication.LinearElasticPlaneStrain2DLaw()

        mp.GetProperties()[1].SetValue(KratosMultiphysics.CONSTITUTIVE_LAW,cl)

    def _create_shell_nodes(self,mp):
        mp.CreateNewNode(1, -0.5, - 0.45,  0.1)
        mp.CreateNewNode(2,  0.7,  -0.5,   0.2)
        mp.CreateNewNode(3,  0.55,  0.6,   0.15)
        mp.CreateNewNode(4, -0.48,  0.65,  0.0)
        mp.CreateNewNode(5,  0.02, -0.01, -0.15)

    def _create_shell_elements(self,mp,element_name = "ShellThinElementCorotational3D3N"):
        mp.CreateNewElement(element_name, 1, [1,2,5], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 2, [2,3,5], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 3, [3,4,5], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 4, [4,1,5], mp.GetProperties()[1])

    def test_nodal_moi(self):
        dim = 3
        nr_nodes = 4
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("structural_part_nodal_masses")
        mp.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = dim

        #create nodes
        dx = 1.2
        for i in range(nr_nodes):
            mp.CreateNewNode(i+1,i*dx,0.00,0.00)
        #add dofs

        #create Element
        elem1 = mp.CreateNewElement("NodalConcentratedElement2D1N", 1, [1], mp.GetProperties()[0])
        elem2 = mp.CreateNewElement("NodalConcentratedElement2D1N", 2, [2], mp.GetProperties()[0])
        elem3 = mp.CreateNewElement("NodalConcentratedElement3D1N", 3, [3], mp.GetProperties()[0])
        elem4 = mp.CreateNewElement("NodalConcentratedElement3D1N", 4, [4], mp.GetProperties()[0])

        elem1.SetValue(KratosMultiphysics.NODAL_MASS,21.234)
        elem2.SetValue(KratosMultiphysics.NODAL_MASS,5.234)
        elem3.SetValue(KratosMultiphysics.NODAL_MASS,112.234)
        elem4.SetValue(KratosMultiphysics.NODAL_MASS,78.234)

        moi_process = StructuralMechanicsApplication.ComputeMomentOfInertiaProcess(mp)
        moi_process.Execute()
        moment_of_inertia = mp.ProcessInfo[StructuralMechanicsApplication.MOMENT_OF_INERTIA]

        # self.assertAlmostEqual(2.5688903639, moment_of_inertia)

    def test_beam_moi(self):
        dim = 3
        nr_nodes = 11
        nr_elements = nr_nodes-1
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("structural_part_beams")
        mp.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = dim
        self._apply_beam_material_properties(mp,dim)

        #create nodes
        dx = 1.20 / nr_elements
        for i in range(nr_nodes):
            mp.CreateNewNode(i+1,i*dx,0.00,0.00)
        #add dofs

        #create Element
        for i in range(nr_elements):
            elem = mp.CreateNewElement("CrLinearBeamElement3D2N", i+1, [i+1,i+2], mp.GetProperties()[0])

        moi_process = StructuralMechanicsApplication.ComputeMomentOfInertiaProcess(mp)
        moi_process.Execute()
        moment_of_inertia = mp.ProcessInfo[StructuralMechanicsApplication.MOMENT_OF_INERTIA]

        # self.assertAlmostEqual(0.6, moment_of_inertia)

    def test_shell_moi(self):
        dim = 3
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("structural_part_shells")
        mp.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = dim
        mp.SetBufferSize(2)

        self._apply_shell_material_properties(mp)
        self._create_shell_nodes(mp)
        self._create_shell_elements(mp)

        p1 = KratosMultiphysics.Point(1.0, 1.0, 2.5)
        p2 = KratosMultiphysics.Point(1.0, 1.0, 3.5)

        moi_process = StructuralMechanicsApplication.ComputeMomentOfInertiaProcess(mp, p1, p2)
        moi_process.Execute()
        moment_of_inertia = mp.ProcessInfo[StructuralMechanicsApplication.MOMENT_OF_INERTIA]

        # self.assertAlmostEqual(0.0723057, moment_of_inertia)

    def test_orthotropic_shell_moi(self):
        dim = 3
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("structural_part_orthotropic_shells")
        mp.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = dim
        mp.SetBufferSize(2)

        self._apply_orthotropic_shell_material_properties(mp)
        self._create_shell_nodes(mp)
        self._create_shell_elements(mp)

        moi_process = StructuralMechanicsApplication.ComputeMomentOfInertiaProcess(mp)
        moi_process.Execute()
        moment_of_inertia = mp.ProcessInfo[StructuralMechanicsApplication.MOMENT_OF_INERTIA]

        # self.assertAlmostEqual(0.0723057, moment_of_inertia)

    def test_solid_moi(self):
        dim = 2
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("structural_part_solids")
        mp.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = dim
        mp.SetBufferSize(2)
        self._apply_solid_material_properties(mp)

        #create nodes
        mp.CreateNewNode(1,0.5,0.5,0.0)
        mp.CreateNewNode(2,0.7,0.2,0.0)
        mp.CreateNewNode(3,0.9,0.8,0.0)
        mp.CreateNewNode(4,0.3,0.7,0.0)
        mp.CreateNewNode(5,0.6,0.6,0.0)

        #create Element
        mp.CreateNewElement("TotalLagrangianElement2D3N", 1, [1,2,5], mp.GetProperties()[1])
        mp.CreateNewElement("TotalLagrangianElement2D3N", 2, [2,3,5], mp.GetProperties()[1])
        mp.CreateNewElement("TotalLagrangianElement2D3N", 3, [3,4,5], mp.GetProperties()[1])
        mp.CreateNewElement("TotalLagrangianElement2D3N", 4, [4,1,5], mp.GetProperties()[1])

        moi_process = StructuralMechanicsApplication.ComputeMomentOfInertiaProcess(mp)
        moi_process.Execute()
        moment_of_inertia = mp.ProcessInfo[StructuralMechanicsApplication.MOMENT_OF_INERTIA]

        # self.assertAlmostEqual(0.6416666667, moment_of_inertia)

if __name__ == '__main__':
    KratosUnittest.main()


