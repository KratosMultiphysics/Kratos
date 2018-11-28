from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest

from math import sqrt, sin, cos, pi, exp, atan

class TestComputeMassMomentOfInertia(KratosUnittest.TestCase):
    # muting the output
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

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
        mp.CreateNewNode(1,0.0,0.0,0.0)
        mp.CreateNewNode(2,0.2,0.0,0.0)
        mp.CreateNewNode(3,0.2,0.1,0.0)
        mp.CreateNewNode(4,0.0,0.1,0.0)
        mp.CreateNewNode(5,0.6,0.0,0.0)
        mp.CreateNewNode(6,0.6,0.1,0.0)
        mp.CreateNewNode(7,1.0,0.0,0.0)
        mp.CreateNewNode(8,1.0,0.1,0.0)
        mp.CreateNewNode(9,0.0,0.4,0.0)
        mp.CreateNewNode(10,0.2,0.4,0.0)
        mp.CreateNewNode(11,0.6,0.4,0.0)
        mp.CreateNewNode(12,1.0,0.4,0.0)
        mp.CreateNewNode(13,1.0,0.5,0.0)
        mp.CreateNewNode(14,0.6,0.5,0.0)
        mp.CreateNewNode(15,0.2,0.5,0.0)
        mp.CreateNewNode(16,0.0,0.5,0.0)

    def _create_shell_elements(self,mp,element_name = "ShellThinElementCorotational3D4N"):
        mp.CreateNewElement(element_name, 1, [1,2,3,4], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 2, [2,5,6,3], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 3, [5,7,8,6], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 4, [4,3,10,9], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 5, [3,6,11,10], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 6, [6,8,12,11], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 7, [9,10,15,16], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 8, [10,11,14,15], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 9, [11,12,13,14], mp.GetProperties()[1])

    def test_nodal_moi(self):
        dim = 3
        nr_nodes = 4
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("structural_part_nodal_masses")
        mp.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = dim

        # create nodes
        dx = 1.2
        for i in range(nr_nodes):
            mp.CreateNewNode(i+1,i*dx,0.00,0.00)

        # create elements
        elem1 = mp.CreateNewElement("NodalConcentratedElement2D1N", 1, [1], mp.GetProperties()[0])
        elem2 = mp.CreateNewElement("NodalConcentratedElement2D1N", 2, [2], mp.GetProperties()[0])
        elem3 = mp.CreateNewElement("NodalConcentratedElement3D1N", 3, [3], mp.GetProperties()[0])
        elem4 = mp.CreateNewElement("NodalConcentratedElement3D1N", 4, [4], mp.GetProperties()[0])

        elem1.SetValue(KratosMultiphysics.NODAL_MASS,21.234)
        elem2.SetValue(KratosMultiphysics.NODAL_MASS,5.234)
        elem3.SetValue(KratosMultiphysics.NODAL_MASS,112.234)
        elem4.SetValue(KratosMultiphysics.NODAL_MASS,78.234)

        p1 = KratosMultiphysics.Point(1.8, 0.0, 0.0)
        p2 = KratosMultiphysics.Point(1.8, 2.0, 0.0)

        moi_process = StructuralMechanicsApplication.ComputeMassMomentOfInertiaProcess(mp, p1, p2)
        moi_process.Execute()
        moment_of_inertia = mp.ProcessInfo[StructuralMechanicsApplication.MASS_MOMENT_OF_INERTIA]
        self.assertAlmostEqual(364.5648, moment_of_inertia)

    def test_beam_moi(self):
        dim = 3
        nr_nodes = 100
        nr_elements = nr_nodes-1
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("structural_part_beams")
        mp.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = dim
        self._apply_beam_material_properties(mp,dim)

        #create nodes
        dx = 1.20 / nr_elements
        for i in range(nr_nodes):
            mp.CreateNewNode(i+1,i*dx,0.00,0.00)

        # create elements
        for i in range(nr_elements):
            elem = mp.CreateNewElement("CrLinearBeamElement3D2N", i+1, [i+1,i+2], mp.GetProperties()[0])
        p1 = KratosMultiphysics.Point(0.6, 0.0, 0.0)
        p2 = KratosMultiphysics.Point(0.6, 2.0, 0.0)
        moi_process = StructuralMechanicsApplication.ComputeMassMomentOfInertiaProcess(mp, p1, p2)
        moi_process.Execute()
        moment_of_inertia = mp.ProcessInfo[StructuralMechanicsApplication.MASS_MOMENT_OF_INERTIA]
        self.assertAlmostEqual(11.3028466, moment_of_inertia)

    def test_shell_moi(self):
        dim = 3
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("structural_part_shells")
        mp.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = dim
        mp.SetBufferSize(2)

        self._apply_shell_material_properties(mp)
        self._create_shell_nodes(mp)
        self._create_shell_elements(mp)

        p1 = KratosMultiphysics.Point(0.5, 0.25, 0.0)
        p2 = KratosMultiphysics.Point(0.5, 0.25, 1.0)

        moi_process = StructuralMechanicsApplication.ComputeMassMomentOfInertiaProcess(mp, p1, p2)
        moi_process.Execute()
        moment_of_inertia = mp.ProcessInfo[StructuralMechanicsApplication.MASS_MOMENT_OF_INERTIA]

        self.assertAlmostEqual(0.044, moment_of_inertia)

    def test_orthotropic_shell_moi(self):
        dim = 3
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("structural_part_orthotropic_shells")
        mp.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = dim
        mp.SetBufferSize(2)

        self._apply_orthotropic_shell_material_properties(mp)
        self._create_shell_nodes(mp)
        self._create_shell_elements(mp)

        p1 = KratosMultiphysics.Point(0.5, 0.25, 0.0)
        p2 = KratosMultiphysics.Point(0.5, 0.25, 1.0)
        moi_process = StructuralMechanicsApplication.ComputeMassMomentOfInertiaProcess(mp, p1, p2)
        moi_process.Execute()
        moment_of_inertia = mp.ProcessInfo[StructuralMechanicsApplication.MASS_MOMENT_OF_INERTIA]

        self.assertAlmostEqual(1.4762, moment_of_inertia)

    def test_solid_moi(self):
        dim = 2
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("structural_part_solids")
        mp.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = dim
        mp.SetBufferSize(2)
        self._apply_solid_material_properties(mp)

        # create nodes
        mp.CreateNewNode(1,0.0,0.0,0.0)
        mp.CreateNewNode(2,0.3,0.0,0.0)
        mp.CreateNewNode(3,0.3,0.2,0.0)
        mp.CreateNewNode(4,0.0,0.2,0.0)
        mp.CreateNewNode(5,0.8,0.0,0.0)
        mp.CreateNewNode(6,0.8,0.2,0.0)
        mp.CreateNewNode(7,0.8,0.5,0.0)
        mp.CreateNewNode(8,0.3,0.5,0.0)
        mp.CreateNewNode(9,0.0,0.5,0.0)

        # create elements
        mp.CreateNewElement("TotalLagrangianElement2D4N", 1, [1,2,3,4], mp.GetProperties()[1])
        mp.CreateNewElement("TotalLagrangianElement2D4N", 2, [2,5,6,3], mp.GetProperties()[1])
        mp.CreateNewElement("TotalLagrangianElement2D4N", 3, [3,6,7,8], mp.GetProperties()[1])
        mp.CreateNewElement("TotalLagrangianElement2D4N", 4, [4,3,8,9], mp.GetProperties()[1])
        p1 = KratosMultiphysics.Point(0.4, 0.25, 0.0)
        p2 = KratosMultiphysics.Point(0.4, 0.25, 1.0)
        moi_process = StructuralMechanicsApplication.ComputeMassMomentOfInertiaProcess(mp,p1,p2)
        moi_process.Execute()
        moment_of_inertia = mp.ProcessInfo[StructuralMechanicsApplication.MASS_MOMENT_OF_INERTIA]

        self.assertAlmostEqual(0.021, moment_of_inertia)

if __name__ == '__main__':
    KratosUnittest.main()


