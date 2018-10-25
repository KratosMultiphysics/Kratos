from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest

from math import sqrt, sin, cos, pi, exp, atan

class TestMassCalculation(KratosUnittest.TestCase):
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

    def _add_dofs(self,mp):
        # Adding dofs AND their corresponding reactions
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_X, KratosMultiphysics.REACTION_MOMENT_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Y, KratosMultiphysics.REACTION_MOMENT_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Z, KratosMultiphysics.REACTION_MOMENT_Z,mp)

    def _add_variables(self,mp):
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_MOMENT)

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

    def _set_and_fill_buffer(self,mp,buffer_size,delta_time):
        # Set buffer size
        mp.SetBufferSize(buffer_size)

        # Fill buffer
        time = mp.ProcessInfo[KratosMultiphysics.TIME]
        time = time - delta_time * (buffer_size)
        mp.ProcessInfo.SetValue(KratosMultiphysics.TIME, time)
        for size in range(0, buffer_size):
            step = size - (buffer_size -1)
            mp.ProcessInfo.SetValue(KratosMultiphysics.STEP, step)
            time = time + delta_time
            #delta_time is computed from previous time in process_info
            mp.CloneTimeStep(time)

        mp.ProcessInfo[KratosMultiphysics.IS_RESTARTED] = False

    def test_nodal_mass(self):
        dim = 3
        nr_nodes = 4
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("structural_part")
        mp.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = dim
        self._add_variables(mp)

        #create nodes
        dx = 1.2
        for i in range(nr_nodes):
            mp.CreateNewNode(i+1,i*dx,0.00,0.00)
        #add dofs
        self._add_dofs(mp)

        #create Element
        elem1 = mp.CreateNewElement("NodalConcentratedElement2D1N", 1, [1], mp.GetProperties()[0])
        elem2 = mp.CreateNewElement("NodalConcentratedElement2D1N", 2, [2], mp.GetProperties()[0])
        elem3 = mp.CreateNewElement("NodalConcentratedElement3D1N", 3, [3], mp.GetProperties()[0])
        elem4 = mp.CreateNewElement("NodalConcentratedElement3D1N", 4, [4], mp.GetProperties()[0])

        expected_elemental_masses_by_id = {
            1 : 21.234,
            2 : 5.234,
            3 : 112.234,
            4 : 78.234
        }
        elem1.SetValue(KratosMultiphysics.NODAL_MASS, expected_elemental_masses_by_id[1])
        elem2.SetValue(KratosMultiphysics.NODAL_MASS, expected_elemental_masses_by_id[2])
        elem3.SetValue(KratosMultiphysics.NODAL_MASS, expected_elemental_masses_by_id[3])
        elem4.SetValue(KratosMultiphysics.NODAL_MASS, expected_elemental_masses_by_id[4])

        mass_process = StructuralMechanicsApplication.TotalStructuralMassProcess(mp)
        mass_process.Execute()

        self.__CheckElementalMasses(mp, expected_elemental_masses_by_id)

        total_mass = mp.ProcessInfo[KratosMultiphysics.NODAL_MASS]
        self.assertAlmostEqual(sum(expected_elemental_masses_by_id.values()), total_mass)

    def test_beam_mass(self):
        dim = 3
        nr_nodes = 11
        nr_elements = nr_nodes-1
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("structural_part")
        mp.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = dim
        self._add_variables(mp)
        self._apply_beam_material_properties(mp,dim)

        #create nodes
        dx = 1.20 / nr_elements
        for i in range(nr_nodes):
            mp.CreateNewNode(i+1,i*dx,0.00,0.00)
        #add dofs
        self._add_dofs(mp)

        #create Element
        for i in range(nr_elements):
            elem = mp.CreateNewElement("CrLinearBeamElement3D2N", i+1, [i+1,i+2], mp.GetProperties()[0])

        mass_process = StructuralMechanicsApplication.TotalStructuralMassProcess(mp)
        mass_process.Execute()

        expected_elemental_masses_by_id = {
            1  : 9.42,
            2  : 9.42,
            3  : 9.42,
            4  : 9.42,
            5  : 9.42,
            6  : 9.42,
            7  : 9.42,
            8  : 9.42,
            9  : 9.42,
            10 : 9.42
        }
        self.__CheckElementalMasses(mp, expected_elemental_masses_by_id)

        total_mass = mp.ProcessInfo[KratosMultiphysics.NODAL_MASS]
        self.assertAlmostEqual(sum(expected_elemental_masses_by_id.values()), total_mass)

    def test_shell_mass(self):
        dim = 3
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("structural_part")
        mp.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = dim
        mp.SetBufferSize(2)

        self._add_variables(mp)
        self._apply_shell_material_properties(mp)
        self._create_shell_nodes(mp)
        self._add_dofs(mp)
        self._create_shell_elements(mp)

        mass_process = StructuralMechanicsApplication.TotalStructuralMassProcess(mp)
        mass_process.Execute()

        expected_elemental_masses_by_id = {
            1 : 0.32856211361019694,
            2 : 0.38250906851995015,
            3 : 0.3509898645260286,
            4 : 0.30527210812650407
        }
        self.__CheckElementalMasses(mp, expected_elemental_masses_by_id)

        total_mass = mp.ProcessInfo[KratosMultiphysics.NODAL_MASS]
        self.assertAlmostEqual(sum(expected_elemental_masses_by_id.values()), total_mass)

    def test_orthotropic_shell_mass(self):
        dim = 3
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("structural_part")
        mp.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = dim
        mp.SetBufferSize(2)

        self._add_variables(mp)
        self._apply_orthotropic_shell_material_properties(mp)
        self._create_shell_nodes(mp)
        self._add_dofs(mp)
        self._create_shell_elements(mp)

        mass_process = StructuralMechanicsApplication.TotalStructuralMassProcess(mp)
        mass_process.Execute()

        expected_elemental_masses_by_id = {
            1 : 11.023258911622108,
            2 : 12.833179248844328,
            3 : 11.77570995484826,
            4 : 10.241879227644212
        }
        self.__CheckElementalMasses(mp, expected_elemental_masses_by_id)

        total_mass = mp.ProcessInfo[KratosMultiphysics.NODAL_MASS]
        self.assertAlmostEqual(sum(expected_elemental_masses_by_id.values()), total_mass)

    def test_solid_mass(self):
        dim = 2
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("structural_part")
        mp.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = dim
        mp.SetBufferSize(2)
        self._add_variables(mp)
        self._apply_solid_material_properties(mp)

        #create nodes
        mp.CreateNewNode(1,0.5,0.5,0.0)
        mp.CreateNewNode(2,0.7,0.2,0.0)
        mp.CreateNewNode(3,0.9,0.8,0.0)
        mp.CreateNewNode(4,0.3,0.7,0.0)
        mp.CreateNewNode(5,0.6,0.6,0.0)

        self._add_dofs(mp)

        #create Element
        mp.CreateNewElement("TotalLagrangianElement2D3N", 1, [1,2,5], mp.GetProperties()[1])
        mp.CreateNewElement("TotalLagrangianElement2D3N", 2, [2,3,5], mp.GetProperties()[1])
        mp.CreateNewElement("TotalLagrangianElement2D3N", 3, [3,4,5], mp.GetProperties()[1])
        mp.CreateNewElement("TotalLagrangianElement2D3N", 4, [4,1,5], mp.GetProperties()[1])

        mass_process = StructuralMechanicsApplication.TotalStructuralMassProcess(mp)
        mass_process.Execute()

        expected_elemental_masses_by_id = {
            1 : 0.025,
            2 : 0.07,
            3 : 0.045,
            4 : 0.02
        }
        self.__CheckElementalMasses(mp, expected_elemental_masses_by_id)

        total_mass = mp.ProcessInfo[KratosMultiphysics.NODAL_MASS]
        self.assertAlmostEqual(sum(expected_elemental_masses_by_id.values()), total_mass)

    def __CheckElementalMasses(self, model_part, expected_elemental_masses_by_id):
        domain_size = model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        for elem in model_part.Elements:
            calculated_elemental_mass = StructuralMechanicsApplication.TotalStructuralMassProcess.CalculateElementMass(elem, domain_size)
            self.assertAlmostEqual(expected_elemental_masses_by_id[elem.Id],
                                   calculated_elemental_mass)



if __name__ == '__main__':
    KratosUnittest.main()


