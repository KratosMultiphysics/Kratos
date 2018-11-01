from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest

from math import sqrt, sin, cos, pi, exp, atan

class TestComputeCenterOfGravity(KratosUnittest.TestCase):
    #KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

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



    def test_beam_mass(self):
        dim = 3
        nr_nodes = 50
        nr_elements = nr_nodes-1
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("structural_part")
        mp.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = dim
        self._add_variables(mp)
        self._apply_beam_material_properties(mp,dim)

        #create nodes
        dx = 6.8 / nr_elements
        for i in range(nr_nodes):
            mp.CreateNewNode(i+1,i*dx,0.00,0.00)
        #add dofs
        self._add_dofs(mp)

        COG_analytical = 0.0

        #create Element
        for i in range(nr_elements):
            elem = mp.CreateNewElement("CrLinearBeamElement3D2N", i+1, [i+1,i+2], mp.GetProperties()[0])
            COG_analytical += (7850 * 0.01 * dx) * ( (dx/2) +i*dx )
            # density = 7850 and cross sectional area = 0.01
        COG_analytical /= (7850 * 0.01 * dx) * (nr_elements)
        #print('cog analy', COG_analytical)

        mass_process = StructuralMechanicsApplication.ComputeCenterOfGravityProcess(mp)
        print("printing cog for beam mass")
        mass_process.Execute()
        COG_test = mp.ProcessInfo[StructuralMechanicsApplication.CENTER_OF_GRAVITY]
        #print('cog test', COG_test[0])
        self.assertAlmostEqual(COG_test[0], COG_analytical,5)
        self.assertAlmostEqual(COG_test[1], 0, 5)
        self.assertAlmostEqual(COG_test[2], 0, 5)




if __name__ == '__main__':
    KratosUnittest.main()


