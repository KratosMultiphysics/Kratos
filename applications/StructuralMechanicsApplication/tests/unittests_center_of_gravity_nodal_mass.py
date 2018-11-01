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
        dim = 5
        nr_nodes = 6
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("structural_part")
        mp.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = dim
        self._add_variables(mp)

        ar_var = 2.12
        COG_analy=0.0
        t_mass = 0.0

        #create nodes
        dx = 2.4
        for i in range(nr_nodes):
            mp.CreateNewNode(i+1,i*dx,0.00,0.00) #every next element has mass raised to the power one more than the previous.
            COG_analy += ( ar_var**(i+1) ) * (i*dx)
            t_mass += ar_var**(i+1)
        COG_analy /= t_mass
       # print('cog analy', COG_analy)

        #add dofs
        self._add_dofs(mp)

        #create Element

        elem1 = mp.CreateNewElement("NodalConcentratedElement2D1N", 1, [1], mp.GetProperties()[0])
        elem2 = mp.CreateNewElement("NodalConcentratedElement2D1N", 2, [2], mp.GetProperties()[0])
        elem3 = mp.CreateNewElement("NodalConcentratedElement3D1N", 3, [3], mp.GetProperties()[0])
        elem4 = mp.CreateNewElement("NodalConcentratedElement3D1N", 4, [4], mp.GetProperties()[0])
        elem5 = mp.CreateNewElement("NodalConcentratedElement3D1N", 5, [5], mp.GetProperties()[0])
        elem6 = mp.CreateNewElement("NodalConcentratedElement3D1N", 6, [6], mp.GetProperties()[0])

        

        elem1.SetValue(KratosMultiphysics.NODAL_MASS,ar_var)
        elem2.SetValue(KratosMultiphysics.NODAL_MASS,ar_var**2)
        elem3.SetValue(KratosMultiphysics.NODAL_MASS,ar_var**3)
        elem4.SetValue(KratosMultiphysics.NODAL_MASS,ar_var**4)
        elem5.SetValue(KratosMultiphysics.NODAL_MASS,ar_var**5)
        elem6.SetValue(KratosMultiphysics.NODAL_MASS,ar_var**6)

    

        mass_process = StructuralMechanicsApplication.ComputeCenterOfGravityProcess(mp)
        print("printing cog for nodal mass")
        mass_process.Execute()
        COG_test = mp.ProcessInfo[StructuralMechanicsApplication.CENTER_OF_GRAVITY]
       # print('cog test', COG_test[0])

        self.assertAlmostEqual(COG_test[0], COG_analy, 5)
        self.assertAlmostEqual(COG_test[1], 0, 5)
        self.assertAlmostEqual(COG_test[2], 0, 5)

    

if __name__ == '__main__':
    KratosUnittest.main()


