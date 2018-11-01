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


   

    def _apply_solid_material_properties(self,mp):
        #define properties
        mp.GetProperties()[1].SetValue(KratosMultiphysics.YOUNG_MODULUS,100e3)
        mp.GetProperties()[1].SetValue(KratosMultiphysics.POISSON_RATIO,0.3)
        mp.GetProperties()[1].SetValue(KratosMultiphysics.THICKNESS,0.23)
        mp.GetProperties()[1].SetValue(KratosMultiphysics.DENSITY,25.7)

        cl = StructuralMechanicsApplication.LinearElasticPlaneStrain2DLaw()

        mp.GetProperties()[1].SetValue(KratosMultiphysics.CONSTITUTIVE_LAW,cl)

    
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



    def test_solid_mass(self):

        #starting with 3 node element testing
        dim = 2
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("structural_part")
        mp.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = dim
        mp.SetBufferSize(2)
        self._add_variables(mp)
        self._apply_solid_material_properties(mp)

         #create nodes for 3 node elements
        mp.CreateNewNode(1,0.5,0.5,0.0)
        mp.CreateNewNode(2,0.7,0.2,0.0)
        mp.CreateNewNode(3,0.9,0.8,0.0)
        mp.CreateNewNode(4,0.3,0.7,0.0)
        mp.CreateNewNode(5,0.6,0.6,0.0)

        self._add_dofs(mp)

        #create Element (3 node)
        mp.CreateNewElement("TotalLagrangianElement2D3N", 1, [1,2,5], mp.GetProperties()[1])
        mp.CreateNewElement("TotalLagrangianElement2D3N", 2, [2,3,5], mp.GetProperties()[1])
        mp.CreateNewElement("TotalLagrangianElement2D3N", 3, [3,4,5], mp.GetProperties()[1])
        mp.CreateNewElement("TotalLagrangianElement2D3N", 4, [4,1,5], mp.GetProperties()[1])

        mass_process = StructuralMechanicsApplication.ComputeCenterOfGravityProcess(mp)
        print("printing cog for solid mass")
        mass_process.Execute()
        COG_test = mp.ProcessInfo[StructuralMechanicsApplication.CENTER_OF_GRAVITY]
      #  print("cog TEST", COG_test)
        self.assertAlmostEqual(0.641663, COG_test[0] , 5)
        self.assertAlmostEqual(0.5729174, COG_test[1] , 5)
        self.assertAlmostEqual(0, COG_test[2], 5)



        #creating another model for 4 node element testing
        dim = 2
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("structural_part")
        mp.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = dim
        mp.SetBufferSize(2)
        self._add_variables(mp)
        self._apply_solid_material_properties(mp)


        #create nodes for 4 node elements
        mp.CreateNewNode(6,0.0,0.0,0.0)
        mp.CreateNewNode(7,4.8,0.0,1.8)
        mp.CreateNewNode(8,5.9,0.0,2.5)
        mp.CreateNewNode(9,0.0,1.2,0.0)
        mp.CreateNewNode(10,4.8,1.2,1.8)
        mp.CreateNewNode(11,5.9,1.2,2.5)
        mp.CreateNewNode(12,0.0,3.7,0.0)
        mp.CreateNewNode(13,4.8,3.7,1.8)
        mp.CreateNewNode(14,5.9,3.7,2.5)
        

        self._add_dofs(mp)

        #create Element
        mp.CreateNewElement("TotalLagrangianElement2D4N", 5, [6,7,10,9], mp.GetProperties()[1])
        mp.CreateNewElement("TotalLagrangianElement2D4N", 6, [7,8,11,10], mp.GetProperties()[1])
        mp.CreateNewElement("TotalLagrangianElement2D4N", 7, [9,10,13,12], mp.GetProperties()[1])
        mp.CreateNewElement("TotalLagrangianElement2D4N", 8, [10,11,14,13], mp.GetProperties()[1])

        mass_process = StructuralMechanicsApplication.ComputeCenterOfGravityProcess(mp)
       # print("printing cog for solid mass")
        mass_process.Execute()
        COG_test_4node = mp.ProcessInfo[StructuralMechanicsApplication.CENTER_OF_GRAVITY]
        #print("cog TEST_4node", COG_test_4node)
        self.assertAlmostEqual(2.99816, COG_test_4node[0], 4)
        self.assertAlmostEqual(1.85, COG_test_4node[1] , 4)
        self.assertAlmostEqual(1.15345, COG_test_4node[2] , 4)

if __name__ == '__main__':
    KratosUnittest.main()


