import KratosMultiphysics


import KratosMultiphysics.CableNetApplication as CableNetApplication
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest


import sys

class PlanarMovementRestriction(KratosUnittest.TestCase):
    def setUp(self):
        pass

    def _add_dofs(self,mp):
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,mp)
    def _add_explicit_variables(self,mp):
        mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.MIDDLE_VELOCITY)
        mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.FRACTIONAL_ACCELERATION)
        mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.FRACTIONAL_ANGULAR_ACCELERATION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_MASS)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.FORCE_RESIDUAL)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.RESIDUAL_VECTOR)
        mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.MIDDLE_ANGULAR_VELOCITY)
        mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.NODAL_INERTIA)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.MOMENT_RESIDUAL)
    def _add_variables(self,mp):
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_LOAD)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        self._add_explicit_variables(mp)

    def _add_constitutive_law(self,mp):
        mp.GetProperties()[0].SetValue(KratosMultiphysics.CONSTITUTIVE_LAW,StructuralMechanicsApplication.TrussConstitutiveLaw())

    def _apply_material_properties(self,mp,dim,is_process_test):
        #define properties
        mp.GetProperties()[0].SetValue(KratosMultiphysics.YOUNG_MODULUS,1e8)
        mp.GetProperties()[0].SetValue(KratosMultiphysics.DENSITY,7850)
        mp.GetProperties()[0].SetValue(StructuralMechanicsApplication.CROSS_AREA,1e-2)

        mp.GetProperties()[1].SetValue(KratosMultiphysics.YOUNG_MODULUS,1e6)
        mp.GetProperties()[1].SetValue(KratosMultiphysics.NORMAL_TO_WALL,[0.0, 1.0, 0.0])
        mp.GetProperties()[1].SetValue(CableNetApplication.WALL_POSITION,[0.0, 2.0, 0.0])
    
    
    def _apply_BCs(self,mp):
            KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_Z, True, mp.Nodes)


    def _create_dynamic_explicit_strategy(self,mp):
        scheme = StructuralMechanicsApplication.ExplicitCentralDifferencesScheme(0.00,0.00,0.00)
        strategy = StructuralMechanicsApplication.MechanicalExplicitStrategy(mp,scheme,0,0,1)
        strategy.SetEchoLevel(0)
        strategy.Check()
        return strategy


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


    def _set_up_standard_test(self,is_process_test):
        dim = 3
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        mp.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, dim)
        self._add_variables(mp)
        self._apply_material_properties(mp,dim,is_process_test)
        self._add_constitutive_law(mp)
        #create nodes
        n1 = mp.CreateNewNode(1,0.0,2.2,0.0)
        n2 = mp.CreateNewNode(2,1.0,2.2,0.0)
        n3 = mp.CreateNewNode(3,1.0,3.0,0.0)
        n4 = mp.CreateNewNode(4,0.0,3.0,0.0)

        #create elements
        mp.CreateNewElement("TrussElement3D2N", 1, [1, 2], mp.GetProperties()[0])
        mp.CreateNewElement("TrussElement3D2N", 2, [2, 3], mp.GetProperties()[0])
        mp.CreateNewElement("TrussElement3D2N", 3, [3, 4], mp.GetProperties()[0])
        mp.CreateNewElement("TrussElement3D2N", 4, [4, 1], mp.GetProperties()[0])
        mp.CreateNewElement("TrussElement3D2N", 5, [1, 3], mp.GetProperties()[0])
        mp.CreateNewElement("TrussElement3D2N", 6, [2, 4], mp.GetProperties()[0])

        gravity = -9.81
        n1.SetSolutionStepValue(KratosMultiphysics.VOLUME_ACCELERATION_Y,0,gravity)
        n2.SetSolutionStepValue(KratosMultiphysics.VOLUME_ACCELERATION_Y,0,gravity)
        n3.SetSolutionStepValue(KratosMultiphysics.VOLUME_ACCELERATION_Y,0,gravity)
        n4.SetSolutionStepValue(KratosMultiphysics.VOLUME_ACCELERATION_Y,0,gravity)

        #add dofs
        self._add_dofs(mp)
        #create condition
        mp.CreateNewCondition("PlanarMovementRestrictionCondition3D1N", 1, [1], mp.GetProperties()[1])
        mp.CreateNewCondition("PlanarMovementRestrictionCondition3D1N", 2, [2], mp.GetProperties()[1])
        mp.CreateNewCondition("PlanarMovementRestrictionCondition3D1N", 3, [3], mp.GetProperties()[1])
        mp.CreateNewCondition("PlanarMovementRestrictionCondition3D1N", 4, [4], mp.GetProperties()[1])

        bcs_z = mp.CreateSubModelPart("Dirichlet_Z")
        bcs_z.AddNodes([1,2,3,4])

        #apply boundary conditions
        self._apply_BCs(bcs_z)

        return mp, bcs_z

    def test_explicit_dynamic(self):

        mp, bcs_z = self._set_up_standard_test(is_process_test=False)


        #loop over time
        time_end = 4.0
        time_delta = 0.01
        time_i = 0.0
        self._set_and_fill_buffer(mp,2,time_delta)

        

        strategy = self._create_dynamic_explicit_strategy(mp)
        while (time_i < time_end):
            time_i += time_delta
            mp.CloneTimeStep(time_i)

            #solve + compare
            strategy.Solve()

            disp_solution = min([bcs_z.Nodes[i+1].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y) for i in range(4)])
            self.assertGreater(disp_solution, -0.26) # allow a tollerance in the results

if __name__ == '__main__':
    KratosUnittest.main()
