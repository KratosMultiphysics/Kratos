from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest

try:
    import KratosMultiphysics.EigenSolversApplication as EigenSolversApplication
    eigen_solvers_is_available = True
except ImportError:
    eigen_solvers_is_available = False

#A simply supported square plate under compressive loading is computed
#The test compares the buckling load/multiplier between one model with symmetry conditions (quarter of the plate) and a full model
#Futhermore the results are checked against a reference solution from abaqus (element S4R)
#http://130.149.89.49:2080/v6.8/books/bmk/default.htm?startat=ch01s02ach17.html


class BaseTestPrebucklingAnalysis(KratosUnittest.TestCase):
    @classmethod
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
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_MOMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)


    def _create_nodes(self,mp, NumOfNodes,length):
        # Create Nodes
        counter = 0
        for y in range(NumOfNodes):
            for x in range(NumOfNodes):
                counter = counter + 1
                mp.CreateNewNode(counter,x/(NumOfNodes-1)*length,y/(NumOfNodes-1)*length,0.0)


    def _create_elements(self,mp,NumOfNodes):
        element_name = "ShellThinElementCorotational3D4N"
        counter = 0
        for y in range( NumOfNodes-1 ):
            for x in range(NumOfNodes - 1):
                # Aligned counter-clockwise
                counter = counter + 1
                node1 = NumOfNodes*(y)+(x+1)
                node2 = NumOfNodes*(y)+(x+2)
                node3 = NumOfNodes*(y+1) + (x+2)
                node4 = NumOfNodes*(y+1) + (x+1)
                mp.CreateNewElement( element_name, counter, [ node1, node2, node3, node4 ],mp.GetProperties()[0] )


    def _apply_material_properties(self,mp):
        # Define properties
        mp.GetProperties()[0].SetValue(KratosMultiphysics.YOUNG_MODULUS,1e8)
        mp.GetProperties()[0].SetValue(KratosMultiphysics.POISSON_RATIO,0.3)
        mp.GetProperties()[0].SetValue(KratosMultiphysics.THICKNESS,0.01)
        mp.GetProperties()[0].SetValue(KratosMultiphysics.DENSITY,1.0)

        cl = StructuralMechanicsApplication.LinearElasticPlaneStress2DLaw()
        mp.GetProperties()[0].SetValue(KratosMultiphysics.CONSTITUTIVE_LAW,cl)

    def _apply_BCs_sym_vertical(self,mp):
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_X, True, mp.Nodes)
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.ROTATION_Y, True, mp.Nodes)
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.ROTATION_Z, True, mp.Nodes)

    def _apply_BCs_sym_horizontal(self,mp):
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_Y, True, mp.Nodes)
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.ROTATION_X, True, mp.Nodes)
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.ROTATION_Z, True, mp.Nodes)

    def _apply_BCs_simple_vertical(self,mp):
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_Z, True, mp.Nodes)
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.ROTATION_X, True, mp.Nodes)

    def _apply_BCs_simple_horizontal(self,mp):
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_Z, True, mp.Nodes)
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.ROTATION_Y, True, mp.Nodes)

    def _solve_prebuckling_problem(self,mp,NumOfNodes,iterations,echo=0):
        eigensolver_settings = KratosMultiphysics.Parameters("""
        {
            "max_iteration"         : 1000,
            "tolerance"             : 1e-6,
            "number_of_eigenvalues" : 2,
            "echo_level"            : 0,
            "normalize_eigenvectors": true
        }
        """)

        eigen_solver = EigenSolversApplication.EigensystemSolver(eigensolver_settings)
        eigen_solver_ = KratosMultiphysics.ResidualBasedEliminationBuilderAndSolver(eigen_solver)
        convergence_criterion = KratosMultiphysics.DisplacementCriteria(1e-4,1e-9)
        scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()
        builder_and_solver = KratosMultiphysics.ResidualBasedEliminationBuilderAndSolver(linear_solver)
        convergence_criterion.SetEchoLevel(echo)

        eig_strategy = StructuralMechanicsApplication.PrebucklingStrategy( mp,
                                                                           scheme,
                                                                           eigen_solver_,
                                                                           builder_and_solver,
                                                                           convergence_criterion,
                                                                           10,
                                                                           1.0,
                                                                           0.0005,
                                                                           0.5,
                                                                           0.005 )
        eig_strategy.SetEchoLevel(echo)
        LoadFactor = []
        for i in range(iterations):
            eig_strategy.Solve()
            if( i%2 == 1):
                LoadFactor.append( mp.ProcessInfo[StructuralMechanicsApplication.EIGENVALUE_VECTOR][0] )
            self._updateConditions(mp,NumOfNodes)
        return LoadFactor

    def _updateConditions(self,mp,NumOfNodes):
        conditions = mp.GetConditions()
        load_multiplier = mp.ProcessInfo[KratosMultiphysics.TIME]
        for i, condition in enumerate(conditions):
            tmp =  condition.GetValue(StructuralMechanicsApplication.POINT_LOAD)
            load = 0.25
            # Corner nodes have different values
            if i < 2 or i == NumOfNodes or i == NumOfNodes+1:
                load = 0.125
            if tmp[0] > 0.0:
                tmp[0] = load_multiplier*load
            else:
                tmp[0] = -1*load_multiplier*load
            condition.SetValue( StructuralMechanicsApplication.POINT_LOAD, tmp )


    def _set_conditions(self,mp,nodes,NumOfNodes,direction,counter):
        # Corner Nodes
        cond1 = mp.CreateNewCondition("PointLoadCondition3D1N",counter,[nodes[0]],mp.GetProperties()[0])
        cond2 = mp.CreateNewCondition("PointLoadCondition3D1N",counter+1,[nodes[1]],mp.GetProperties()[0])
        counter = counter + 1
        load_on_cond1 = KratosMultiphysics.Vector(3)
        load_on_cond1[0] = direction*0.125
        load_on_cond1[1] = 0.0
        load_on_cond1[2] = 0.0
        cond1.SetValue(StructuralMechanicsApplication.POINT_LOAD, load_on_cond1)
        cond2.SetValue(StructuralMechanicsApplication.POINT_LOAD, load_on_cond1)
        load_on_cond2 = KratosMultiphysics.Vector(3)
        load_on_cond2[0] = direction*0.25
        load_on_cond2[1] = 0.0
        load_on_cond2[2] = 0.0
        # Center Nodes
        if( direction > 0):
            max_ = NumOfNodes*(NumOfNodes-1)
        else:
            max_ = NumOfNodes**2-1

        for i in range(nodes[0]+NumOfNodes,max_,NumOfNodes):
            counter = counter + 1
            cond_tmp = mp.CreateNewCondition("PointLoadCondition3D1N",counter,[i],mp.GetProperties()[0])
            cond_tmp.SetValue(StructuralMechanicsApplication.POINT_LOAD, load_on_cond2)

    def _apply_Bcs_Symmetry(self,mp,NumOfNodes):
        # Create a submodelpart for dirichlet boundary conditions
        bcs_dirichlet_left = mp.CreateSubModelPart("BoundaryCondtionsDirichlet_left")
        bcs_dirichlet_right = mp.CreateSubModelPart("BoundaryCondtionsDirichlet_right")
        bcs_dirichlet_lower = mp.CreateSubModelPart("BoundaryCondtionsDirichlet_lower")
        bcs_dirichlet_upper = mp.CreateSubModelPart("BoundaryCondtionsDirichlet_upper")
        bcs_dirichlet_left.AddNodes(range(1,NumOfNodes*NumOfNodes,NumOfNodes))
        self._apply_BCs_sym_vertical(bcs_dirichlet_left)
        # Right Edge
        bcs_dirichlet_right.AddNodes(range(NumOfNodes, NumOfNodes*NumOfNodes+1, NumOfNodes))
        self._apply_BCs_simple_vertical(bcs_dirichlet_right)
        # Lower Edge
        bcs_dirichlet_lower.AddNodes(range(1,NumOfNodes+1,1))
        self._apply_BCs_sym_horizontal( bcs_dirichlet_lower )
        # Upper Edge
        bcs_dirichlet_upper.AddNodes( range(( NumOfNodes-1)*NumOfNodes+1, NumOfNodes*NumOfNodes+1, 1) )
        self._apply_BCs_simple_horizontal(bcs_dirichlet_upper )

    def _apply_Bcs_Full(self,mp,NumOfNodes):
        # Create a submodelpart for dirichlet boundary conditions
        bcs_dirichlet_left = mp.CreateSubModelPart("BoundaryCondtionsDirichlet_left")
        bcs_dirichlet_right = mp.CreateSubModelPart("BoundaryCondtionsDirichlet_right")
        bcs_dirichlet_lower = mp.CreateSubModelPart("BoundaryCondtionsDirichlet_lower")
        bcs_dirichlet_upper = mp.CreateSubModelPart("BoundaryCondtionsDirichlet_upper")
        bcs_dirichlet_left.AddNodes(range(1,NumOfNodes*NumOfNodes,NumOfNodes))
        self._apply_BCs_simple_vertical(bcs_dirichlet_left)
        # Right Edge
        bcs_dirichlet_right.AddNodes(range(NumOfNodes, NumOfNodes*NumOfNodes+1, NumOfNodes))
        self._apply_BCs_simple_vertical(bcs_dirichlet_right)
        # Lower Edge
        bcs_dirichlet_lower.AddNodes(range(1,NumOfNodes+1,1))
        self._apply_BCs_simple_horizontal( bcs_dirichlet_lower )
        # Upper Edge
        bcs_dirichlet_upper.AddNodes( range(( NumOfNodes-1)*NumOfNodes+1, NumOfNodes**2+1, 1) )
        self._apply_BCs_simple_horizontal(bcs_dirichlet_upper )
        bcs_dirichlet_center_up = mp.CreateSubModelPart("BoundaryCondtionsDirichlet_center_up")
        bcs_dirichlet_center_low = mp.CreateSubModelPart("BoundaryCondtionsDirichlet_center_low")
        bcs_dirichlet_center = mp.CreateSubModelPart("BoundaryCondtionsDirichlet_center")
        # Apply additional boundary conditions at the center nodes to statically determinate the plate
        # Uppder center Node
        bcs_dirichlet_center_up.AddNodes([int(NumOfNodes**2-(NumOfNodes-1)/2)])
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_X, True, bcs_dirichlet_center_up.Nodes)
        # Lower Center Node
        bcs_dirichlet_center_low.AddNodes([int((NumOfNodes+1)/2)])
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_X, True, bcs_dirichlet_center_low.Nodes)
        # Center Node
        bcs_dirichlet_center.AddNodes([int((NumOfNodes**2)/2 + 1)])
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_X, True, bcs_dirichlet_center.Nodes)
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_Y, True, bcs_dirichlet_center.Nodes)


    def _set_up_system(self,current_model,NumOfNodes,length,symmetry):
        if( symmetry ):
            mp = current_model.CreateModelPart("Structure_Symmetry")
        else:
            mp = current_model.CreateModelPart("Structure_Full")
        #mp.SetBufferSize(2)

        self._add_variables(mp)
        self._apply_material_properties(mp)
        self._create_nodes(mp,NumOfNodes,length)
        self._add_dofs(mp)
        self._create_elements(mp,NumOfNodes)

        #LeftEdge
        if( symmetry ):
            self._apply_Bcs_Symmetry(mp,NumOfNodes)
            self._set_conditions(mp,[NumOfNodes,NumOfNodes**2],NumOfNodes,-1,1)
        else:
            self._apply_Bcs_Full(mp,NumOfNodes)
            #Loads on right edge
            self._set_conditions(mp,[NumOfNodes,NumOfNodes**2],NumOfNodes,-1,1)
            #Loads on left edge
            self._set_conditions(mp,[1,NumOfNodes*(NumOfNodes-1)+1],NumOfNodes,1,NumOfNodes+1)

        return mp


    def _check_load_multiplier(self,load_multiplier1, load_multiplier2, reference):
        #Check if value stays the same in the first and last loadstep
        self.assertLess( abs(1-load_multiplier1[0]/load_multiplier1[8]), 1.0e-4)
        #Check if both models give same values
        self.assertAlmostEqual(load_multiplier1[0], load_multiplier2[0], 5)
        #Compare value against reference from abaqus
        self.assertLess( abs(1-load_multiplier1[0]/reference), 1.0e-2)

class TestPrebucklingAnalysis(BaseTestPrebucklingAnalysis):
    @KratosUnittest.skipUnless(eigen_solvers_is_available,"EigenSolversApplication not available")
    def test_dynamic_eigenvalue_analysis(self):
        reference_value = 92.80
        #Construct model with symmetry conditions (quarter of the full plate 1x1)
        NumOfNodesPerSide = 5
        LoadSteps = 18
        Length = 1
        current_model_sym = KratosMultiphysics.Model()
        mp_sym = self._set_up_system(current_model_sym,NumOfNodesPerSide,Length,True)
        loadFactor_sym = self._solve_prebuckling_problem(mp_sym,NumOfNodesPerSide,LoadSteps)

        #Construct full model (whole plate 2x2)
        NumOfNodesPerSide = 9
        LoadSteps = 2
        Length = 2
        current_model_full = KratosMultiphysics.Model()
        mp_full = self._set_up_system(current_model_full,NumOfNodesPerSide,Length,False)
        loadFactor_full = self._solve_prebuckling_problem(mp_full,NumOfNodesPerSide,LoadSteps)
        #print(loadFactor_sym)
        #print(loadFactor_full)
        self._check_load_multiplier(loadFactor_sym,loadFactor_full, reference_value )


if __name__ == '__main__':
    KratosUnittest.main()
