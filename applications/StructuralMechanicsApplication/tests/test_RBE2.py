# Kratos Imports
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics as KM
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import pathlib
import sys

from KratosMultiphysics.StructuralMechanicsApplication.RBE2_process import ApplyRbe2Process

#from pathlib import Path
#python_scripts_path = (Path(__file__).resolve().parents[1] / "python_scripts")
#sys.path.insert(0, str(python_scripts_path))
#import RBE2_process


def GetFullPathToFile(fileName):
    return pathlib.Path(__file__).absolute().parent / fileName

class TestRBE2(KratosUnittest.TestCase):
    def setUp(self):
        pass

    def _add_variables(self,mp):
        mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KM.ROTATION)
        mp.AddNodalSolutionStepVariable(KM.REACTION)
        mp.AddNodalSolutionStepVariable(KM.REACTION_MOMENT)


    def _add_dofs(self,mp):
        # Adding the dofs AND their corresponding reaction!
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_X, KM.REACTION_X,mp)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y,mp)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z,mp)

        KM.VariableUtils().AddDof(KM.ROTATION_X, KM.REACTION_MOMENT_X,mp)
        KM.VariableUtils().AddDof(KM.ROTATION_Y, KM.REACTION_MOMENT_Y,mp)
        KM.VariableUtils().AddDof(KM.ROTATION_Z, KM.REACTION_MOMENT_Z,mp)

    def _create_nodes(self,mp):
        mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2,  0.5, 0.0, 0.0)
        mp.CreateNewNode(3,  1.0, 0.0, 0.0)
        mp.CreateNewNode(4,  0.0, 0.0, -2.0)
        mp.CreateNewNode(5,  0.5, 0.0, -2.0)
        mp.CreateNewNode(6,  1.0, 0.0, -2.0)
        mp.CreateNewNode(7,  0.0, 0.0, -4.0)
        mp.CreateNewNode(8,  0.5, 0.0, -4.0)
        mp.CreateNewNode(9,  1.0, 0.0, -4.0)
        mp.CreateNewNode(10,  0.0, 0.0, -6.0)
        mp.CreateNewNode(11,  0.5, 0.0, -6.0)
        mp.CreateNewNode(12,  1.0, 0.0, -6.0)
        mp.CreateNewNode(13,  2.0, 0.0, -12.0)

    def _create_elements(self,mp,element_name):
        mp.CreateNewElement(element_name, 1, [1,2,5,4], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 2, [2,3,6,5], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 3, [4,5,8,7], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 4, [5,6,9,8], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 5, [7,8,11,10], mp.GetProperties()[1])
        mp.CreateNewElement(element_name, 6, [8,9,12,11], mp.GetProperties()[1])


    def _apply_dirichlet_BCs(self,mp):
        KM.VariableUtils().ApplyFixity(KM.DISPLACEMENT_X, True, mp.Nodes)
        KM.VariableUtils().ApplyFixity(KM.DISPLACEMENT_Y, True, mp.Nodes)
        KM.VariableUtils().ApplyFixity(KM.DISPLACEMENT_Z, True, mp.Nodes)
        KM.VariableUtils().ApplyFixity(KM.ROTATION_X, True, mp.Nodes)
        KM.VariableUtils().ApplyFixity(KM.ROTATION_Y, True, mp.Nodes)
        KM.VariableUtils().ApplyFixity(KM.ROTATION_Z, True, mp.Nodes)

    def _apply_displacement(self, nodes, displacement):
        value = KM.Vector(displacement)
        for node in nodes:
            node.SetSolutionStepValue(
                KM.DISPLACEMENT,
                0,
                value
            )

    def _apply_rotation(self, nodes, rotation):
        value = KM.Vector(rotation)
        for node in nodes:
            node.SetSolutionStepValue(
                KM.ROTATION,
                0,
                value
            )

    def _apply_material_properties(self,mp):
        #define properties
        mp.CreateNewProperties(1)
        mp.GetProperties()[1].SetValue(KM.YOUNG_MODULUS,210E+3)
        mp.GetProperties()[1].SetValue(KM.POISSON_RATIO,0.3)
        mp.GetProperties()[1].SetValue(KM.THICKNESS,1.0)
        mp.GetProperties()[1].SetValue(KM.DENSITY,7.850E-9)
        cl = StructuralMechanicsApplication.ReissnerMindlinShellElasticConstitutiveLaw()
        mp.GetProperties()[1].SetValue(KM.CONSTITUTIVE_LAW,cl)


    def _solve(self,mp):
        linear_solver = KM.SkylineLUFactorizationSolver()
        builder_and_solver = KM.ResidualBasedEliminationBuilderAndSolverWithConstraints(
            linear_solver
        )
        scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
        convergence_criterion = KM.ResidualCriteria(1e-14,1e-20)

        max_iters = 20
        compute_reactions = True
        reform_step_dofs = True
        calculate_norm_dx = False
        move_mesh_flag = True
        strategy = KM.ResidualBasedLinearStrategy(mp,
                                                  scheme,
                                                  builder_and_solver,
                                                  compute_reactions,
                                                  reform_step_dofs,
                                                  calculate_norm_dx,
                                                  move_mesh_flag)
        strategy.SetEchoLevel(0)
        strategy.Initialize()
        strategy.Check()
        strategy.Solve()

    def _check_slave_disp(self,node,displacement_results, rotation_results):
        #check that the results are exact on the node
        disp = node.GetSolutionStepValue(KM.DISPLACEMENT)
        rot = node.GetSolutionStepValue(KM.ROTATION)
        axis = ["x", "y", "z"]
        for component in range(3):
            self.assertIsClose(disp[component], 
                               displacement_results[component],
                               abs_tol=0.001, 
                               rel_tol=0.0001, 
                               msg= None #wird nur ausgegeben bei Fehler
            )
            #diff = disp[component] - displacement_results[component]
            #print(f"Node {node.Id}, displacement[{axis[component]}]: diff = {diff}")
            self.assertIsClose(
                rot[component],
                rotation_results[component],
                abs_tol=0.001, 
                rel_tol=0.0001, 
                msg= None
            )
            #diff_rot = rot[component] - rotation_results[component]
            #print(f"Node {node.Id}, rotation[{axis[component]}]: diff_rot = {diff_rot}")

    def execute_RBE2_test(self, current_model, element_name, applied_disp, applied_rot, displacement_results, rotation_results):
        mp = current_model.CreateModelPart("Structure")
        model = current_model
        mp.SetBufferSize(2)

        self._add_variables(mp)
        self._apply_material_properties(mp)

        if element_name == "MITCThickShellElement3D4N":
            cl = KM.StructuralMechanicsApplication.ReissnerMindlinShellElasticConstitutiveLaw()
            mp.GetProperties()[1].SetValue(KM.CONSTITUTIVE_LAW, cl)

        self._create_nodes(mp)
        self._add_dofs(mp)
        self._create_elements(mp,element_name)

        #create a submodelpart for slaves
        slave = mp.CreateSubModelPart("RBE2_slaves") #wichtig, dass diese namen mit dem Prozess übereinstimmen
        slave.AddNodes([10,11,12])

        master = mp.CreateSubModelPart("RBE2_master")
        master.AddNodes([13])

        #create a submodelpart for dirichlet boundary conditions
        bcs_dirichlet = mp.CreateSubModelPart("BoundaryCondtionsDirichlet")
        bcs_dirichlet.AddNodes([1,2,3])

        self._apply_dirichlet_BCs(bcs_dirichlet)

        rbe2_settings = KM.Parameters(r"""
        {
            "model_part_name": "Structure",
            "master_sub_model_part": "Structure.RBE2_master",
            "slave_sub_model_part": "Structure.RBE2_slaves",
            "constraint_id_start": 1,
            "constrained_dofs": ["DISPLACEMENT_X", "DISPLACEMENT_Y", "DISPLACEMENT_Z", "ROTATION_X", "ROTATION_Y", "ROTATION_Z"]
        }""")

        rbe2_process = ApplyRbe2Process(
             current_model,
             rbe2_settings
        )

        rbe2_process.ExecuteInitialize()

        master_node = mp.GetNode(13)
        for variable in (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z,
                     KM.ROTATION_X, KM.ROTATION_Y, KM.ROTATION_Z):
            master_node.Fix(variable)

        self._apply_displacement([master_node], applied_disp)
        self._apply_rotation([master_node], applied_rot)
        self._solve(mp)


        for node_id, expected_displacement, expected_rotation in zip(
            [10, 11, 12], displacement_results, rotation_results
        ):
            self._check_slave_disp(
                mp.GetNode(node_id),
                expected_displacement,
                expected_rotation
            )

    def test_RBE2_panel_ux(self):
        self._test_RBE2_panel(
            applied_disp=[-2.0, 0.0, 0.0],
            applied_rot=[0.0, 0.0, 0.0],
            displacement_results=[
                [-2.0, 0.0, 0.0],  # Slave 10
                [-2.0, 0.0, 0.0],  # Slave 11
                [-2.0,0.0, 0.0]   # Slave 12
            ],   
            rotation_results=[                
                [0.0, 0.0, 0.0],  # Slave 10
                [0.0, 0.0, 0.0],  # Slave 11
                [0.0,0.0, 0.0]   # Slave 12],
            ]
        )

    def test_RBE2_panel_uy(self):
        self._test_RBE2_panel(
            applied_disp=[0.0, -2.0, 0.0],
            applied_rot=[0.0, 0.0, 0.0],
            displacement_results=[
                [0.0, -2.0, 0.0],  # Slave 10
                [0.0, -2.0, 0.0],  # Slave 11
                [0.0, -2.0, 0.0]   # Slave 12
            ],   
            rotation_results=[                
                [0.0, 0.0, 0.0],  # Slave 10
                [0.0, 0.0, 0.0],  # Slave 11
                [0.0,0.0, 0.0]   # Slave 12],
            ]
        )  

    def test_RBE2_panel_uz(self):
        self._test_RBE2_panel(
            applied_disp=[0.0, 0.0, -2.0],
            applied_rot=[0.0, 0.0, 0.0],
            displacement_results=[
                [0.0, 0.0, -2.0],  # Slave 10
                [0.0, 0.0, -2.0],  # Slave 11
                [0.0, 0.0, -2.0]   # Slave 12
            ],   
            rotation_results=[                
                [0.0, 0.0, 0.0],  # Slave 10
                [0.0, 0.0, 0.0],  # Slave 11
                [0.0,0.0, 0.0]   # Slave 12],
            ]
        )  


    def test_RBE2_panel_rx(self):
        self._test_RBE2_panel(
            applied_disp=[0.0, 0.0, 0.0],
            applied_rot=[0.01, 0.0, 0.0],
            displacement_results=[
                [0.0,-0.06, 0.0],  # Slave 10
                [0.0,-0.06, 0.0],  # Slave 11
                [0.0,-0.06, 0.0]   # Slave 12
            ],   
            rotation_results=[                
                [0.01, 0.0, 0.0],  # Slave 10
                [0.01, 0.0, 0.0],  # Slave 11
                [0.01,0.0, 0.0]   # Slave 12],
            ]
        )   

    def test_RBE2_panel_ry(self):
        self._test_RBE2_panel(
            applied_disp=[0.0, 0.0, 0.0],
            applied_rot=[0.0, 0.01, 0.0],
            displacement_results=[
                [0.06, 0.0, 0.02],  # Slave 10
                [0.06, 0.0, 0.015],  # Slave 11
                [0.06, 0.0, 0.01]   # Slave 12
            ],   
            rotation_results=[                
                [0.0, 0.01, 0.0],  # Slave 10
                [0.0, 0.01, 0.0],  # Slave 11
                [0.0,0.01, 0.0]   # Slave 12],
            ]
        )   

    def test_RBE2_panel_rz(self):
        self._test_RBE2_panel(
            applied_disp=[0.0, 0.0, 0.0],
            applied_rot=[0.0, 0.0, 0.01],
            displacement_results=[
                [0.0, -0.02, 0.0],  # Slave 10
                [0.0, -0.015, 0.0],  # Slave 11
                [0.0, -0.01, 0.0]   # Slave 12
            ],   
            rotation_results=[                
                [0.0, 0.0, 0.01],  # Slave 10
                [0.0, 0.0, 0.01],  # Slave 11
                [0.0,0.0, 0.01]   # Slave 12],
            ]
        )


    def _test_RBE2_panel(self, applied_disp, applied_rot, displacement_results, rotation_results):
        element_name = "MITCThickShellElement3D4N"
        current_model = KM.Model()
        self.execute_RBE2_test(current_model,
                                element_name,
                                applied_disp,
                                applied_rot,
                                displacement_results,
                                rotation_results) 

if __name__ == "__main__":
    KratosUnittest.main()


