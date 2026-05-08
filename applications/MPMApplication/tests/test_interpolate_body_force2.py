import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.MPMApplication as KratosMPM
from KratosMultiphysics.MPMApplication.mpm_analysis import MpmAnalysis

class TestBodyForceInterpolationMPM(KratosUnittest.TestCase):

    def setUp(self):
        super().setUp()
        self.model = KratosMultiphysics.Model()
        
        # 1. Crear ModelParts principales
        mpm_mp = self.model.CreateModelPart("MPM_Material")
        initial_mesh_mp = self.model.CreateModelPart("Initial_MPM_Material")
        grid_mp = self.model.CreateModelPart("Background_Grid")
            
        # 3. Configurar Malla Inicial de Material (Cuerpo)
        initial_mesh_mp.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)
        initial_mesh_sub_model_part = initial_mesh_mp.CreateSubModelPart("SubInitialMesh")

        
        # Propiedades del material
        prop = initial_mesh_sub_model_part.GetProperties()[1]
        prop.SetValue(KratosMPM.MATERIAL_POINTS_PER_ELEMENT, 1)
        prop.SetValue(KratosMultiphysics.DENSITY, 7850.0)
        prop.SetValue(KratosMultiphysics.YOUNG_MODULUS, 206900000000.0)
        prop.SetValue(KratosMultiphysics.POISSON_RATIO, 0.29)
        
        # Nodos del cuerpo
        initial_mesh_sub_model_part.CreateNewNode(5, 1.0, 0.0, 0.0)
        initial_mesh_sub_model_part.CreateNewNode(7, 1.0, 1.0, 0.0)
        initial_mesh_sub_model_part.CreateNewNode(10, 0.0, 0.0, 0.0)
        initial_mesh_sub_model_part.CreateNewNode(12, 0.0, 1.0, 0.0)
        
        # Elemento del cuerpo
        initial_mesh_sub_model_part.CreateNewElement("Element2D4N", 4, [5, 7, 12, 10], prop)
        
        # 4. Configurar Malla de Fondo (Grid)
        grid_mp.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)
        grid_sub_model_part = grid_mp.CreateSubModelPart("SubBackgroundGrid")
        
        # Nodos del grid
        grid_sub_model_part.CreateNewNode(1, 3.0, 0.0, 0.0)
        grid_sub_model_part.CreateNewNode(2, 3.0, 1.0, 0.0)
        grid_sub_model_part.CreateNewNode(3, 2.0, 0.0, 0.0)
        grid_sub_model_part.CreateNewNode(4, 2.0, 1.0, 0.0)
        grid_sub_model_part.CreateNewNode(6, 1.0, 0.0, 0.0)
        grid_sub_model_part.CreateNewNode(8, 1.0, 1.0, 0.0)
        grid_sub_model_part.CreateNewNode(9, 0.0, 0.0, 0.0)
        grid_sub_model_part.CreateNewNode(11, 0.0, 1.0, 0.0)
        
        # Elementos del grid
        grid_prop = grid_sub_model_part.GetProperties()[1]
        grid_sub_model_part.CreateNewElement("Element2D4N", 1, [6, 8, 11, 9], grid_prop)
        grid_sub_model_part.CreateNewElement("Element2D4N", 2, [3, 4, 8, 6], grid_prop)
        grid_sub_model_part.CreateNewElement("Element2D4N", 3, [1, 2, 4, 3], grid_prop)
        
        # 5. Generación de Puntos Materiales y Condiciones
        mpm_mp.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)
        mpm_mp.SetNodes(grid_mp.GetNodes())
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.ACTIVE, True, initial_mesh_mp.Elements)
        
        KratosMPM.GenerateMaterialPointElement(grid_mp, initial_mesh_mp, mpm_mp, False)
        KratosMPM.GenerateMaterialPointCondition(grid_mp, initial_mesh_mp, mpm_mp)
        
        
        print(f"\nSubModelParts de {mpm_mp.Name}:")
        for sub_mp in mpm_mp.SubModelParts:
            print(f" - {sub_mp.Name}")
        print("\n--- INSPECCIÓN DE MPM_MATERIAL ---")
        print(f"Número de nodos: {mpm_mp.NumberOfNodes()}")
        print(f"Número de elementos (Material Points): {mpm_mp.NumberOfElements()}")
        print(f"Número de condiciones: {mpm_mp.NumberOfConditions()}")
 
    def test_body_force_interpolation(self):
        # 1. Definir parámetros internamente (sustituye a ProjectParameters.json)
        project_parameters = KratosMultiphysics.Parameters("""{
            "problem_data"     : {
                "problem_name"  : "test_body_force",
                "parallel_type" : "OpenMP",
                "echo_level"    : 0,
                "start_time"    : 0.0,
                "end_time"      : 0.1
            },
            "solver_settings"  : {
                "time_stepping"                   : { "time_step" : 0.1 },
                "solver_type"                     : "Dynamic",
                "model_part_name"                 : "MPM_Material",
                "domain_size"                     : 2,
                "echo_level"                      : 0,
                "analysis_type"                   : "non_linear",
                "time_integration_method"         : "implicit",
                "scheme_type"                     : "newmark",
                "compute_reactions"               : false,
                "convergence_criterion"           : "residual_criterion",
                "displacement_relative_tolerance" : 0.0001,
                "displacement_absolute_tolerance" : 1e-9,
                "residual_relative_tolerance"     : 0.0001,
                "residual_absolute_tolerance"     : 1e-9,
                "max_iteration"                   : 10,
                "linear_solver_settings"          : { "solver_type" : "LinearSolversApplication.sparse_lu" },
                "auxiliary_variables_list"        : ["NORMAL","IS_STRUCTURE","BODY_FORCE"]
            },
            "processes" : { "gravity" : [] }
        }""")


        # 3. Ejecutar simulación
        simulation = BodyForceTestSimulation(self.model, project_parameters)
        simulation.Run()

        # --- Get MP model part ---
        mp_model_part = self.model.GetModelPart("MPM_Material")

        self.assertEqual(len(mp_model_part.Elements), 1)

        # Take first MP
        mp = next(iter(mp_model_part.Elements))
        variable_disp = KratosMultiphysics.KratosGlobals.GetVariable( "MP_DISPLACEMENT" )
        variable_acce = KratosMultiphysics.KratosGlobals.GetVariable( "MP_VOLUME_ACCELERATION" )
        displacement = mp.CalculateOnIntegrationPoints(variable_disp,mp_model_part.ProcessInfo)[0]
        acceleration = mp.CalculateOnIntegrationPoints(variable_acce,mp_model_part.ProcessInfo)[0]

        print(acceleration)
        print(displacement)


        # Expected (constant field)
        self.assertAlmostEqual(displacement[0], 0.005, places=12)
        self.assertAlmostEqual(displacement[1], -0.0075, places=12)
        self.assertAlmostEqual(displacement[2], 0.0, places=12)


class BodyForceTestSimulation(MpmAnalysis):

    def __init__(self, model, parameters):
        self.model = model
        self.project_parameters = parameters
   
        self.project_parameters["problem_data"]["end_time"].SetDouble(0.1)
        self.project_parameters["solver_settings"]["time_stepping"]["time_step"].SetDouble(0.1)

        super().__init__(self.model, self.project_parameters)

    def ModifyBeforeSolutionLoop(self):
        super().ModifyBeforeSolutionLoop()
        grid_mp = self.model.CreateModelPart("Background_Grid")
        # Aplicar fuerza de cuerpo constante en los nodos del grid
        for node in grid_mp.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.BODY_FORCE_X, 2.0)
            node.SetSolutionStepValue(KratosMultiphysics.BODY_FORCE_Y, -3.0)
            node.SetSolutionStepValue(KratosMultiphysics.BODY_FORCE_Z, 0.0)

if __name__ == '__main__':
    KratosUnittest.main()
