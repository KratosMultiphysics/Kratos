from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Importing post-process
from KratosMultiphysics.vtk_output_process import VtkOutputProcess
from KratosMultiphysics.gid_output_process import GiDOutputProcess

class TestPatchTestSmallDisplacementMixedVolumetricStrain(KratosUnittest.TestCase):
    def setUp(self):
        self.tolerance = 1.0e-6
        self.print_output = True

    def _add_variables(self, ModelPart):
        ModelPart.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        ModelPart.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUMETRIC_STRAIN)
        ModelPart.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)
        ModelPart.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        ModelPart.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_FLUX)

    def _apply_BCs(self, ModelPart, A, b):
        for node in ModelPart.Nodes:
            node.Fix(KratosMultiphysics.DISPLACEMENT_X)
            node.Fix(KratosMultiphysics.DISPLACEMENT_Y)
            node.Fix(KratosMultiphysics.DISPLACEMENT_Z)

            x_vec = KratosMultiphysics.Vector(3)
            x_vec[0] = node.X0
            x_vec[1] = node.Y0
            x_vec[2] = node.Z0

            u = A*x_vec
            u += b

            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, 0, u)

    def _apply_material_properties(self, ModelPart, Dimension):
        # Define material properties
        ModelPart.GetProperties()[1].SetValue(KratosMultiphysics.YOUNG_MODULUS, 200.0e9)
        ModelPart.GetProperties()[1].SetValue(KratosMultiphysics.POISSON_RATIO, 0.4)
        ModelPart.GetProperties()[1].SetValue(KratosMultiphysics.DENSITY, 1.0)

        # Define body force
        g = [0,0,0]
        ModelPart.GetProperties()[1].SetValue(KratosMultiphysics.VOLUME_ACCELERATION, g)

        # Define constitutive law
        if(Dimension == 2):
            cons_law = StructuralMechanicsApplication.LinearElasticPlaneStrain2DLaw()
        else:
            cons_law = StructuralMechanicsApplication.LinearElastic3DLaw()
        ModelPart.GetProperties()[1].SetValue(KratosMultiphysics.CONSTITUTIVE_LAW, cons_law)

    def _define_movement(self, Dimension):
        if(Dimension == 2):
            #define the applied motion - the idea is that the displacement is defined as u = A*xnode + b
            #so that the displcement is linear and the exact F = I + A
            A = KratosMultiphysics.Matrix(3,3)
            A[0,0] = 1.0e-10;  A[0,1] = 2.0e-10; A[0,2] = 0.0
            A[1,0] = 0.5e-10;  A[1,1] = 0.7e-10; A[1,2] = 0.0
            A[2,0] = 0.0;      A[2,1] = 0.0;     A[2,2] = 0.0

            b = KratosMultiphysics.Vector(3)
            b[0] = 0.5e-10
            b[1] = -0.2e-10
            b[2] = 0.0
        else:
            #define the applied motion - the idea is that the displacement is defined as u = A*xnode + b
            #so that the displcement is linear and the exact F = I + A
            A = KratosMultiphysics.Matrix(3,3)
            A[0,0] = 1.0e-10;   A[0,1] = 2.0e-10; A[0,2] = 0.0
            A[1,0] = 0.5e-10;   A[1,1] = 0.7e-10; A[1,2] = 0.1e-10
            A[2,0] = -0.2e-10;  A[2,1] = 0.0;     A[2,2] = -0.3e-10

            b = KratosMultiphysics.Vector(3)
            b[0] = 0.5e-10
            b[1] = -0.2e-10
            b[2] = 0.7e-10

        return A,b

    def _solve(self, ModelPart):
        # Define a linear strategy to solve the problem
        linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()
        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
        scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        compute_reactions = True
        reform_step_dofs = True
        calculate_norm_dx = False
        move_mesh_flag = True

        strategy = KratosMultiphysics.ResidualBasedLinearStrategy(
            ModelPart,
            scheme,
            linear_solver,
            builder_and_solver,
            compute_reactions,
            reform_step_dofs,
            calculate_norm_dx,
            move_mesh_flag)
        strategy.SetEchoLevel(0)
        strategy.Check()

        # Solve the problem
        strategy.Solve()

    def _check_results(self, ModelPart, A, b):
        # Check that the results are exact on  the nodes
        for node in ModelPart.Nodes:
            x_vec = KratosMultiphysics.Vector(3)
            x_vec[0] = node.X0
            x_vec[1] = node.Y0
            x_vec[2] = node.Z0

            u = A*x_vec
            u += b

            coor_list = ["X","Y","Z"]

            d = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)
            for i in range(3):
                if abs(u[i]) > 0.0:
                    error = (d[i] - u[i])/u[i]
                    if error > self.tolerance:
                       print("NODE ", node.Id,": Component ", coor_list[i],":\t",u[i],"\t",d[i], "\tError: ", error)
                    self.assertLess(error, self.tolerance)

    def testSmallDisplacementMixedVolumetricStrainElement2DTriangle(self):
        dimension = 2
        current_model = KratosMultiphysics.Model()
        model_part = current_model.CreateModelPart("MainModelPartTriangle")
        self._add_variables(model_part)
        self._apply_material_properties(model_part, dimension)

        # Create nodes
        model_part.CreateNewNode(1,0.5,0.5,0.0)
        model_part.CreateNewNode(2,0.7,0.2,0.0)
        model_part.CreateNewNode(3,0.9,0.8,0.0)
        model_part.CreateNewNode(4,0.3,0.7,0.0)
        model_part.CreateNewNode(5,0.6,0.6,0.0)

        # Add DOFs
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X, model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y, model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z, model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.VOLUMETRIC_STRAIN, KratosMultiphysics.REACTION_FLUX, model_part)

        # Create a submodelpart for boundary conditions
        boundary_model_part = model_part.CreateSubModelPart("BoundaryCondtions")
        boundary_model_part.AddNodes([1,2,3,4])

        # Create elements
        model_part.CreateNewElement("SmallDisplacementMixedVolumetricStrainElement2D3N", 1, [1,2,5], model_part.GetProperties()[1])
        model_part.CreateNewElement("SmallDisplacementMixedVolumetricStrainElement2D3N", 2, [2,3,5], model_part.GetProperties()[1])
        model_part.CreateNewElement("SmallDisplacementMixedVolumetricStrainElement2D3N", 3, [3,4,5], model_part.GetProperties()[1])
        model_part.CreateNewElement("SmallDisplacementMixedVolumetricStrainElement2D3N", 4, [4,1,5], model_part.GetProperties()[1])

        A,b = self._define_movement(dimension)

        self._apply_BCs(boundary_model_part, A, b)
        self._solve(model_part)
        self._check_results(model_part, A, b)
        if self.print_output:
            self.__post_process(model_part)

    def testSmallDisplacementMixedVolumetricStrainElement3DTetrahedra(self):
        dimension = 3
        current_model = KratosMultiphysics.Model()
        model_part = current_model.CreateModelPart("MainModelPartTetrahedra")
        self._add_variables(model_part)
        self._apply_material_properties(model_part, dimension)

        #create nodes
        model_part.CreateNewNode(1,0.0, 1.0, 0.0)
        model_part.CreateNewNode(2,0.0, 1.0, 0.1)
        model_part.CreateNewNode(3, 0.28739360416666665, 0.27808503701741405, 0.05672979583333333)
        model_part.CreateNewNode(4, 0.0, 0.1, 0.0)
        model_part.CreateNewNode(5, 0.1, 0.1, 0.1)
        model_part.CreateNewNode(6, 1.0, 0.0, 0.0)
        model_part.CreateNewNode(7, 1.2, 0.0, 0.1)

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.VOLUMETRIC_STRAIN, KratosMultiphysics.REACTION_FLUX, model_part)

        #create a submodelpart for boundary conditions
        boundary_model_part = model_part.CreateSubModelPart("BoundaryCondtions")
        boundary_model_part.AddNodes([1,2,4,5,6,7])

        #create Element
        model_part.CreateNewElement("SmallDisplacementMixedVolumetricStrainElement3D4N", 1,[5,3,1,2], model_part.GetProperties()[1])
        model_part.CreateNewElement("SmallDisplacementMixedVolumetricStrainElement3D4N", 2,[3,1,2,6], model_part.GetProperties()[1])
        model_part.CreateNewElement("SmallDisplacementMixedVolumetricStrainElement3D4N", 3,[6,4,7,3], model_part.GetProperties()[1])
        model_part.CreateNewElement("SmallDisplacementMixedVolumetricStrainElement3D4N", 4,[5,4,1,3], model_part.GetProperties()[1])
        model_part.CreateNewElement("SmallDisplacementMixedVolumetricStrainElement3D4N", 5,[4,1,3,6], model_part.GetProperties()[1])
        model_part.CreateNewElement("SmallDisplacementMixedVolumetricStrainElement3D4N", 6,[5,4,3,7], model_part.GetProperties()[1])
        model_part.CreateNewElement("SmallDisplacementMixedVolumetricStrainElement3D4N", 7,[3,5,7,2], model_part.GetProperties()[1])
        model_part.CreateNewElement("SmallDisplacementMixedVolumetricStrainElement3D4N", 8,[6,7,2,3], model_part.GetProperties()[1])

        A,b = self._define_movement(dimension)

        self._apply_BCs(boundary_model_part, A, b)
        self._solve(model_part)
        self._check_results(model_part, A, b)
        if self.print_output:
            self.__post_process(model_part)

    def __post_process(self, main_model_part, post_type = "gid"):
        if post_type == "gid":
            self.gid_output = GiDOutputProcess(
                main_model_part,
                main_model_part.Name,
                KratosMultiphysics.Parameters(r"""
                {
                    "result_file_configuration" : {
                    "gidpost_flags": {
                        "GiDPostMode": "GiD_PostBinary",
                        "WriteDeformedMeshFlag": "WriteUndeformed",
                        "WriteConditionsFlag": "WriteConditions",
                        "MultiFileFlag": "SingleFile"
                    },
                    "nodal_results"       : ["DISPLACEMENT","VOLUMETRIC_STRAIN"],
                    "gauss_point_results" : []
                    }
                }"""))

            self.gid_output.ExecuteInitialize()
            self.gid_output.ExecuteBeforeSolutionLoop()
            self.gid_output.ExecuteInitializeSolutionStep()
            self.gid_output.PrintOutput()
            self.gid_output.ExecuteFinalizeSolutionStep()
            self.gid_output.ExecuteFinalize()

        elif post_type == "vtk":
            vtk_output_parameters = KratosMultiphysics.Parameters(r"""
            {
                "model_part_name": "",
                "extrapolate_gauss_points": false,
                "nodal_solution_step_data_variables" : ["DISPLACEMENT","VOLUMETRIC_STRAIN"],
                "gauss_point_variables": []
            }""")
            vtk_output_parameters["model_part_name"].SetString(main_model_part.Name)
            self.vtk_output_process = VtkOutputProcess(
                main_model_part.GetModel(),
                vtk_output_parameters)

            self.vtk_output_process.ExecuteInitialize()
            self.vtk_output_process.ExecuteBeforeSolutionLoop()
            self.vtk_output_process.ExecuteInitializeSolutionStep()
            self.vtk_output_process.PrintOutput()
            self.vtk_output_process.ExecuteFinalizeSolutionStep()
            self.vtk_output_process.ExecuteFinalize()

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()