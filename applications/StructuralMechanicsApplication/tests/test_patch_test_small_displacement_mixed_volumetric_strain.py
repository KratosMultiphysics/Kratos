import KratosMultiphysics

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Importing post-process
from KratosMultiphysics.vtk_output_process import VtkOutputProcess
from KratosMultiphysics.gid_output_process import GiDOutputProcess

class TestPatchTestSmallDisplacementMixedVolumetricStrain(KratosUnittest.TestCase):
    def setUp(self):
        self.tolerances = {
            "relative" : 1e-6,
            "absolute" : {
                "displacement" : 1e-6,
                "stress"       : 1e2,
                "strain"       : 1e-2
            }
        }
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

    def _apply_user_provided_material_properties(self, ModelPart, Dimension):
        # Define body force
        g = [0,0,0]
        ModelPart.GetProperties()[1].SetValue(KratosMultiphysics.VOLUME_ACCELERATION, g)

        # Define constitutive law
        if(Dimension == 2):
            cons_law = StructuralMechanicsApplication.UserProvidedLinearElastic2DLaw()
        else:
            cons_law = StructuralMechanicsApplication.UserProvidedLinearElastic3DLaw()
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
        scheme = StructuralMechanicsApplication.StructuralMechanicsStaticScheme(KratosMultiphysics.Parameters("{}"))
        compute_reactions = True
        reform_step_dofs = True
        calculate_norm_dx = False
        move_mesh_flag = True

        strategy = KratosMultiphysics.ResidualBasedLinearStrategy(
            ModelPart,
            scheme,
            builder_and_solver,
            compute_reactions,
            reform_step_dofs,
            calculate_norm_dx,
            move_mesh_flag)
        strategy.SetEchoLevel(0)
        strategy.Initialize()
        strategy.Check()

        # Solve the problem
        strategy.Solve()

    def _solve_non_linear(self, model_part, projection_variables_list = []):
        # Define a Newton-Raphson strategy to solve the problem
        linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()
        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
        scheme_settings = KratosMultiphysics.Parameters("""{
            "projection_variables_list" : []
        }""")
        scheme_settings["projection_variables_list"].SetStringArray(projection_variables_list)
        scheme = StructuralMechanicsApplication.StructuralMechanicsStaticScheme(scheme_settings)
        convergence_criterion = KratosMultiphysics.MixedGenericCriteria(
            [(KratosMultiphysics.DISPLACEMENT, 1.0e-6, 1.0e-8),
            (KratosMultiphysics.VOLUMETRIC_STRAIN, 1.0e-6, 1.0e-8)])
        max_iteration = 10
        compute_reactions = False
        reform_step_dofs = False
        move_mesh_flag = True
        strategy = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(
            model_part,
            scheme,
            convergence_criterion,
            builder_and_solver,
            max_iteration,
            compute_reactions,
            reform_step_dofs,
            move_mesh_flag)
        strategy.SetEchoLevel(0)
        strategy.Initialize()
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
                    error = abs((d[i] - u[i])/u[i])
                    self.assertLess(error, self.tolerances["relative"], msg=f"NODE {node.Id}: Component {coor_list[i]}: {u[i]} {d[i]} Error: {error}")
                else:
                    self.assertLess(abs(d[i]), self.tolerances["absolute"]["displacement"])

    def _calculate_reference_strain(self, A, dim):
        # Given the matrix A, the analytic deformation gradient is F+I
        F = A
        for i in range(3):
            F[i,i] += 1.0

        # Here compute the Cauchy Green strain tensor
        cauchy_green_strain_tensor = KratosMultiphysics.Matrix(3,3)

        for i in range(3):
            for j in range(3):
                cauchy_green_strain_tensor[i,j] = 0.0

        for i in range(3):
            for j in range(3):
                for k in range(3):
                    cauchy_green_strain_tensor[i,j] += F[k,i]*F[k,j]

        for i in range(3):
            cauchy_green_strain_tensor[i,i] -= 1.0

        for i in range(3):
            for j in range(3):
                cauchy_green_strain_tensor[i,j] = 0.5*cauchy_green_strain_tensor[i,j]

        # Cauchy Green strain tensor in Voigt notation
        if(dim == 2):
            reference_strain = KratosMultiphysics.Vector(3)
            reference_strain[0] = cauchy_green_strain_tensor[0,0]
            reference_strain[1] = cauchy_green_strain_tensor[1,1]
            reference_strain[2] = 2.0*cauchy_green_strain_tensor[0,1]
        else:
            reference_strain = KratosMultiphysics.Vector(6)
            reference_strain[0] = cauchy_green_strain_tensor[0,0]
            reference_strain[1] = cauchy_green_strain_tensor[1,1]
            reference_strain[2] = cauchy_green_strain_tensor[2,2]
            reference_strain[3] = 2.0*cauchy_green_strain_tensor[0,1]
            reference_strain[4] = 2.0*cauchy_green_strain_tensor[1,2]
            reference_strain[5] = 2.0*cauchy_green_strain_tensor[0,2]

        return reference_strain

    def _check_stress(self, model_part, A, dim):
        # Calculate the reference strain
        reference_strain = self._calculate_reference_strain(A, dim)

        young = model_part.GetProperties()[1].GetValue(KratosMultiphysics.YOUNG_MODULUS)
        poisson = model_part.GetProperties()[1].GetValue(KratosMultiphysics.POISSON_RATIO)

        # Finally compute stress
        if(dim == 2):
            #here assume plane stress
            c1 = young / (1.00 - poisson*poisson)
            c2 = c1 * poisson
            c3 = 0.5* young / (1 + poisson)
            reference_stress = KratosMultiphysics.Vector(3)
            reference_stress[0] = c1*reference_strain[0] + c2 * (reference_strain[1])
            reference_stress[1] = c1*reference_strain[1] + c2 * (reference_strain[0])
            reference_stress[2] = c3*reference_strain[2]
        else:
            c1 = young / (( 1.00 + poisson ) * ( 1 - 2 * poisson ) )
            c2 = c1 * ( 1 - poisson )
            c3 = c1 * poisson
            c4 = c1 * 0.5 * ( 1 - 2 * poisson )
            reference_stress = KratosMultiphysics.Vector(6)
            reference_stress[0] = c2*reference_strain[0] + c3 * (reference_strain[1] + reference_strain[2])
            reference_stress[1] = c2*reference_strain[1] + c3 * (reference_strain[0] + reference_strain[2])
            reference_stress[2] = c2*reference_strain[2] + c3 * (reference_strain[0] + reference_strain[1])
            reference_stress[3] = c4*reference_strain[3]
            reference_stress[4] = c4*reference_strain[4]
            reference_stress[5] = c4*reference_strain[5]

        for elem in model_part.Elements:
            out = elem.CalculateOnIntegrationPoints(KratosMultiphysics.PK2_STRESS_VECTOR, model_part.ProcessInfo)
            for stress in out:
                for i in range(len(reference_stress)):
                    if abs(stress[i]) > 0.0:
                        self.assertLess((reference_stress[i] - stress[i])/stress[i], self.tolerances["relative"])
                    else:
                        self.assertLess(abs(stress[i]), self.tolerances["absolute"]["stress"])

    def _check_stress_user_provided(self, model_part, A, dim):
        # Calculate the reference strain
        reference_strain = self._calculate_reference_strain(A, dim)

        # Finally compute stress
        elasticity_tensor = model_part.GetProperties()[1].GetValue(StructuralMechanicsApplication.ELASTICITY_TENSOR)
        reference_stress = elasticity_tensor * reference_strain

        for elem in model_part.Elements:
            out = elem.CalculateOnIntegrationPoints(KratosMultiphysics.PK2_STRESS_VECTOR, model_part.ProcessInfo)
            for stress in out:
                for i in range(len(reference_stress)):
                    if abs(stress[i]) > 0.0:
                        self.assertLess((reference_stress[i] - stress[i])/stress[i], self.tolerances["relative"])
                    else:
                        self.assertLess(abs(stress[i]), self.tolerances["absolute"]["stress"])

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
        self._check_stress(model_part, A, dimension)
        if self.print_output:
            self.__post_process(model_part)

    def testSmallDisplacementMixedVolumetricStrainOssElement2DTriangle(self):
        dimension = 2
        current_model = KratosMultiphysics.Model()
        model_part = current_model.CreateModelPart("MainModelPartTriangle")
        self._add_variables(model_part)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT_PROJECTION)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUMETRIC_STRAIN_PROJECTION)
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
        model_part.CreateNewElement("SmallDisplacementMixedVolumetricStrainOssElement2D3N", 1, [1,2,5], model_part.GetProperties()[1])
        model_part.CreateNewElement("SmallDisplacementMixedVolumetricStrainOssElement2D3N", 2, [2,3,5], model_part.GetProperties()[1])
        model_part.CreateNewElement("SmallDisplacementMixedVolumetricStrainOssElement2D3N", 3, [3,4,5], model_part.GetProperties()[1])
        model_part.CreateNewElement("SmallDisplacementMixedVolumetricStrainOssElement2D3N", 4, [4,1,5], model_part.GetProperties()[1])

        A,b = self._define_movement(dimension)

        self._apply_BCs(boundary_model_part, A, b)
        self._solve_non_linear(model_part, ["DISPLACEMENT_PROJECTION", "VOLUMETRIC_STRAIN_PROJECTION"])
        self._check_results(model_part, A, b)
        self._check_stress(model_part, A, dimension)
        if self.print_output:
            self.__post_process(model_part)

    def testSmallDisplacementMixedVolumetricStrainOssNonLinearElement2DTriangle(self):
        dimension = 2
        current_model = KratosMultiphysics.Model()
        model_part = current_model.CreateModelPart("MainModelPartTriangle")
        self._add_variables(model_part)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT_PROJECTION)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUMETRIC_STRAIN_PROJECTION)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT_PROJECTION_REACTION)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUMETRIC_STRAIN_PROJECTION_REACTION)
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
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_PROJECTION_X, KratosMultiphysics.DISPLACEMENT_PROJECTION_REACTION_X, model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_PROJECTION_Y, KratosMultiphysics.DISPLACEMENT_PROJECTION_REACTION_Y, model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_PROJECTION_Z, KratosMultiphysics.DISPLACEMENT_PROJECTION_REACTION_Z, model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.VOLUMETRIC_STRAIN_PROJECTION, KratosMultiphysics.VOLUMETRIC_STRAIN_PROJECTION_REACTION, model_part)

        # Create a submodelpart for boundary conditions
        boundary_model_part = model_part.CreateSubModelPart("BoundaryCondtions")
        boundary_model_part.AddNodes([1,2,3,4])

        # Create elements
        model_part.CreateNewElement("SmallDisplacementMixedVolumetricStrainOssNonLinearElement2D3N", 1, [1,2,5], model_part.GetProperties()[1])
        model_part.CreateNewElement("SmallDisplacementMixedVolumetricStrainOssNonLinearElement2D3N", 2, [2,3,5], model_part.GetProperties()[1])
        model_part.CreateNewElement("SmallDisplacementMixedVolumetricStrainOssNonLinearElement2D3N", 3, [3,4,5], model_part.GetProperties()[1])
        model_part.CreateNewElement("SmallDisplacementMixedVolumetricStrainOssNonLinearElement2D3N", 4, [4,1,5], model_part.GetProperties()[1])

        A,b = self._define_movement(dimension)

        self._apply_BCs(boundary_model_part, A, b)
        self._solve_non_linear(model_part)
        self._check_results(model_part, A, b)
        self._check_stress(model_part, A, dimension)
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
        self._check_stress(model_part, A, dimension)
        if self.print_output:
            self.__post_process(model_part)

    def testSmallDisplacementMixedVolumetricStrainElement3DTetrahedraAnisotropic(self):
        dimension = 3
        current_model = KratosMultiphysics.Model()
        model_part = current_model.CreateModelPart("MainModelPartTetrahedra")
        self._add_variables(model_part)
        self._apply_user_provided_material_properties(model_part, dimension)

        #set the user-provided anisotropic elasticity tensor
        elasticity_tensor = KratosMultiphysics.Matrix(6,6)
        # for i in range(6):
        #     for j in range(6):
        #         elasticity_tensor[i,j] = 0.0
        # E = 2.0e11
        # NU = 0.3
        # c1 = E / (( 1.00 + NU ) * ( 1 - 2 * NU ) );
        # c2 = c1 * ( 1 - NU );
        # c3 = c1 * NU;
        # c4 = c1 * 0.5 * ( 1 - 2 * NU );
        # elasticity_tensor[0, 0] = c2;
        # elasticity_tensor[0, 1] = c3;
        # elasticity_tensor[0, 2] = c3;
        # elasticity_tensor[1, 0] = c3;
        # elasticity_tensor[1, 1] = c2;
        # elasticity_tensor[1, 2] = c3;
        # elasticity_tensor[2, 0] = c3;
        # elasticity_tensor[2, 1] = c3;
        # elasticity_tensor[2, 2] = c2;
        # elasticity_tensor[3, 3] = c4;
        # elasticity_tensor[4, 4] = c4;
        # elasticity_tensor[5, 5] = c4;
        aux_elasticity_tensor = [
            [5.99E+11,5.57E+11,5.34E+11,0,0,4.44E+09],
            [5.57E+11,5.71E+11,5.34E+11,0,0,-3.00E+09],
            [5.34E+11,5.34E+11,5.37E+11,0,0,9.90E+05],
            [0,0,0,1.92E+09,9.78E+06,0],
            [0,0,0,9.78E+06,2.12E+09,0],
            [4.44E+09,-3.00E+09,9.90E+05,0,0,2.56E+10]]
        for i in range(6):
            for j in range(6):
                elasticity_tensor[i,j] = aux_elasticity_tensor[i][j]

        # T = KratosMultiphysics.Matrix(6,6)
        # T.fill(0.0)
        # aux_E = 0.0
        # aux_G = 0.0
        # for i in range(6):
        #     if i < 3:
        #         aux_E += elasticity_tensor[i,i]
        #     else:
        #         aux_G += elasticity_tensor[i,i]
        # aux_E /= 3.0
        # aux_G /= 3.0
        # for i in range(6):
        #     if i < 3:
        #         T[i,i] = math.sqrt(aux_E/elasticity_tensor[i,i])
        #     else:
        #         T[i,i] = math.sqrt(aux_G/elasticity_tensor[i,i])

        # elasticity_tensor = T * elasticity_tensor * T

        model_part.GetProperties()[1].SetValue(StructuralMechanicsApplication.ELASTICITY_TENSOR, elasticity_tensor)

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
        if self.print_output:
            self.__post_process(model_part)
        self._check_results(model_part, A, b)
        self._check_stress_user_provided(model_part, A, dimension)

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
                    "gauss_point_results" : ["CAUCHY_STRESS_VECTOR"]
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
                "gauss_point_variables": ["CAUCHY_STRESS_VECTOR"]
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