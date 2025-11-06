import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartOperation
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import IsSameContainerExpression
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import ContainerExpressionTypes, SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.filtering.filter import Factory as FilterFactory
import numpy as np
from scipy.spatial import KDTree


def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> Control:
    if not parameters.Has("name"):
        raise RuntimeError(f"LevelSetControlHJ instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"LevelSetControlHJ instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return LevelSetControlShape(parameters["name"].GetString(), model, parameters["settings"], optimization_problem)


class LevelSetControlShape(Control):
    def __init__(self, name: str, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> None:
        super().__init__(name)

        self.parameters = parameters
        self.model = model
        self.optimization_problem = optimization_problem

        default_settings = Kratos.Parameters("""{
            "controlled_model_part_names"       : [""],
            "output_all_fields"                 : true,
            "echo_level"                        : 0,
            "consider_recursive_property_update": false,
            "filter_settings"                   : {},
            "density"                           : 1.0,
            "young_modulus"                     : 1.0,
            "k"                                 : 400.0,
            "initial_phi_value"                 : 0.0,
            "mininum_modulus"                   : 1000000,
            "reinit_interval"                   : 20
        }""")


        #fill in missing info from defaults
        parameters.ValidateAndAssignDefaults(default_settings)

        self.controlled_physical_variables = [Kratos.DENSITY, Kratos.YOUNG_MODULUS]
        self.output_all_fields = parameters["output_all_fields"].GetBool()
        self.echo_level = parameters["echo_level"].GetInt()
        self.density = parameters["density"].GetDouble()
        self.young_modulus = parameters["young_modulus"].GetDouble()
        self.k = parameters["k"].GetDouble()
        self.initial_phi = parameters["initial_phi_value"].GetDouble()
        self.epsilon = parameters["mininum_modulus"].GetDouble()
        self.reinit_interval = parameters["reinit_interval"].GetInt()

        controlled_model_names_parts = parameters["controlled_model_part_names"].GetStringArray()
        if len(controlled_model_names_parts) == 0:
            raise RuntimeError(f"No model parts are provided for LevelSetControl. [ control name = \"{self.GetName()}\"]")
        self.model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"control_{self.GetName()}", controlled_model_names_parts, False)
        # self.model_part: 'typing.Optional[Kratos.ModelPart]' = None

        self.consider_recursive_property_update = parameters["consider_recursive_property_update"].GetBool()

        self.filter = FilterFactory(self.model, self.model_part_operation.GetModelPartFullName(), Kratos.NODAL_AREA, Kratos.Globals.DataLocation.NodeNonHistorical, parameters["filter_settings"])

    def Initialize(self) -> None:
        self.un_buffered_data = ComponentDataView(self, self.optimization_problem).GetUnBufferedData()

        # init model parts
        self.model_part = self.model_part_operation.GetModelPart()

        # # Creating element specific properties for Youngs Modulus and Density
        if not KratosOA.OptAppModelPartUtils.CheckModelPartStatus(self.model_part, "element_specific_properties_created"):
            KratosOA.OptimizationUtils.CreateEntitySpecificPropertiesForContainer(self.model_part, self.model_part.Elements, self.consider_recursive_property_update)
            KratosOA.OptAppModelPartUtils.LogModelPartStatus(self.model_part, "element_specific_properties_created")

        self.filter.SetComponentDataView(ComponentDataView(self, self.optimization_problem))
        self.filter.Initialize()

        # check if the density or Young's modulus is defined
        element: Kratos.Element
        for element in self.model_part.Elements:
            break
        is_density_defined = element.Properties.Has(Kratos.DENSITY)
        is_youngs_modulus_defined = element.Properties.Has(Kratos.YOUNG_MODULUS)
        if is_density_defined:
            if is_youngs_modulus_defined:
                Kratos.Logger.PrintWarning(self.__class__.__name__, f"Elements of {self.model_part.FullName()} defines both DENSITY and YOUNG_MODULUS. Using DENSITY for initial field calculation and ignoring YOUNG_MODULUS.")
            density = Kratos.Expression.ElementExpression(self.model_part)
            KratosOA.PropertiesVariableExpressionIO.Read(density, Kratos.DENSITY)
        elif is_youngs_modulus_defined:
            young_modulus = Kratos.Expression.ElementExpression(self.model_part)
            KratosOA.PropertiesVariableExpressionIO.Read(young_modulus, Kratos.YOUNG_MODULUS)
        else:
            raise RuntimeError(f"Elements of {self.model_part.FullName()} does not define either DENSITY or YOUNG_MODULUS.")


        # # get the control field
        # self.control_phi = self.GetEmptyField()
        # Kratos.Expression.LiteralExpressionIO.SetData(self.control_phi, self.initial_phi)
        # initialize as SDF
        self._InitializeSignedDistanceKDTree(min_dirac=0.1)

        # get element sizes (no remeshing)
        self.element_len = self._GetElementLength()

        # initialize physical fields 
        self._UpdateAndOutputFields(self.GetEmptyField())


    def _InitializeSignedDistanceKDTree(self, interface_points=None, min_dirac=0.1):
        """
        Build a signed distance function φ using a KDTree for unstructured FEM nodes.
        If interface_points=None, assumes interface at x=0 plane for demonstration.
        """
        interface_points = np.array([[node.X, node.Y, node.Z] for node in self.model_part.GetSubModelPart("BOUNDARY_NODES").Nodes])
        nodes = np.array([[node.X, node.Y, node.Z] for node in self.model_part.Nodes])
        # if interface_points is None:
        #     # Example: interface plane at x = 0
        #     interface_points = np.array([[0.0, y, z] for y, z in zip(nodes[:, 1], nodes[:, 2])])

        # Build KDTree from interface points
        tree = KDTree(interface_points)

        # Compute distances from each node to nearest interface point
        distances, _ = tree.query(nodes)

        # Determine inside/outside using current phi or coordinate sign (x>0 => material)
        signs = np.sign(nodes[:, 0])  # <-- Replace with your material/void criterion
        max_phi = 1/(2*self.k) * np.log((self.k - min_dirac + np.sqrt(self.k*(self.k - min_dirac)))/(2*min_dirac))
        phi_values = signs * distances

        # ---- SCALE AND CLIP ----
        phi_values /= np.max(np.abs(phi_values))  # normalize to [-1, 1]
        phi_values *= max_phi                  # rescale to [-phi_max, phi_max]
        phi_values = np.clip(phi_values, -max_phi, max_phi)

        # Set BC nodes to material
        BC_points = [node.Id for node in self.model_part.GetSubModelPart("DISPLACEMENT_fixed_support").Nodes] + \
                    [node.Id for node in self.model_part.GetSubModelPart("DISPLACEMENT_roller_support").Nodes] + \
                    [node.Id for node in self.model_part.GetSubModelPart("PointLoad3D_load").Nodes]
        for id in BC_points:
            phi_values[id-1] = max_phi

        # Store as nodal expression
        phi_exp = Kratos.Expression.NodalExpression(self.model_part)
        Kratos.Expression.CArrayExpressionIO.Read(phi_exp, phi_values)

        self.control_phi = phi_exp.Clone()

        print(f"Initialized signed distance φ in range: [{phi_values.min():.3f}, {phi_values.max():.3f}]")

    def _ReinitializePhiSussman(self, num_iterations=10, dt=0.01):
        """
        Reinitializes the level set function φ using the Sussman PDE:
            ∂φ/∂τ + S(φ₀)(|∇φ| - 1) = 0
        to maintain φ as a signed distance function.
        
        Args:
            num_iterations (int): Number of pseudo-time iterations
            dt (float): Pseudo-time step (should be small, e.g. 0.01)
        """

        # Clone original φ as φ₀
        phi_0 = self.control_phi.Clone()

        # Precompute sign(phi₀) smoothed: S(phi₀) = phi₀ / sqrt(phi₀² + ε²)
        eps = 1e-6
        phi0_np = phi_0.Evaluate()
        sign_phi = phi0_np / np.sqrt(phi0_np**2 + eps**2)
        sign_phi_exp = Kratos.Expression.NodalExpression(self.model_part)
        Kratos.Expression.CArrayExpressionIO.Read(sign_phi_exp, sign_phi)

        # Pseudo-time marching loop
        for it in range(num_iterations):
            # Compute |∇φ|
            grad_phi_elem = Kratos.Expression.ElementExpression(self.model_part)
            KratosOA.ExpressionUtils.GetGradientExpression(grad_phi_elem, self.control_phi, self.model_part.ProcessInfo[Kratos.DOMAIN_SIZE])
            grad_norm = np.linalg.norm(grad_phi_elem.Evaluate(), axis=1)
            grad_norm_exp = Kratos.Expression.ElementExpression(self.model_part)
            Kratos.Expression.CArrayExpressionIO.Read(grad_norm_exp, grad_norm)

            # Collapse to nodal field for φ update
            grad_norm_nodal = Kratos.Expression.NodalExpression(self.model_part)
            KratosOA.ExpressionUtils.ProjectElementalToNodalViaShapeFunctions(grad_norm_nodal, grad_norm_exp)

            # Compute update: ∂φ/∂τ = -S(φ₀)(|∇φ| - 1)
            update = -sign_phi_exp * (grad_norm_nodal - 1.0)
            # Time integration (explicit Euler)
            self.control_phi += update * dt
            self.control_phi = Kratos.Expression.Utils.Collapse(self.control_phi)

            if it % 5 == 0 and self.echo_level > 0:
                print(f" [Reinit] Iter {it}/{num_iterations} | mean(|∇φ|-1) = {np.mean(np.abs(grad_norm - 1)):.3e}")

        print(f"Sussman reinitialization done. Final φ range: [{min(self.control_phi.Evaluate()):.3f}, {max(self.control_phi.Evaluate()):.3f}]")

    def _GetElementLength(self) -> ContainerExpressionTypes:
        domain_size = []
        size_container = Kratos.Expression.ElementExpression(self.model_part)
        element: Kratos.Element
        for element in self.model_part.Elements:
            geometry: Kratos.Geometry = element.GetGeometry()
            domain_size.append(geometry.DomainSize())
        Kratos.Expression.CArrayExpressionIO.Read(size_container, np.array(domain_size))
        return Kratos.Expression.Utils.Pow(size_container, 1/self.model_part.ProcessInfo[Kratos.DOMAIN_SIZE])
        
    def Check(self) -> None:
        self.filter.Check()
        # pass

    def Finalize(self) -> None:
        self.filter.Finalize()
        # pass

    def GetPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return [Kratos.DENSITY, Kratos.YOUNG_MODULUS]

    def GetEmptyField(self, variable = None) -> ContainerExpressionTypes:
        if (variable in self.GetPhysicalKratosVariables()):
            field = Kratos.Expression.ElementExpression(self.model_part)
            Kratos.Expression.LiteralExpressionIO.SetData(field, 0.0)
            # print("element var", variable)
            return field
        else:
            field = Kratos.Expression.NodalExpression(self.model_part)
            Kratos.Expression.LiteralExpressionIO.SetData(field, 0.0)
            # print("nodal var", variable)
            return field
    
    def GetControlField(self) -> ContainerExpressionTypes:
        return self.control_phi.Clone()

    def MapGradient(self, physical_gradient_variable_container_expression_map: 'dict[SupportedSensitivityFieldVariableTypes, ContainerExpressionTypes]') -> ContainerExpressionTypes:
        # Get physical variables
        keys = physical_gradient_variable_container_expression_map.keys()
        print(physical_gradient_variable_container_expression_map)
        if len(keys) != 2:
            raise RuntimeError(f"Not provided required gradient fields for control \"{self.GetName()}\". Following are the variables:\n\t" + "\n\t".join([k.Name() for k in keys]))
        if  Kratos.DENSITY not in keys or Kratos.YOUNG_MODULUS not in keys:
            raise RuntimeError(f"The required gradient for control \"{self.GetName()}\" w.r.t. DENSITY or YOUNG_MODULUS not found. Followings are the variables:\n\t" + "\n\t".join([k.Name() for k in keys]))

        # first calculate the density partial sensitivity of the response function (elemental)
        d_j_d_density = physical_gradient_variable_container_expression_map[Kratos.DENSITY]

        # second calculate the young modulus partial sensitivity of the response function (elemental)
        d_j_d_youngs = physical_gradient_variable_container_expression_map[Kratos.YOUNG_MODULUS]

        # now compute response function total sensitivity w.r.t. elemental phi
        d_j_d_phi = d_j_d_density * self.d_density_d_phi + d_j_d_youngs * self.d_young_modulus_d_phi

        # get derivative w.r.t. nodal phi
        d_j_d_phi_nodal = Kratos.Expression.NodalExpression(self.model_part)
        KratosOA.ExpressionUtils.ProjectElementalToNodalViaShapeFunctions(d_j_d_phi_nodal, d_j_d_phi)

        return d_j_d_phi_nodal

    def Update(self, control_field: ContainerExpressionTypes) -> bool:
        
        if not IsSameContainerExpression(control_field, self.GetEmptyField()):
             raise RuntimeError(f"Updates for the required element container not found for control \"{self.GetName()}\". [ required model part name: {self.model_part.FullName()}, given model part name: {control_field.GetModelPart().FullName()} ]")

        # Nodal updated control field
        update = control_field - self.control_phi

        if Kratos.Expression.Utils.NormL2(update) > 1e-15:
            # Calculate nodal phi update
            gradient_norm = self.GetEmptyField()
            gradient_norm = Kratos.Expression.Utils.Collapse(self._ComputePhiGradientNorm3(self.control_phi+update))
            nodal_grad = Kratos.Expression.NodalExpression(self.model_part)
            KratosOA.ExpressionUtils.ProjectElementalToNodalViaShapeFunctions(nodal_grad, gradient_norm)
            # compute timestep from CFL
            time_step = self._ComputeCFLCondition(update)
            print(f"Time_step: {time_step}")

            new_update = update * nodal_grad * time_step
            # self.filter.ForwardFilterField(new_update)
            self.control_phi += new_update
            self.control_phi = Kratos.Expression.Utils.Collapse(self.control_phi)

            # Reinitialize the LSF as Signed Distance Function
            # if self.optimization_problem.GetStep() % self.reinit_interval == 0:
            #     self._ReinitializePhiSussman(num_iterations=15, dt = Kratos.Expression.Utils.EntityMin(self.element_len).Evaluate()[0] * time_step)

            self._UpdateAndOutputFields(new_update)
            self.filter.Update()
            return True

        return False

    def _UpdateAndOutputFields(self, update: ContainerExpressionTypes) -> None:
        # compute elemental phi
        elemental_phi = Kratos.Expression.ElementExpression(self.model_part)
        KratosOA.ExpressionUtils.ProjectNodalToElementalViaShapeFunctions(elemental_phi, self.control_phi)
        
        # get elemental heaviside values
        heaviside = Kratos.Expression.ElementExpression(self.model_part)
        KratosOA.ExpressionUtils.Heaviside(heaviside, elemental_phi, 100, self.k)
        
        # update physical variable fields: density = H(phi)*density_0
        density = Kratos.Expression.ElementExpression(self.model_part)
        density = heaviside * self.density
        KratosOA.PropertiesVariableExpressionIO.Write(density, Kratos.DENSITY)
        if self.consider_recursive_property_update:
            KratosOA.OptimizationUtils.UpdatePropertiesVariableWithRootValueRecursively(density.GetContainer(), Kratos.DENSITY)
        self.un_buffered_data.SetValue("DENSITY", density.Clone(), overwrite=True)

        # E = H(phi)*E_0
        youngs_modulus = Kratos.Expression.ElementExpression(self.model_part)
        youngs_modulus = heaviside * (self.young_modulus - self.epsilon) + self.epsilon
        KratosOA.PropertiesVariableExpressionIO.Write(youngs_modulus, Kratos.YOUNG_MODULUS)
        if self.consider_recursive_property_update:
            KratosOA.OptimizationUtils.UpdatePropertiesVariableWithRootValueRecursively(youngs_modulus.GetContainer(), Kratos.YOUNG_MODULUS)
        self.un_buffered_data.SetValue("YOUNG_MODULUS", youngs_modulus.Clone(), overwrite=True)

        # Calculate elemental Dirac Delta of phi (dH/d_phi)
        dirac_delta = Kratos.Expression.ElementExpression(self.model_part)
        KratosOA.ExpressionUtils.DiracDelta(dirac_delta, elemental_phi, 100, self.k)

        # now calculate the total sensitivities of density and E w.r.t. phi (elemental)
        self.d_density_d_phi = dirac_delta * self.density
        self.d_young_modulus_d_phi = dirac_delta * self.young_modulus

        # now output the fields
        un_buffered_data = ComponentDataView(self, self.optimization_problem).GetUnBufferedData()
        un_buffered_data.SetValue(f"{self.GetName()}", self.control_phi.Clone(),overwrite=True)
        if self.output_all_fields:
            un_buffered_data.SetValue(f"{self.GetName()}_update", update.Clone(), overwrite=True)
            un_buffered_data.SetValue(f"dDENSITY_d{self.GetName()}", self.d_density_d_phi.Clone(), overwrite=True)
            un_buffered_data.SetValue(f"dYOUNG_MODULUS_d{self.GetName()}", self.d_young_modulus_d_phi.Clone(), overwrite=True)

        if True in (np.isnan(np.array(heaviside.Evaluate())) | np.isinf(np.array(heaviside.Evaluate()))):
            np.savetxt("density.txt", density.Evaluate())
            np.savetxt("youngs.txt", youngs_modulus.Evaluate())
            np.savetxt("heaviside.txt", heaviside.Evaluate())
            np.savetxt("d_young_modulus_d_phi.txt", self.d_young_modulus_d_phi.Evaluate())
            np.savetxt("elemental_phi.txt", elemental_phi.Evaluate())
            raise RuntimeError(f"Nan in dirac")

    def _ComputePhiGradientNorm3(self, delta_phi: ContainerExpressionTypes) -> ContainerExpressionTypes:
        # delta_phi has to be nodal expression
        gradient_elemental = Kratos.Expression.ElementExpression(self.model_part)
        KratosOA.ExpressionUtils.GetGradientExpression(gradient_elemental, delta_phi, self.model_part.ProcessInfo[Kratos.DOMAIN_SIZE])
        # np.savetxt("grad_direct.txt", gradient_elemental.Evaluate())
        gradient_norm = np.linalg.norm(gradient_elemental.Evaluate(), axis=1)
        # np.savetxt("grad_norm.txt", gradient_norm)
        # convert nupmy array back to elemental expression
        norm_exp = Kratos.Expression.ElementExpression(self.model_part)
        Kratos.Expression.CArrayExpressionIO.Read(norm_exp, gradient_norm)
        return norm_exp

    def _ComputeCFLCondition(self, design_velocity: ContainerExpressionTypes) -> float:
        elemental_velocity = Kratos.Expression.ElementExpression(self.model_part)
        KratosOA.ExpressionUtils.ProjectNodalToElementalViaShapeFunctions(elemental_velocity, design_velocity)
        time = self.element_len / Kratos.Expression.Utils.Abs(elemental_velocity)
        return min(time.Evaluate()) * 0.1    

    @staticmethod
    def ComputeHeavisideValue(phi: ContainerExpressionTypes, k: float) -> ContainerExpressionTypes:  
        e = phi.Clone()
        Kratos.Expression.LiteralExpressionIO.SetData(e, np.exp(1))
        return Kratos.Expression.Utils.Pow(Kratos.Expression.Utils.Pow(e, phi * k * (-2)) + 1, -1)
    @staticmethod
    def ComputeHevisideGradientValue(phi: ContainerExpressionTypes, k: float) -> ContainerExpressionTypes:
        e = phi.Clone()
        Kratos.Expression.LiteralExpressionIO.SetData(e, np.exp(1))
        return Kratos.Expression.Utils.Pow(Kratos.Expression.Utils.Pow(e, phi * k * 2) + 1, -2) * Kratos.Expression.Utils.Pow(e, phi * k * 2) * k * 2

    
if __name__ == "__main__":
    import numpy
    # Prep for data
    model = Kratos.Model()
    model_part = model.CreateModelPart("test")
    k = 0.1
    for i in range(10):
        model_part.CreateNewNode(i + 1, 0, 0, 0).SetValue(Kratos.DENSITY, -2 + i / 2.5)

    phi = Kratos.Expression.NodalExpression(model_part)
    Kratos.Expression.VariableExpressionIO.Read(phi, Kratos.DENSITY, is_historical=False)

    # Calculate heaviside with fuction
    haviside_phi = LevelSetControlHJ.ComputeHeavisideValue(phi, k)

    print(f"Heaviside function: \n{haviside_phi.Evaluate()}")

    # manual calculation of heaviside
    numpy_phi = phi.Evaluate()
    numpy_heaviside = 1 / (1 + numpy.exp(numpy_phi * k * (-2)))

    print(f"Numpy Heaviside: \n{numpy_heaviside}")

    # Calculate Heaviside gradient with function
    heaviside_gradient_phi = LevelSetControlHJ.ComputeHevisideGradientValue(phi, k)

    print(f"Heaviside gradient function: \n{heaviside_gradient_phi.Evaluate()}")

    # Calculate Heaviside gradient manually with finite differences
    perturbation = 1e-5
    FD_heaviside_gradient = (LevelSetControlHJ.ComputeHeavisideValue(phi + perturbation, k) - LevelSetControlHJ.ComputeHeavisideValue(phi, k))/perturbation

    print(f"Heaviside gradient FD: \n{FD_heaviside_gradient.Evaluate()}")

    # Calculate Heaviside gradient with numpy
    numpy_heaviside_gradient = numpy.exp(numpy_phi * k * 2) * k * 2 / (1 + numpy.exp(numpy_phi * k * 2))**2

    print(f"Heaviside gradient numpy: \n{numpy_heaviside_gradient}")