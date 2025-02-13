import KratosMultiphysics as KM
import KratosMultiphysics.GeoMechanicsApplication as KGM
import KratosMultiphysics.StructuralMechanicsApplication as KSM
import matplotlib.pyplot as plt
import KratosMultiphysics.GeoMechanicsApplication.geomechanics_analysis as analysis

class element_test():

    def __init__(self):
        self.model = KM.Model()
        self.model_part = self.model.CreateModelPart("test")
        self._add_variables()
        self.model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 2)

        self._set_constitutive_law()
        self._set_material_properties()
        self._create_geometry()
        self._addDoF()
        self._applyBC()
        self._applyStress()
        self._setupSolver()
        print("Initialized")

    def _create_geometry(self):
        self.node_1 = self.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        self.node_2 = self.model_part.CreateNewNode(2, 0.0, 1.0, 0.0)
        self.node_3 = self.model_part.CreateNewNode(3, 1.0, 1.0, 0.0)
        self.node_4 = self.model_part.CreateNewNode(4, 1.0, 0.0, 0.0)
        self.nnodes = 4
        self.element = self.model_part.CreateNewElement('UPwSmallStrainAxisymmetricElement2D4N', 1, [4,3,2,1], self.model_part.GetProperties()[1])
        self.baseNodes = self.model_part.CreateSubModelPart("BaseNodes")
        self.baseNodes.AddNodes([1])
        self.axisNodes = self.model_part.CreateSubModelPart("AxisNodes")
        self.axisNodes.AddNodes([1, 2])
        self.topNodes = self.model_part.CreateSubModelPart("TopNodes")
        self.topNodes.AddNodes([2, 3])
        self.bottomRightNodes = self.model_part.CreateSubModelPart("BottomRightNodes")
        self.bottomRightNodes.AddNodes([4])
        self.rightNodes = self.model_part.CreateSubModelPart("RightNodes")
        self.rightNodes.AddNodes([3, 4])

    def _set_constitutive_law(self):
        self.cl = KSM.LinearElasticPlaneStrain2DLaw()

    def _set_material_properties(self):
        print("Setting material properties")
        element_properties = self.model_part.GetProperties()[1]
        element_properties.SetValue(KGM.IGNORE_UNDRAINED, False)
        element_properties.SetValue(KGM.YOUNG_MODULUS, 2e9)
        element_properties.SetValue(KGM.POISSON_RATIO, 0.0)
        element_properties.SetValue(KGM.DENSITY_SOLID, 2650)
        element_properties.SetValue(KGM.DENSITY_WATER, 1000)
        element_properties.SetValue(KGM.POROSITY, 0.3)
        element_properties.SetValue(KGM.BULK_MODULUS_SOLID, 1.0e9)
        element_properties.SetValue(KGM.BULK_MODULUS_FLUID, 2.0e6)
        element_properties.SetValue(KGM.PERMEABILITY_XX, 4.5e-13)
        element_properties.SetValue(KGM.PERMEABILITY_YY, 4.5e-13)
        element_properties.SetValue(KGM.PERMEABILITY_XY, 0.0)
        element_properties.SetValue(KGM.DYNAMIC_VISCOSITY, 8.90e-7)
        element_properties.SetValue(KGM.THICKNESS, 1.0)
        element_properties.SetValue(KGM.BIOT_COEFFICIENT, 1.0)
        element_properties.SetValue(KGM.RETENTION_LAW, "SaturatedLaw")
        element_properties.SetValue(KGM.SATURATED_SATURATION, 1.0)
        element_properties.SetValue(KGM.CONSTITUTIVE_LAW, self.cl)
        g = [0, 0, 0]
        element_properties.SetValue(KM.VOLUME_ACCELERATION, g)
        print("Setting material properties - Completed")

    def _add_variables(self):
        self.model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        self.model_part.AddNodalSolutionStepVariable(KM.REACTION)
        self.model_part.AddNodalSolutionStepVariable(KM.VOLUME_ACCELERATION)
        self.model_part.AddNodalSolutionStepVariable(KM.WATER_PRESSURE)
        self.model_part.AddNodalSolutionStepVariable(KGM.DT_WATER_PRESSURE)
        self.model_part.AddNodalSolutionStepVariable(KM.NORMAL_CONTACT_STRESS)
        self.model_part.AddNodalSolutionStepVariable(KM.REACTION_WATER_PRESSURE)

    def _addDoF(self):
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_X, KM.REACTION_X, self.model_part)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y, self.model_part)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z, self.model_part)
        KM.VariableUtils().AddDof(KM.WATER_PRESSURE, KM.REACTION_WATER_PRESSURE, self.model_part)
        KM.VariableUtils().AddDof(KGM.DT_WATER_PRESSURE, self.model_part)

    def _applyBC(self):
        # completely fixed
        KM.VariableUtils().SetVariable(KM.DISPLACEMENT_X, 0.0, self.baseNodes.Nodes)
        KM.VariableUtils().SetVariable(KM.DISPLACEMENT_Y, 0.0, self.baseNodes.Nodes)
        KM.VariableUtils().ApplyFixity(KM.DISPLACEMENT_X, True, self.baseNodes.Nodes)
        KM.VariableUtils().ApplyFixity(KM.DISPLACEMENT_Y, True, self.baseNodes.Nodes)

        # slide fixity Y
        KM.VariableUtils().SetVariable(KM.DISPLACEMENT_X, 0.0, self.axisNodes.Nodes)
        KM.VariableUtils().ApplyFixity(KM.DISPLACEMENT_X, True, self.axisNodes.Nodes)

        # slide fixity X
        KM.VariableUtils().SetVariable(KM.DISPLACEMENT_Y, 0.0, self.bottomRightNodes.Nodes)
        KM.VariableUtils().ApplyFixity(KM.DISPLACEMENT_Y, True, self.bottomRightNodes.Nodes)

    def _applyStress(self):
        KM.VariableUtils().SetVariable(KSM.POINT_LOAD_X, -1.0, self.rightNodes.Nodes)

    def _applyStrain(self, strain):
        KM.VariableUtils().SetVariable(KM.DISPLACEMENT_Y, strain, self.topNodes.Nodes)
        KM.VariableUtils().ApplyFixity(KM.DISPLACEMENT_Y, True, self.topNodes.Nodes)

    def _setupSolver(self):





        # self.model_part.SetBufferSize(2)
        #
        # time_scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
        # linear_solver = KM.SkylineLUFactorizationSolver()
        # relative_tolerance = 1e-7
        # absolute_tolerance = 1e-7
        # conv_criteria = KM.ResidualCriteria(relative_tolerance, absolute_tolerance)
        # conv_criteria.SetEchoLevel(0)
        # maximum_iterations = 100
        # compute_reactions = True
        # reform_dofs_at_each_iteration = True
        # move_mesh_flag = False
        # self.solver = KM.ResidualBasedNewtonRaphsonStrategy(
        #     self.model_part,
        #     time_scheme,
        #     linear_solver,
        #     conv_criteria,
        #     maximum_iterations,
        #     compute_reactions,
        #     reform_dofs_at_each_iteration,
        #     move_mesh_flag
        # )

    def solve(self, strain):
        self.model_part.CloneTimeStep(1)
        self._applyStrain(strain)
        self.solver.Solve()

    def lab_experiment(self, steps, max_strain):
        self.stress_out = []
        self.strain_out = []
        for i in range(steps):
            strain = i*max_strain/steps
            self.solve(strain)
            self.strain_out.append(self.node_2.GetSolutionStepValue(KM.DISPLACEMENT_Y))
            self.stress_out.append(self.element.CalculateOnIntegrationPoints(KM.CAUCHY_STRESS_VECTOR, self.model_part.ProcessInfo)[0][1])

    def plot_results(self):
        print(len(self.strain_out), self.strain_out[-1])
        print(len(self.stress_out), self.stress_out[-1])
        plt.plot(self.strain_out, self.stress_out, label='Stress-Strain')
        plt.show()


if __name__ == "__main__":
    test = element_test()
    test.lab_experiment(100, -0.1)
    test.plot_results()
