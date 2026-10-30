import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.kratos_utilities as KratosUtilities
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid

have_convection_diffusion = KratosUtilities.CheckIfApplicationsAvailable("ConvectionDiffusionApplication")
if have_convection_diffusion:
    from KratosMultiphysics.ConvectionDiffusionApplication import coupled_fluid_thermal_solver

@UnitTest.skipIfApplicationsNotAvailable("ConvectionDiffusionApplication")
class CoupledFluidThermalSolverTest(UnitTest.TestCase):

    def testOutletInflowContribution(self):
        # Unit square with two triangles and its four edges as conditions (outlet on the right edge, x = 1)
        fluid_density = 1.2
        thermal_density = 0.0
        settings = KratosMultiphysics.Parameters("""{
            "solver_type" : "ThermallyCoupled",
            "domain_size" : 2,
            "echo_level" : 0,
            "fluid_solver_settings" : {
                "solver_type" : "Monolithic",
                "model_part_name" : "FluidModelPart",
                "domain_size" : 2,
                "model_import_settings" : {"input_type" : "use_input_model_part"},
                "material_import_settings" : {"materials_filename" : ""},
                "formulation" : {"element_type" : "qsvms"},
                "volume_model_part_name" : "Fluid",
                "skin_parts" : ["Walls", "Outlet"],
                "no_skin_parts" : [],
                "time_stepping" : {"automatic_time_step" : false, "time_step" : 0.1}
            },
            "thermal_solver_settings" : {
                "solver_type" : "transient",
                "analysis_type" : "linear",
                "model_part_name" : "ThermalModelPart",
                "domain_size" : 2,
                "model_import_settings" : {"input_type" : "use_input_model_part"},
                "material_import_settings" : {"materials_filename" : ""},
                "problem_domain_sub_model_part_list" : ["Fluid"],
                "processes_sub_model_part_list" : ["Walls", "Outlet"],
                "time_stepping" : {"time_step" : 0.1}
            }
        }""")
        model = KratosMultiphysics.Model()
        fluid_model_part = model.CreateModelPart("FluidModelPart")
        solver = coupled_fluid_thermal_solver.CreateSolver(model, settings)
        solver.AddVariables()

        fluid_properties = fluid_model_part.CreateNewProperties(1)
        fluid_properties.SetValue(KratosMultiphysics.DENSITY, fluid_density)
        fluid_properties.SetValue(KratosMultiphysics.DYNAMIC_VISCOSITY, 1.0e-3)
        fluid_properties.SetValue(KratosMultiphysics.CONSTITUTIVE_LAW, KratosFluid.Newtonian2DLaw())
        fluid_model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        fluid_model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        fluid_model_part.CreateNewNode(3, 1.0, 1.0, 0.0)
        fluid_model_part.CreateNewNode(4, 0.0, 1.0, 0.0)
        fluid_sub_model_part = fluid_model_part.CreateSubModelPart("Fluid")
        fluid_sub_model_part.AddNodes([1, 2, 3, 4])
        fluid_sub_model_part.CreateNewElement("Element2D3N", 1, [1, 2, 3], fluid_properties)
        fluid_sub_model_part.CreateNewElement("Element2D3N", 2, [1, 3, 4], fluid_properties)
        walls_model_part = fluid_model_part.CreateSubModelPart("Walls")
        walls_model_part.AddNodes([1, 2, 3, 4])
        walls_model_part.CreateNewCondition("LineCondition2D2N", 1, [1, 2], fluid_properties)
        walls_model_part.CreateNewCondition("LineCondition2D2N", 2, [3, 4], fluid_properties)
        walls_model_part.CreateNewCondition("LineCondition2D2N", 3, [4, 1], fluid_properties)
        outlet_model_part = fluid_model_part.CreateSubModelPart("Outlet")
        outlet_model_part.AddNodes([2, 3])
        outlet_model_part.CreateNewCondition("LineCondition2D2N", 4, [2, 3], fluid_properties)

        # The thermal elements get their own properties (e.g. from a thermal materials file), without the fluid DENSITY
        solver.ImportModelPart()
        thermal_model_part = model.GetModelPart("ThermalModelPart")
        thermal_properties = thermal_model_part.CreateNewProperties(0)
        thermal_properties.SetValue(KratosMultiphysics.DENSITY, thermal_density)
        for element in thermal_model_part.Elements:
            element.Properties = thermal_properties
        solver.PrepareModelPart()
        solver.AddDofs()
        solver.Initialize()

        # Inflow through the outlet with the outlet inflow contribution activated
        inflow_velocity = -2.0
        for node in fluid_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, [inflow_velocity, 0.0, 0.0])
        process_info = fluid_model_part.ProcessInfo
        process_info.SetValue(KratosFluid.CHARACTERISTIC_VELOCITY, 1.0)
        outlet_condition = fluid_model_part.GetCondition(4)
        outlet_condition.Set(KratosMultiphysics.OUTLET, True)

        # Outlet inflow contribution: difference between the RHS with and without it
        rhs_with, rhs_without = KratosMultiphysics.Vector(), KratosMultiphysics.Vector()
        process_info.SetValue(KratosFluid.OUTLET_INFLOW_CONTRIBUTION_SWITCH, True)
        outlet_condition.CalculateRightHandSide(rhs_with, process_info)
        process_info.SetValue(KratosFluid.OUTLET_INFLOW_CONTRIBUTION_SWITCH, False)
        outlet_condition.CalculateRightHandSide(rhs_without, process_info)

        # The contribution 0.5 * rho * |v|^2 * S_0 * n (S_0 = 1 for this inflow) integrated over the unit length outlet
        # must use the fluid density (zero if the conditions took the thermal elements as parents)
        expected_force_x = 0.5 * fluid_density * inflow_velocity**2
        self.assertAlmostEqual((rhs_with[0] - rhs_without[0]) + (rhs_with[3] - rhs_without[3]), expected_force_x, 12)
        self.assertAlmostEqual((rhs_with[1] - rhs_without[1]) + (rhs_with[4] - rhs_without[4]), 0.0, 12)

if __name__ == '__main__':
    UnitTest.main()
