# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.IgaApplication as KratosIGA

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_coupling_operation import CoSimulationCouplingOperation

# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools


def Create(*args):
    return ComputeInterfaceTractionShell3pIGA(*args)


class ComputeInterfaceTractionShell3pIGA(CoSimulationCouplingOperation):
    """Computes and stores interface tractions on Shell 3p IGA interface conditions."""

    def __init__(self, settings, solver_wrappers, process_info, data_communicator):
        super().__init__(settings, process_info, data_communicator)

        model = solver_wrappers[self.settings["solver"].GetString()].model

        self.model_part_name = self.settings["model_part_name"].GetString()
        self.model_part = model[self.model_part_name]

        self.interval = KM.IntervalUtility(self.settings)

        self.traction_variable = KM.KratosGlobals.GetVariable(
            self.settings["traction_variable_name"].GetString()
        )

    def InitializeCouplingIteration(self):
        if self.interval.IsInInterval(self.model_part.ProcessInfo[KM.TIME]):
            self._ComputeInterfaceTractionShell3pIGA()

            if self.echo_level > 0:
                cs_tools.cs_print_info(
                    self._ClassName(),
                    "Interface tractions calculated in ModelPart: " + self.model_part_name
                )
    
    def _CopyMissingMaterialPropertiesToConditions(self):
        root_model_part = self.model_part.GetRootModelPart()

        elem_props = root_model_part.Elements[1].Properties

        vars_to_copy = [
            KM.DENSITY,
            KM.YOUNG_MODULUS,
            KM.POISSON_RATIO,
            KM.THICKNESS,
            KM.CONSTITUTIVE_LAW
        ]

        for cond in self.model_part.Conditions:
            cond_props = cond.Properties

            for var in vars_to_copy:
                if elem_props.Has(var) and not cond_props.Has(var):
                    cond_props.SetValue(var, elem_props[var])


    def _ComputeInterfaceTractionShell3pIGA(self):
        self._CopyMissingMaterialPropertiesToConditions()
            
        KratosIGA.ComputeInterfaceTractionShell3pUtility.ComputeAndSetInterfaceTraction(
            self.model_part,
            self.traction_variable
        )

    @classmethod
    def _GetDefaultParameters(cls):
        this_defaults = KM.Parameters("""{
            "solver"                 : "UNSPECIFIED",
            "model_part_name"        : "",
            "traction_variable_name" : "INTERFACE_TRACTION",
            "interval"               : [0.0, 1e30]
        }""")
        this_defaults.AddMissingParameters(super()._GetDefaultParameters())
        return this_defaults