import numpy as np
from pathlib import Path
import KratosMultiphysics
import KratosMultiphysics.RomApplication as KratosROM

def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object.")
    return RomResidualsOutputProcess(model, settings["Parameters"])

class RomResidualsOutputProcess(KratosMultiphysics.OutputProcess):
    """Outputs HROM projected residuals directly using the RomResidualsUtility."""

    def __init__(self, model, settings):
        super().__init__()
        settings.ValidateAndAssignDefaults(self.GetDefaultParameters())
        self.settings = settings
        self.model_part = model[settings["model_part_name"].GetString()]

        self.output_path = Path(self.settings["output_path"].GetString())
        self.output_path.mkdir(parents=True, exist_ok=True)

        controller_settings = KratosMultiphysics.Parameters("""{}""")
        controller_settings.AddString("model_part_name", self.model_part.FullName())
        controller_settings.AddValue("output_control_type", settings["output_control_type"])
        controller_settings.AddValue("output_interval", settings["output_interval"])
        self.controller = KratosMultiphysics.OutputController(self.model_part.GetModel(), controller_settings)

        self.solver = None

    @classmethod
    def GetDefaultParameters(cls):
        return KratosMultiphysics.Parameters("""{
            "model_part_name": "",
            "output_control_type": "step",
            "output_interval": 1,
            "output_path": "rom_data/Residuals",
            "sub_solver_name": "",
            "rom_basis_output_name": "RomParameters",
            "rom_basis_output_folder": "rom_data",
            "range_of_elements_to_export": [0, -1]
        }""")

    def SetDependencies(self, main_solver, rom_parameters):
        """Inject the Python solver object and extract required projection settings."""
        sub_solver_name = self.settings["sub_solver_name"].GetString()
        self.solver = getattr(main_solver, sub_solver_name) if sub_solver_name else main_solver

        self.rom_settings = rom_parameters["rom_settings"]
        self.projection_strategy = rom_parameters["projection_strategy"].GetString()
        self.num_of_right_rom_dofs = self.rom_settings["number_of_rom_dofs"].GetInt()
        self.rom_settings.RemoveValue("rom_bns_settings")


    def _GetJacobianPhiMultiplication(self, computing_model_part):
        """Assembles the Jacobian and multiplies it by the ROM basis for LSPG."""
        jacobian_matrix = KratosMultiphysics.CompressedMatrix()
        builder_and_solver = self.solver._GetBuilderAndSolver()
        system_size = builder_and_solver.GetEquationSystemSize()

        residual_vector = KratosMultiphysics.Vector(system_size)
        delta_x_vector = KratosMultiphysics.Vector(system_size)

        builder_and_solver.BuildAndApplyDirichletConditions(self.solver._GetScheme(), computing_model_part, jacobian_matrix, residual_vector, delta_x_vector)

        right_rom_basis = KratosMultiphysics.Matrix(system_size, self.num_of_right_rom_dofs)
        builder_and_solver.GetRightROMBasis(computing_model_part, right_rom_basis)

        jacobian_scipy_format = KratosMultiphysics.scipy_conversion_tools.to_csr(jacobian_matrix)
        return jacobian_scipy_format @ right_rom_basis

    def _GetCurrentResidualsProjected(self):
        """Calculates the projected residuals using the C++ RomResidualsUtility."""
        computing_model_part = self.solver.GetComputingModelPart()

        if not hasattr(self, '__rom_residuals_utility'):
            self.__rom_residuals_utility = KratosROM.RomResidualsUtility(
                computing_model_part,
                self.rom_settings,
                self.solver._GetScheme())

        if self.projection_strategy == "galerkin":
            res_mat = self.__rom_residuals_utility.GetProjectedResidualsOntoPhi()
        elif self.projection_strategy == "lspg":
            jacobian_phi_product = self._GetJacobianPhiMultiplication(computing_model_part)
            res_mat = self.__rom_residuals_utility.GetProjectedResidualsOntoJPhi(jacobian_phi_product)
        elif self.projection_strategy == "petrov_galerkin":
            res_mat = self.__rom_residuals_utility.GetProjectedResidualsOntoPsi()
        else:
            raise Exception(f"Projection strategy '{self.projection_strategy}' is not supported.")

        return np.asarray(res_mat)

    def IsOutputStep(self):
        return self.controller.Evaluate()


    def CaptureResiduals(self):
        """Called manually by RomAnalysis before the solver clears the database."""
        # Only perform the heavy matrix math if we are actually going to output this step
        if self.IsOutputStep():
            self._buffered_res_mat = self._GetCurrentResidualsProjected()


    def PrintOutput(self):
        if self.solver is None:
            return

        if hasattr(self, "_buffered_res_mat") and self._buffered_res_mat is not None:
            # 1. Slice based on JSON settings
            export_range = self.settings["range_of_elements_to_export"].GetVector()
            start = int(export_range[0])
            end = int(export_range[1]) if int(export_range[1]) >= 0 else None
            sliced_res_mat = self._buffered_res_mat[start:end, :]

            # 2. Label formatting
            if self.settings["output_control_type"].GetString() == "time":
                time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
                file_label = f"{time:.7f}"
            else:
                step = self.model_part.ProcessInfo[KratosMultiphysics.STEP]
                file_label = f"{step}"

            # 3. Perform I/O
            file_path = self.output_path / f"Residual_{file_label}.npy"
            np.save(file_path, sliced_res_mat)

            # 4. Clear the memory buffer immediately
            self._buffered_res_mat = None

        self.controller.Update()
