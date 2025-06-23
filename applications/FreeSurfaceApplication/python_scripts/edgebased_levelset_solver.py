# Kratos imports
import KratosMultiphysics
from KratosMultiphysics.python_solver import PythonSolver
import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory

import KratosMultiphysics.FreeSurfaceApplication as FreeSurface

# STL imports
import enum
import pathlib


class EdgeBasedLevelSetSolver(PythonSolver):

    class PorousResistanceComputation(enum.IntEnum):
        NONE = 0
        ERGUN = 1
        CUSTOM = 2

    def __init__(self,
                 model: KratosMultiphysics.Model,
                 parameters: KratosMultiphysics.Parameters):
        """
        Default parameters.
        {
            "model_part_name"                       : "",
            "domain_size"                           : 3,
            "density"                               : 0.0,
            "viscosity"                             : 0.0,
            "wall_law_y"                            : 0.0,
            "stabdt_pressure_factor"                : 1.0,
            "stabdt_convection_factor"              : 0.1,
            "use_mass_correction"                   : true,
            "redistance_frequency"                  : 5,
            "extrapolation_layers"                  : 5,
            "tau2_factor"                           : 0.0,
            "max_safety_factor"                     : 0.5,
            "max_time_step_size"                    : 1e-2,
            "initial_time_step_size"                : 1e-5,
            "number_of_initial_time_steps"          : 10,
            "reduction_on_failure"                  : 0.3,
            "assume_constant_pressure"              : true,
            "use_parallel_distance_calculation"     : false,
            "compute_porous_resistance_law"         : "NONE",    ["NONE", "ERGUN", "CUSTOM"]
            "echo_level"                            : 0,
            "solver_type"                           : "EdgebasedLevelset",
            "linear_solver_settings"      : {
                "solver_type"                       : "amgcl"
            },
            "model_import_settings"                 : { }
        }
        """
        # Base class defines:
        # - self.model
        # - self.settings
        PythonSolver.__init__(self, model, parameters)

        # Get root ModelPart
        model_part_name = self.settings["model_part_name"].GetString()
        if self.model.HasModelPart(model_part_name):
            self.model_part = self.model.GetModelPart(model_part_name)
        else:
            self.model_part = self.model.CreateModelPart(model_part_name)

        # Parse numeric parameters
        self.domain_size = self.settings["domain_size"].GetInt()
        self.density = self.settings["density"].GetDouble()
        self.viscosity = self.settings["viscosity"].GetDouble()
        self.wall_law_y = self.settings["wall_law_y"].GetDouble()
        self.stabdt_pressure_factor = self.settings["stabdt_pressure_factor"].GetDouble()
        self.stabdt_convection_factor = self.settings["stabdt_convection_factor"].GetDouble()
        self.redistance_frequency = self.settings["redistance_frequency"].GetInt()
        self.extrapolation_layers = self.settings["extrapolation_layers"].GetInt()
        self.tau2_factor = self.settings["tau2_factor"].GetDouble()
        self.max_safety_factor = self.settings["max_safety_factor"].GetDouble()
        self.max_time_step_size = self.settings["max_time_step_size"].GetDouble()
        self.initial_time_step_size = self.settings["initial_time_step_size"].GetDouble()
        self.number_of_initial_time_steps = self.settings["number_of_initial_time_steps"].GetInt()
        self.reduction_on_failure = self.settings["reduction_on_failure"].GetDouble()

        # Parse options
        self.use_mass_correction = self.settings["use_mass_correction"].GetBool()
        self.assume_constant_pressure = self.settings["assume_constant_pressure"].GetBool()
        self.use_parallel_distance_calculation = self.settings["use_parallel_distance_calculation"].GetBool()
        self.compute_porous_resistance_law = EdgeBasedLevelSetSolver.PorousResistanceComputation[self.settings["compute_porous_resistance_law"].GetString()]
        self.pressure_linear_solver = self._CreateLinearSolver()

        # Preparations before the model part is read
        # (apart from adding variables)
        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, self.domain_size)

        self.distance_size = 0
        self.matrix_container = self.__MakeMatrixContainer()
        self.distance_utils = self.__MakeDistanceUtilities()

        # Declare other members
        self.fluid_solver = None # initialized by __MakeEdgeBasedLevelSet after reading the ModelPart
        self.safety_factor = self.max_safety_factor
        self.current_step_size = self.initial_time_step_size
        self.current_max_step_size = self.max_time_step_size
        self.timer = KratosMultiphysics.Timer()

    @staticmethod
    def GetDefaultParameters() -> KratosMultiphysics.Parameters:
        return KratosMultiphysics.Parameters("""{
            "model_part_name"                       : "",
            "domain_size"                           : 3,
            "density"                               : 0.0,
            "viscosity"                             : 0.0,
            "wall_law_y"                            : 0.0,
            "stabdt_pressure_factor"                : 1.0,
            "stabdt_convection_factor"              : 0.1,
            "use_mass_correction"                   : true,
            "redistance_frequency"                  : 5,
            "extrapolation_layers"                  : 5,
            "tau2_factor"                           : 0.0,
            "max_safety_factor"                     : 0.5,
            "max_time_step_size"                    : 1e-2,
            "initial_time_step_size"                : 1e-5,
            "number_of_initial_time_steps"          : 10,
            "reduction_on_failure"                  : 0.3,
            "assume_constant_pressure"              : true,
            "use_parallel_distance_calculation"     : false,
            "compute_porous_resistance_law"         : "NONE",
            "echo_level"                            : 0,
            "solver_type"                           : "EdgebasedLevelset",
            "linear_solver_settings"      : {
                "solver_type"                       : "amgcl"
            },
            "model_import_settings"                 : { }
        }""")

    @staticmethod
    def GetMinimumBufferSize() -> int:
        return 2

    def AddVariables(self) -> None:
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.EXTERNAL_PRESSURE)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.AUX_INDEX)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.BODY_FORCE)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESS_PROJ)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.POROSITY)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.LIN_DARCY_COEF)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NONLIN_DARCY_COEF)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DIAMETER)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.STRUCTURE_VELOCITY)

    def AddDofs(self) -> None:
        variable_utils = KratosMultiphysics.VariableUtils()
        variable_utils.AddDof(KratosMultiphysics.PRESSURE, self.model_part)
        variable_utils.AddDof(KratosMultiphysics.VELOCITY_X, self.model_part)
        variable_utils.AddDof(KratosMultiphysics.VELOCITY_Y, self.model_part)
        variable_utils.AddDof(KratosMultiphysics.VELOCITY_Z, self.model_part)

    def ImportModelPart(self) -> None:
        self._ImportModelPart(self.model_part, self.settings["model_import_settings"].Clone())

    def PrepareModelPart(self) -> None:
        self.model_part.SetBufferSize(self.GetMinimumBufferSize())

        # Set ProcessInfo variables
        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0)
        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.STEP, 0)

    def Check(self) -> None:
        super().Check()

    def Initialize(self) -> None:
        # Get rid of isolated nodes
        KratosMultiphysics.FindGlobalNodalNeighboursProcess(self.model_part).Execute()
        KratosMultiphysics.FindGlobalNodalElementalNeighboursProcess(self.model_part).Execute()
        KratosMultiphysics.EliminateIsolatedNodesProcess(self.model_part).Execute()

        # Build edge data structure
        self.matrix_container.ConstructCSRVector(self.model_part)
        self.matrix_container.BuildCSRData(self.model_part)

        # Initialize solver
        self.fluid_solver = self.__MakeEdgeBasedLevelSet()

        if self.use_parallel_distance_calculation:
            self.distance_size = 3.0 * self.distance_utils.FindMaximumEdgeSize()
        else:
            self.distance_size = 3.0 * self.distance_utils.FindMaximumEdgeSize(self.model_part)
        self.fluid_solver.SetShockCapturingCoefficient(0.0)

        # Note:
        # original script prints the number of edges with positive/negative DISTANCE values
        # before and after calling self.fluid_solver.Initialize
        self.fluid_solver.Initialize()

        if 1e-10 < self.wall_law_y:
            self.fluid_solver.ActivateWallResistance(self.wall_law_y)

    def InitializeSolutionStep(self) -> None:
        # Transfer density to the nodes
        for node in self.model_part.Nodes:
            if node.GetSolutionStepValue(KratosMultiphysics.DISTANCE) <= 0.0:
                node.SetSolutionStepValue(KratosMultiphysics.DENSITY, self.density)
            else:
                node.SetSolutionStepValue(KratosMultiphysics.DENSITY, 0.0)

    def AdvanceInTime(self, current_time: float) -> float:
        """
        Compute a bounded variable time step and create new step data.
        Note: this function modifies
            -> self.safety_factor
            -> self.current_max_step_size
            -> self.current_step_size
        """
        # Bound time step size
        self.current_max_step_size = self.initial_time_step_size
        self.model_part.ProcessInfo[KratosMultiphysics.STEP] += 1

        if self.number_of_initial_time_steps < self.model_part.ProcessInfo[KratosMultiphysics.STEP]:
            self.current_max_step_size = self.max_time_step_size

            self.safety_factor *= 1.2
            if self.max_safety_factor < self.safety_factor:
                self.safety_factor = self.max_safety_factor

        # Advance
        self.current_step_size = self.__EstimateTimeStep(self.safety_factor, self.current_max_step_size)
        new_time = current_time + self.current_step_size
        self.model_part.CloneTimeStep(new_time)

        return new_time

    def SolveSolutionStep(self) -> bool:
        """Perform a local solution loop until the time step size shrinks to the target bounds."""
        # Note:
        # To stick to the original implementation, 3 initial steps are taken
        # without solving the system, though I'd assume we only need 2.
        while True and (3 < self.model_part.ProcessInfo[KratosMultiphysics.STEP]):
            self.__SolveLocalSolutionStep()

            # Check step size and revert if it's too large
            # Note: does hardcoding a 0.95 safety factor make sense here
            # or should it be variable?
            step_size_estimate = self.__EstimateTimeStep(0.95, self.current_max_step_size)

            if step_size_estimate < self.current_step_size:
                # Log failure
                self.__Log("Reducing time step ...")
                time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
                self.__Log("Failed time: {}".format(time))
                self.__Log("Failed step size: {}".format(self.current_step_size))

                # Clear step data
                self.fluid_solver.ReduceTimeStep(self.model_part, time)

                # Get reduced step
                time -= self.current_step_size
                self.safety_factor *= self.reduction_on_failure
                self.current_step_size = self.__EstimateTimeStep(
                    self.safety_factor,
                    self.current_max_step_size)
                time += self.current_step_size

                self.__Log("Reduced time: {}".format(time))
                self.__Log("Reduced step size: {}".format(self.current_step_size))

                self.fluid_solver.ReduceTimeStep(self.model_part, time)
            else:
                break

        return True

    def GetComputingModelPart(self) -> KratosMultiphysics.ModelPart:
        return self.model_part

    def ExportModelPart(self) -> None:
        file_name = pathlib.Path(self.settings["model_import_settings"]["input_file_name"].GetString() + ".out")
        KratosMultiphysics.ModelPartIO(
            file_name,
            KratosMultiphysics.IO.WRITE).WriteModelPart(self.model_part)

        self.__Log("Model part written to '{}'".format(file_name))

    def Finalize(self):
        super().Finalize()
        self.fluid_solver.Clear()

    ## EdgebasedLevelSetSolver specific methods.

    def _Redistance(self) -> None:
        if self.use_parallel_distance_calculation:
            self.distance_utils.Execute()
        else:
            self.distance_utils.CalculateDistances(
                self.model_part,
                KratosMultiphysics.DISTANCE,
                self.distance_size)

        # Make sure that at least one node has negative distance
        active_node_count = False
        for node in self.model_part.Nodes:
            if node.GetSolutionStepValue(KratosMultiphysics.DISTANCE) < 0.0:
                active_node_count = True
                break

        if not active_node_count:
            raise RuntimeError("At least 1 node must have negative DISTANCE")

    def _CreateLinearSolver(self):
        # Create the pressure linear solver
        linear_solver_configuration = self.settings["linear_solver_settings"]
        linear_solver = linear_solver_factory.ConstructSolver(linear_solver_configuration)
        return linear_solver

    def __MakeMatrixContainer(self) -> FreeSurface.MatrixContainer3D:
        if self.domain_size == 2:
            return FreeSurface.MatrixContainer2D()
        elif self.domain_size == 3:
            return FreeSurface.MatrixContainer3D()
        else:
            raise ValueError("Invalid domain size: {}".format(self.domain_size))

    def __MakeDistanceUtilities(self) -> KratosMultiphysics.ParallelDistanceCalculationProcess3D:
        if self.use_parallel_distance_calculation:
            parallel_distance_settings = KratosMultiphysics.Parameters("""{
                "max_levels" : 25,
                "max_distance" : 1.0
            }""")
            parallel_distance_settings["max_levels"].SetInt(self.extrapolation_layers)
            parallel_distance_settings["max_levels"].SetDouble(self.distance_size)

            if self.domain_size == 2:
                return KratosMultiphysics.ParallelDistanceCalculationProcess2D(
                    self.model_part,
                    parallel_distance_settings)

            elif self.domain_size == 3:
                return KratosMultiphysics.ParallelDistanceCalculationProcess3D(
                    self.model_part,
                    parallel_distance_settings)
            else:
                raise ValueError("Invalid domain size: {}".format(self.domain_size))
        else:
            if self.domain_size == 2:
                return FreeSurface.EdgebasedLevelsetAuxiliaryUtils2D()
            elif self.domain_size == 3:
                return FreeSurface.EdgebasedLevelsetAuxiliaryUtils3D()
            else:
                raise ValueError("Invalid domain size: {}".format(self.domain_size))

    def __MakeEdgeBasedLevelSet(self) -> FreeSurface.EdgeBasedLevelSet3D:
        if self.domain_size == 2:
            EdgeBasedLevelSet = FreeSurface.EdgeBasedLevelSet2D
        elif self.domain_size == 3:
            # Conditions' neighbours must be evaluated in 3D
            KratosMultiphysics.FindConditionsNeighboursProcess(
                self.model_part, self.domain_size, 10).Execute()
            EdgeBasedLevelSet = FreeSurface.EdgeBasedLevelSet3D
        else:
            raise ValueError("Invalid domain size: {}".format(self.domain_size))

        return EdgeBasedLevelSet(
            self.matrix_container,
            self.model_part,
            self.viscosity,
            self.density,
            self.use_mass_correction,
            self.stabdt_pressure_factor,
            self.stabdt_convection_factor,
            self.tau2_factor,
            self.assume_constant_pressure)

    def __EstimateTimeStep(self, safety_factor: float, max_time_step_size: float) -> float:
        step_size = self.fluid_solver.ComputeTimeStep(safety_factor, max_time_step_size)
        if max_time_step_size < step_size:
            step_size = max_time_step_size

        return step_size

    def __SolveLocalSolutionStep(self):
        """Solve the system in its current state without checking the time step."""
        if self.extrapolation_layers < 3:
            raise ValueError("Insufficient number of extrapolation layers: {} (minimum is 3)".format(self.extrapolation_layers))

        self.timer.Start("Calculate porous resistance law")
        self.fluid_solver.CalculatePorousResistanceLaw(int(self.compute_porous_resistance_law))
        self.timer.Stop("Calculate porous resistance law")

        self.timer.Start("Update fixed velocity values")
        self.fluid_solver.UpdateFixedVelocityValues()
        self.timer.Stop("Update fixed velocity values")

        self.timer.Start("Extrapolate values")
        self.fluid_solver.ExtrapolateValues(self.extrapolation_layers)
        self.timer.Stop("Extrapolate values")

        self.timer.Start("Convect distance")
        self.fluid_solver.ConvectDistance()
        self.timer.Stop("Convect distance")

        # Note: the original script kept track of two different step counts,
        # leading to some quirks with performing redistancing. The current
        # behaviour mimics the original one to satisfy the tests, but it
        # would probably make more sense to check (STEP % frequency == 0) instead.
        step = self.model_part.ProcessInfo[KratosMultiphysics.STEP]
        if 0 < (step - 4) and ((step - 4) % self.redistance_frequency) == 0:
            self.timer.Start("Redistance")
            self._Redistance()
            self.timer.Stop("Redistance")

        # Solve fluid
        self.timer.Start("Solve step 1")
        self.fluid_solver.SolveStep1()
        self.timer.Stop("Solve step 1")

        self.timer.Start("Solve step 2")
        self.fluid_solver.SolveStep2(self.pressure_linear_solver)
        self.timer.Stop("Solve step 2")

        self.timer.Start("Solve step 3")
        self.fluid_solver.SolveStep3()
        self.timer.Stop("Solve step 3")

    def __Log(self, message: str) -> None:
        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, message)


def CreateSolver(model: KratosMultiphysics.Model, parameters: KratosMultiphysics.Parameters) -> EdgeBasedLevelSetSolver:
    if not isinstance(model, KratosMultiphysics.Model):
        raise TypeError("Expecting a KratosMultiphysics.Model object, but got {}".format(type(model)))

    if not isinstance(parameters, KratosMultiphysics.Parameters):
        raise TypeError("Expecting a KratosMultiphysics.Parameters object, but got {}".format(type(parameters)))

    return EdgeBasedLevelSetSolver(model, parameters)
