import KratosMultiphysics

import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory

class DistanceReinitialization:
    def __init__(self,model_part,params):
        self.params = params
        self.model_part = model_part
        self._CreateDistanceInitializationProcess()
         
    def _CreateDistanceInitializationProcess(self):
        if not self.params.Has("distance_reinitialization_type"):
            self.params.AddString("distance_reinitialization_type", "variational")

        self._reinitialization_type = self.params["distance_reinitialization_type"].GetString()
        if (self._reinitialization_type == "variational"):
            self._ValidateVariationalDefaultSettings()
            maximum_iterations = self.params["max_iterations"].GetInt()
            linear_solver = self._GetLinearSolver()
            process_flag = self._GetVariationalCalculateExactDistancesFlag()
            self.distance_reinitialization_process = self._ConstructVariationalProcess(maximum_iterations, linear_solver, process_flag)

        elif (self._reinitialization_type == "parallel"):
            self._ValidateParallelDefaultSettings()
            self.distance_reinitialization_process = self._ConstructParallelProcess()
            
        elif (self._reinitialization_type == "none"):
            KratosMultiphysics.Logger.PrintInfo("DistanceReinitialization", "Redistancing is turned off.")
        else:
            raise Exception("Please use a valid distance reinitialization type or set it as \'none\'. Valid types are: \'variational\' and \'parallel\'.")

    def _ValidateVariationalDefaultSettings(self):
        default_settings = KratosMultiphysics.Parameters("""{
            "distance_reinitialization_type" : "variational",
            "max_iterations" : 2,
            "linear_solver_settings" : {
                "solver_type" : "amgcl"
            },
            "calculate_exact_distances_to_plane" : true
        }""")
        self.params.ValidateAndAssignDefaults(default_settings)

    def _ValidateParallelDefaultSettings(self):
        default_settings = KratosMultiphysics.Parameters("""{
            "distance_reinitialization_type" : "parallel",
            "max_layers" : 25,
            "max_distance": 1.0,
            "calculate_exact_distances_to_plane" : true
        }""")
        self.params.ValidateAndAssignDefaults(default_settings)

    def _GetVariationalCalculateExactDistancesFlag(self):
        if self.params["calculate_exact_distances_to_plane"].GetBool():
            return KratosMultiphysics.VariationalDistanceCalculationProcess2D.CALCULATE_EXACT_DISTANCES_TO_PLANE
        else:
            return KratosMultiphysics.VariationalDistanceCalculationProcess2D.CALCULATE_EXACT_DISTANCES_TO_PLANE.AsFalse()

    def _GetParallelCalculateExactDistancesFlag(self):
        if self.params["calculate_exact_distances_to_plane"].GetBool():
            return KratosMultiphysics.ParallelDistanceCalculator2D.CALCULATE_EXACT_DISTANCES_TO_PLANE
        else:
            return KratosMultiphysics.ParallelDistanceCalculator3D.CALCULATE_EXACT_DISTANCES_TO_PLANE.AsFalse()

    def _GetLinearSolver(self):
        linear_solver_configuration = self.params["linear_solver_settings"]
        return linear_solver_factory.ConstructSolver(linear_solver_configuration)

    def _ConstructVariationalProcess(self, maximum_iterations, linear_solver, process_flag):
        if self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
            return KratosMultiphysics.VariationalDistanceCalculationProcess2D(
                self.model_part,
                linear_solver,
                maximum_iterations,
                process_flag)
        else:
            return KratosMultiphysics.VariationalDistanceCalculationProcess3D(
                self.model_part,
                linear_solver,
                maximum_iterations,
                process_flag)
            
    def _ConstructParallelProcess(self):
        if self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
            return KratosMultiphysics.ParallelDistanceCalculator2D()
        else:
            return KratosMultiphysics.ParallelDistanceCalculator3D()

    def ReinitializeDistances(self):
        if (self._reinitialization_type == "variational"):
            self.distance_reinitialization_process.Execute()

        elif (self._reinitialization_type == "parallel"):
            layers = self.params["max_layers"].GetInt()
            max_distance = self.params["max_distance"].GetDouble()
            process_flag = self._GetParallelCalculateExactDistancesFlag()
            self.distance_reinitialization_process.CalculateDistances(
                self.model_part,
                KratosMultiphysics.DISTANCE,
                KratosMultiphysics.NODAL_AREA,
                layers,
                max_distance,
                process_flag)

        if (self._reinitialization_type != "none"):
            KratosMultiphysics.Logger.PrintInfo("DistanceReinitialization", "Redistancing process is finished.")
