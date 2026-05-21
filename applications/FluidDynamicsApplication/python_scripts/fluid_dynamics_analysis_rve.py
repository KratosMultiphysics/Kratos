import KratosMultiphysics

from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis

# Importing other libraries
import math

class FluidDynamicsAnalysisRve(FluidDynamicsAnalysis):
    def __init__(self, model, project_parameters):
        # Validate RVE settings
        input_rve_settings = project_parameters["rve_settings"]
        input_rve_settings.ValidateAndAssignDefaults(self.GetFluidRVEDefaultParameters())

        # Check input domain size
        self.domain_size = project_parameters["solver_settings"]["domain_size"].GetInt()
        if self.domain_size != 2 and self.domain_size != 3:
            err_msg = "Wrong 'domain_size' value {}. Expected 2 or 3.".format(self.domain_size)
            raise ValueError(err_msg)

        # Input parameters of the analysis
        self.fix_pressure = project_parameters["rve_settings"]["fix_pressure"].GetBool()
        self.boundary_mp_name = project_parameters["rve_settings"]["boundary_mp_name"].GetString()
        if self.boundary_mp_name == "":
            raise Exception("'boundary_mp_name' in 'rve_settings' is expected to be provided by the user.")
        self.averaging_mp_name = project_parameters["rve_settings"]["averaging_mp_name"].GetString()
        if self.averaging_mp_name == "":
            KratosMultiphysics.Logger.PrintWarning("'averaging_mp_name' is not provided. Taking 'model_part_name' in 'solver_settings'.")
            self.averaging_mp_name = project_parameters["solver_settings"]["model_part_name"].GetString()

        # Dictionary containing the jump data for periodicity
        self.strains = {}

        self.strains["00"] = project_parameters["rve_settings"]["jump_XX"].GetDouble()
        self.strains["11"] = project_parameters["rve_settings"]["jump_YY"].GetDouble()
        self.strains["01"] = project_parameters["rve_settings"]["jump_XY"].GetDouble()
        if(self.domain_size == 3):
            self.strains["22"] = project_parameters["rve_settings"]["jump_ZZ"].GetDouble()
            self.strains["02"] = project_parameters["rve_settings"]["jump_XZ"].GetDouble()
            self.strains["12"] = project_parameters["rve_settings"]["jump_YZ"].GetDouble()

        self.strain = KratosMultiphysics.Matrix(self.domain_size, self.domain_size)
        self.strain.fill(0.0)

        # Create and fill "strain" matrix
        for i in range(self.domain_size):
            for j in range(i ,self.domain_size):
                strainIndex = str(i) + str(j)

                self.strain[i, j] = self.strains[strainIndex]
                self.strain[j, i] = self.strains[strainIndex]

        self.populate_search_eps = 1e-4 ##tolerance in finding which conditions belong to the surface (will be multiplied by the length of the diagonal)
        self.geometrical_search_tolerance = 1e-4 #tolerance to be used in the search of the condition it falls into

        super().__init__(model, project_parameters)

    @classmethod
    def GetFluidRVEDefaultParameters(cls):
        # Settings string in json format
        default_settings = KratosMultiphysics.Parameters("""{
            "boundary_mp_name": "",
            "averaging_mp_name": "",
            "fix_pressure": true,
            "jump_XX": 0,
            "jump_XY": 1e-6,
            "jump_YY": 0,
            "jump_XZ": 1e-6,
            "jump_ZZ": 0,
            "jump_YZ": 1e-6
        }""")

        return default_settings

    def ModifyInitialGeometry(self):
        ''' Here populate the submodelparts to be used for periodicity'''
        super().ModifyInitialGeometry()
        boundary_mp = self.model[self.boundary_mp_name]
        self.averaging_mp = self.model[self.averaging_mp_name]

        # Construct auxiliary modelparts
        self.min_corner, self.max_corner = self._DetectBoundingBox(self.averaging_mp)
        self._ConstructFaceModelParts(self.min_corner, self.max_corner, boundary_mp)

    def ModifyAfterSolverInitialize(self):
        super().ModifyAfterSolverInitialize()

        self._ApplyPeriodicity(self.strain, self.averaging_mp)

        self._ComputeNodalArea() # Needed for average velocity correction

    def ApplyBoundaryConditions(self):
        super().ApplyBoundaryConditions()

        # Here the pressure is fixed in one node to ease the solver. Can be turned off in parameters.json
        if self.fix_pressure:
            model_part = self._GetSolver().GetComputingModelPart()

            for node in model_part.Nodes:
                x = node.X
                y = node.Y
                z = node.Z

                if x < self.min_corner[0] + self.geometrical_search_tolerance and y < self.min_corner[1] + self.geometrical_search_tolerance:
                    if self.domain_size == 3 and z < self.min_corner[2] + self.geometrical_search_tolerance:
                        node.Fix(KratosMultiphysics.PRESSURE)
                        node.SetSolutionStepValue(KratosMultiphysics.PRESSURE, 0.0)
                        break
                    else: #have to think about a better way of doing this...
                        node.Fix(KratosMultiphysics.PRESSURE)
                        node.SetSolutionStepValue(KratosMultiphysics.PRESSURE, 0.0)
                        break

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()

        # Ensures that the avg velocity in the RVE domain remains 0
        self._AvgVelocityCorrection()

    def _DetectBoundingBox(self, mp):
        min_corner = KratosMultiphysics.Array3()
        min_corner[0] = 1e20
        min_corner[1] = 1e20
        min_corner[2] = 1e20

        max_corner = KratosMultiphysics.Array3()
        max_corner[0] = -1e20
        max_corner[1] = -1e20
        max_corner[2] = -1e20

        for node in mp.Nodes:
            x = node.X
            min_corner[0] = min(min_corner[0], x)
            max_corner[0] = max(max_corner[0], x)

            y = node.Y
            min_corner[1] = min(min_corner[1], y)
            max_corner[1] = max(max_corner[1], y)

            z = node.Z
            min_corner[2] = min(min_corner[2], z)
            max_corner[2] = max(max_corner[2], z)

        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Bounding box detected")
        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Min. corner = ", min_corner)
        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Max. corner = ", max_corner)

        return min_corner, max_corner

    @classmethod
    def __PopulateMp(self, face_name, coordinate, component, eps, mp):
        if mp.NumberOfConditions() == 0:
            raise Exception("Boundary_mp is expected to have conditions and has none")

        mp = mp.GetRootModelPart()

        if not mp.HasSubModelPart(face_name):
            mp.CreateSubModelPart(face_name)
        face_mp = mp.GetSubModelPart(face_name)

        for cond in mp.Conditions:
            xc = cond.GetGeometry().Center()
            if abs(xc[component] - coordinate) < eps:
                face_mp.AddCondition(cond)

        node_ids = set()
        for cond in face_mp.Conditions:
            for node in cond.GetNodes():
                if(not node.Is(KratosMultiphysics.SLAVE)):
                    node_ids.add(node.Id)
                    node.Set(KratosMultiphysics.SLAVE)

        face_mp.AddNodes(list(node_ids))
        return face_mp

    def _ConstructFaceModelParts(self, min_corner, max_corner, mp):
        diag_vect = max_corner - min_corner
        diag_lenght = math.sqrt(diag_vect[0]**2 + diag_vect[1]**2 + diag_vect[2]**2 )
        eps = self.populate_search_eps * diag_lenght

        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.SLAVE, False, mp.Nodes)
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.MASTER, False, mp.Nodes)

        # Populate the slave faces
        self.max_x_face = self.__PopulateMp("max_x_face", max_corner[0], 0, eps, mp)
        self.max_y_face = self.__PopulateMp("max_y_face", max_corner[1], 1, eps, mp)
        if self.domain_size == 3:
            self.max_z_face = self.__PopulateMp("max_z_face", max_corner[2], 2, eps, mp)

        # First populate the master faces (min)
        self.min_x_face = self.__PopulateMp("min_x_face", min_corner[0], 0, eps, mp)
        self.min_y_face = self.__PopulateMp("min_y_face", min_corner[1], 1, eps, mp)
        if self.domain_size == 3:
            self.min_z_face = self.__PopulateMp("min_z_face", min_corner[2], 2, eps, mp)

        if self.min_x_face.NumberOfConditions() == 0:
            raise Exception("min_x_face has 0 conditions")
        if self.min_y_face.NumberOfConditions() == 0:
            raise Exception("min_y_face has 0 conditions")
        if self.domain_size == 3 and self.min_z_face.NumberOfConditions() == 0:
            raise Exception("min_z_face has 0 conditions")

    def _ApplyPeriodicity(self, strain, volume_mp):
        # clear
        for constraint in volume_mp.GetRootModelPart().MasterSlaveConstraints:
            constraint.Set(KratosMultiphysics.TO_ERASE)
        volume_mp.GetRootModelPart().RemoveMasterSlaveConstraintsFromAllLevels(
            KratosMultiphysics.TO_ERASE)

        dx = self.max_corner[0] - self.min_corner[0]
        dy = self.max_corner[1] - self.min_corner[1]
        if self.domain_size == 3:
            dz = self.max_corner[2] - self.min_corner[2]

        periodicity_utility = KratosMultiphysics.RVEPeriodicityUtility(self._GetSolver().GetComputingModelPart())

        # assign periodicity to faces
        direction_x = KratosMultiphysics.Vector([dx, 0.0, 0.0]) if self.domain_size == 3 else KratosMultiphysics.Vector([dx, 0.0])
        periodicity_utility.AssignPeriodicity(self.min_x_face, self.max_x_face, strain, direction_x, self.geometrical_search_tolerance)
        direction_y = KratosMultiphysics.Vector([0.0, dy, 0.0]) if self.domain_size == 3 else KratosMultiphysics.Vector([0.0, dy])
        periodicity_utility.AssignPeriodicity(self.min_y_face, self.max_y_face, strain, direction_y, self.geometrical_search_tolerance)
        if self.domain_size == 3:
            direction_z = KratosMultiphysics.Vector([0.0, 0.0, dz])
            periodicity_utility.AssignPeriodicity(self.min_z_face, self.max_z_face, strain, direction_z, self.geometrical_search_tolerance)

        periodicity_utility.Finalize(KratosMultiphysics.VELOCITY)

    def _ComputeNodalArea(self):
        nodal_area_process = KratosMultiphysics.CalculateNodalAreaProcess(
            self.averaging_mp,
            self.domain_size)
        nodal_area_process.Execute()

    def _AvgVelocityCorrection(self):
        v_avg = KratosMultiphysics.Array3()
        v_avg[0] = 0.0
        v_avg[1] = 0.0
        v_avg[2] = 0.0
        w_tot = 0.0

        for node in self.averaging_mp.Nodes:
            w = node.GetSolutionStepValue(KratosMultiphysics.NODAL_AREA)
            v_avg += w * node.GetSolutionStepValue(KratosMultiphysics.VELOCITY)
            w_tot += w
        v_avg /= w_tot

        for node in self.averaging_mp.Nodes:
            v_new = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY) - v_avg
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, v_new)

class FluidDynamicsAnalysisRVE(FluidDynamicsAnalysisRve):
    def __init__(self, model, project_parameters):
        KratosMultiphysics.Logger.PrintWarning("'FluidDynamicsAnalysisRVE' is deprecated. Use 'FluidDynamicsAnalysisRve' instead.")
        super().__init__(model, project_parameters)