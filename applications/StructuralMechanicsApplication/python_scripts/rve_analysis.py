import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

# Importing other libraries
import math

class RVEAnalysis(StructuralMechanicsAnalysis):
    def __init__(self, model, project_parameters):

        # input parameters of the analysis
        self.boundary_mp_name = project_parameters["rve_settings"]["boundary_mp_name"].GetString()
        self.averaging_mp_name = project_parameters["rve_settings"]["averaging_mp_name"].GetString()
        self.print_rve_post = project_parameters["rve_settings"]["print_rve_post"].GetBool()
        self.perturbation = project_parameters["rve_settings"]["perturbation"].GetDouble()

        self.averaging_volume = -1.0  # it will be computed in initialize
        domain_size = project_parameters["solver_settings"]["domain_size"].GetInt()
        if(domain_size == 2):
            self.strain_size = 3
        else:
            self.strain_size = 6

        # Pseudo time to be used for output
        self.time = 0.0

        self.populate_search_eps = 1e-4 ##tolerance in finding which conditions belong to the surface (will be multiplied by the lenght of the diagonal)
        self.geometrical_search_tolerance = 1e-4 #tolerance to be used in the search of the condition it falls into

        super().__init__(model, project_parameters)

    # Here populate the submodelparts to be used for periodicity
    def ModifyInitialGeometry(self):
        super().ModifyInitialGeometry()

        boundary_mp = self.model[self.boundary_mp_name]
        averaging_mp = self.model[self.averaging_mp_name]

        # Construct auxiliary modelparts
        self.min_corner, self.max_corner = self._DetectBoundingBox(averaging_mp)
        self._ConstructFaceModelParts(self.min_corner, self.max_corner, boundary_mp)

        self.averaging_volume = (self.max_corner[0]-self.min_corner[0]) * (self.max_corner[1]-self.min_corner[1]) * (self.max_corner[2]-self.min_corner[2])
        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "RVE undeformed averaging volume = ", self.averaging_volume)

    def InitializeSolutionStep(self):
        raise Exception("Should use the _CustomInitializeSolutionStep instead of this")

    def __CustomInitializeSolutionStep(self, strain, boundary_mp, averaging_mp):
        #reset position
        KratosMultiphysics.VariableUtils().UpdateCurrentToInitialConfiguration(averaging_mp.Nodes)

        self.ApplyBoundaryConditions()  # here the processes are called

        # construct MPCs according to the provided strain
        self._ApplyPeriodicity(strain, averaging_mp, boundary_mp)

        # apply BCs for RVE according to the provided strain
        self._ApplyMinimalConstraints(
            averaging_mp, strain, self.min_corner, self.max_corner)

        self.ChangeMaterialProperties()  # this is normally empty

        self._GetSolver().InitializeSolutionStep()

        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "STEP: ", self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.STEP])
        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "TIME: ", self.time)

    def RunSolutionLoop(self):

        perturbation = self.perturbation

        boundary_mp = self.model[self.boundary_mp_name]
        averaging_mp = self.model[self.averaging_mp_name]

        stress_and_strain = []
        if(self.strain_size == 3):  # 2D case - ordering s00 s11 s01
            stress_and_strain.append(self._ComputeEquivalentStress(0, 0, perturbation, boundary_mp, averaging_mp))
            stress_and_strain.append(self._ComputeEquivalentStress(1, 1, perturbation, boundary_mp, averaging_mp))
            stress_and_strain.append(self._ComputeEquivalentStress(0, 1, perturbation, boundary_mp, averaging_mp))
        elif(self.strain_size == 6):  # 3D case - ordering:  s00 s11 s22 s01 s12 s02
            stress_and_strain.append(self._ComputeEquivalentStress(0, 0, perturbation, boundary_mp, averaging_mp))
            stress_and_strain.append(self._ComputeEquivalentStress(1, 1, perturbation, boundary_mp, averaging_mp))
            stress_and_strain.append(self._ComputeEquivalentStress(2, 2, perturbation, boundary_mp, averaging_mp))
            stress_and_strain.append(self._ComputeEquivalentStress(0, 1, perturbation, boundary_mp, averaging_mp))
            stress_and_strain.append(self._ComputeEquivalentStress(1, 2, perturbation, boundary_mp, averaging_mp))
            stress_and_strain.append(self._ComputeEquivalentStress(0, 2, perturbation, boundary_mp, averaging_mp))

        C = self._ComputeEquivalentElasticTensor(stress_and_strain, perturbation)
        averaging_mp.SetValue(KratosMultiphysics.StructuralMechanicsApplication.ELASTICITY_TENSOR, C)
        self._MatrixOutput(C)

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

        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Boundng box detected")
        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Min. corner = ", min_corner)
        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Max. corner = ", max_corner)

        return min_corner, max_corner

    def __PopulateMp(self, face_name, coordinate, component, eps, mp):
        if mp.NumberOfConditions() == 0:
            raise Exception("Boundary_mp is expected to have conditions and has none")

        mp = mp.GetRootModelPart()

        if not mp.HasSubModelPart(face_name):
            mp.CreateSubModelPart(face_name)
        face_mp = mp.GetSubModelPart(face_name)

        for cond in mp.Conditions:
            xc = cond.GetGeometry().Center()
            if abs(xc[component]-coordinate) < eps:
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
        diag_lenght = math.sqrt(diag_vect[0]**2+diag_vect[1]**2+diag_vect[2]**2 )
        eps = self.populate_search_eps*diag_lenght

        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.SLAVE, False, mp.Nodes)
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.MASTER, False, mp.Nodes)

        # Populate the slave faces
        self.max_x_face = self.__PopulateMp("max_x_face", max_corner[0], 0, eps, mp)
        self.max_y_face = self.__PopulateMp("max_y_face", max_corner[1], 1, eps, mp)
        self.max_z_face = self.__PopulateMp("max_z_face", max_corner[2], 2, eps, mp)

        # First populate the master faces (min)
        self.min_x_face = self.__PopulateMp("min_x_face", min_corner[0], 0, eps, mp)
        self.min_y_face = self.__PopulateMp("min_y_face", min_corner[1], 1, eps, mp)
        self.min_z_face = self.__PopulateMp("min_z_face", min_corner[2], 2, eps, mp)

        if self.min_x_face.NumberOfConditions() == 0:
            raise Exception("min_x_face has 0 conditions")
        if self.min_y_face.NumberOfConditions() == 0:
            raise Exception("min_y_face has 0 conditions")
        if self.min_z_face.NumberOfConditions() == 0:
            raise Exception("min_z_face has 0 conditions")

    def _SelectClosestNode(self, mp, coords):
        min_distance = 1e30
        selected_node = 0
        for node in mp.Nodes:
            dx = node.X0 - coords[0]
            dy = node.Y0 - coords[1]
            dz = node.Z0 - coords[2]
            d = dx**2 + dy**2 + dz**2

            if(d < min_distance):
                selected_node = node
                min_distance = d

        return selected_node

    # prescribed conditions to avoid rigid body motions
    def _ApplyMinimalConstraints(self, mp, strain, min_corner, max_corner):
        aux = KratosMultiphysics.Array3()

        # point coinciding with the min_corner
        node = self._SelectClosestNode(mp, min_corner)
        node.Fix(KratosMultiphysics.DISPLACEMENT_X)
        node.Fix(KratosMultiphysics.DISPLACEMENT_Y)
        node.Fix(KratosMultiphysics.DISPLACEMENT_Z)

        coords_min_corner = KratosMultiphysics.Array3(node)
        coords_min_corner[0] = node.X0
        coords_min_corner[1] = node.Y0
        coords_min_corner[2] = node.Z0

        disp_min_corner = strain*coords_min_corner
        node.SetSolutionStepValue(
            KratosMultiphysics.DISPLACEMENT, 0, disp_min_corner)

    def _ComputeEquivalentElasticTensor(self, stress_and_strain, perturbation):
        C = KratosMultiphysics.Matrix(self.strain_size, self.strain_size)

        for j in range(len(stress_and_strain)):
            stress = stress_and_strain[j][0]
            for i in range(self.strain_size):
                C[i, j] = stress[i]

        inverse_perturbation = KratosMultiphysics.Matrix(self.strain_size, self.strain_size)
        inverse_perturbation.fill(0.0)
        if(self.strain_size == 3):
            inverse_perturbation[0, 0] = 1.0/perturbation
            inverse_perturbation[1, 1] = 1.0/perturbation
            inverse_perturbation[2, 2] = 0.5/perturbation
        else:
            inverse_perturbation[0, 0] = 1.0/perturbation
            inverse_perturbation[1, 1] = 1.0/perturbation
            inverse_perturbation[2, 2] = 1.0/perturbation
            inverse_perturbation[3, 3] = 0.5/perturbation
            inverse_perturbation[4, 4] = 0.5/perturbation
            inverse_perturbation[5, 5] = 0.5/perturbation

        C = C*inverse_perturbation

        return C

    def _MatrixOutput(self, C, filename="rve_elasticity_tensor.txt"):
        f = open(filename, 'w')

        if(self.strain_size == 3):  # 2D
            f.write(str(C[0, 0]) + " " + str(C[0, 1]) +
                    " " + str(C[0, 2]) + "\n")
            f.write(str(C[1, 0]) + " " + str(C[1, 1]) +
                    " " + str(C[1, 2]) + "\n")
            f.write(str(C[2, 0]) + " " + str(C[2, 1]) +
                    " " + str(C[2, 2]) + "\n")

        elif(self.strain_size == 6):
            for i in range(6):
                f.write(str(C[i, 0]) + " " + str(C[i, 1]) + " " + str(C[i, 2]) +
                        " " + str(C[i, 3]) + " " + str(C[i, 4]) + " " + str(C[i, 5]) + "\n")

        f.close()

    def _ComputeEquivalentStress(self, i, j, perturbation, boundary_mp, averaging_mp):
        # Here use a pseudotime for output
        self.time = self.time + 1.0
        averaging_mp.GetRootModelPart().CloneTimeStep(self.time)

        strain = KratosMultiphysics.Matrix(3, 3)
        strain.fill(0.0)

        strain[i, j] = perturbation
        strain[j, i] = perturbation

        strain_vector = KratosMultiphysics.Vector(self.strain_size)
        if(self.strain_size == 2):
            strain_vector[0] = strain[0, 0]
            strain_vector[1] = strain[1, 1]
            strain_vector[2] = 2.0*strain[1, 2]
        elif(self.strain_size == 6):
            strain_vector[0] = strain[0, 0]
            strain_vector[1] = strain[1, 1]
            strain_vector[2] = strain[2, 2]
            strain_vector[3] = 2.0*strain[0, 1]
            strain_vector[4] = 2.0*strain[1, 2]
            strain_vector[5] = 2.0*strain[0, 2]

        self.__CustomInitializeSolutionStep(strain, boundary_mp, averaging_mp)

        self._GetSolver().Predict()

        self._GetSolver().SolveSolutionStep()
        process_info = averaging_mp.ProcessInfo
        avg_stress = KratosMultiphysics.Vector(self.strain_size)
        avg_stress.fill(0.0)
        measured_volume = 0.0

        for elem in averaging_mp.Elements:
            tmp = elem.CalculateOnIntegrationPoints(KratosMultiphysics.PK2_STRESS_VECTOR, process_info)
            ngauss = len(tmp)
            A = elem.GetGeometry().Area()
            measured_volume += A
            # TODO: this is only valid for gauss points with the same weight. should be generalized
            Agauss = A/ngauss
            for item in tmp:
                avg_stress = avg_stress + item*Agauss

        self._GetSolver().Clear()

        KratosMultiphysics.Logger.PrintInfo(
            self._GetSimulationName(), "Measured volume = ", measured_volume)

        avg_stress /= self.averaging_volume

        KratosMultiphysics.Logger.PrintInfo(
            self._GetSimulationName(), "Applied strain = ", strain_vector)
        KratosMultiphysics.Logger.PrintInfo(
            self._GetSimulationName(), "Average stress = ", avg_stress)

        if self.print_rve_post:
            self.OutputSolutionStep()

        # Reset position of nodes
        KratosMultiphysics.VariableUtils().UpdateCurrentToInitialConfiguration(averaging_mp.Nodes)
        zero = KratosMultiphysics.Vector(3)
        zero[0] = 0.0
        zero[1] = 0.0
        zero[2] = 0.0
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.DISPLACEMENT, zero, averaging_mp.Nodes)
        return avg_stress, strain_vector

    def _ApplyPeriodicity(self, strain, volume_mp, boundary_mp):
        # clear
        for constraint in volume_mp.GetRootModelPart().MasterSlaveConstraints:
            constraint.Set(KratosMultiphysics.TO_ERASE)
        volume_mp.GetRootModelPart().RemoveMasterSlaveConstraintsFromAllLevels(
            KratosMultiphysics.TO_ERASE)

        dx = self.max_corner[0] - self.min_corner[0]
        dy = self.max_corner[1] - self.min_corner[1]
        dz = self.max_corner[2] - self.min_corner[2]

        periodicity_utility = KratosMultiphysics.RVEPeriodicityUtility(self._GetSolver().GetComputingModelPart())

        # assign periodicity to faces
        search_tolerance = self.geometrical_search_tolerance
        periodicity_utility.AssignPeriodicity(self.min_x_face, self.max_x_face, strain, KratosMultiphysics.Vector([dx, 0.0, 0.0]),search_tolerance)
        periodicity_utility.AssignPeriodicity(self.min_y_face, self.max_y_face, strain, KratosMultiphysics.Vector([0.0, dy, 0.0]),search_tolerance)
        periodicity_utility.AssignPeriodicity(self.min_z_face, self.max_z_face, strain, KratosMultiphysics.Vector([0.0, 0.0, dz]),search_tolerance)

        periodicity_utility.Finalize(KratosMultiphysics.DISPLACEMENT)

        # start from the exact solution in the case of a constant strain
        x = KratosMultiphysics.Array3()
        for node in volume_mp.Nodes:
            x[0] = node.X0
            x[1] = node.Y0
            x[2] = node.Z0
            d = strain*x
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, 0, d)
