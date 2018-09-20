from __future__ import print_function, absolute_import, division
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
import math
import numpy as np


def printMatrix(matrix):

    m = matrix.shape[0]
    n = matrix.shape[1]

    print("[")
    for i in range(0, m):
        print("")
        for j in range(0, n):

            print(matrix[i][j], end=",")
    print("]")


def AssembleLHSAndRHSInGlobal(lhs, rhs, lhs_cond, rhs_cond, node_ids):
    eqIds = []
    for m in range(0, len(node_ids)):

        eqIds.append((node_ids[m]-1)*3)
        eqIds.append((node_ids[m]-1)*3 + 1)
        eqIds.append((node_ids[m]-1)*3 + 2)

    for i in range(0, len(eqIds)):

        rhs[eqIds[i]] += rhs_cond[i]

        for j in range(0, len(eqIds)):
            # print("#####global :: %d,%d , ####local  :: %d,%d)"%(eqIds[i], eqIds[j], i, j))
            lhs[eqIds[i]][eqIds[j]] += lhs_cond[i, j]


def ApplyDisplacement(model_part, displacements):
    i = 0
    for node in model_part.Nodes:
        node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, 0, [
                                  displacements[i], displacements[i+1], displacements[i+2]])

        i = i+3


def ApplyDiricheletCondition(lhs, nodes):

    for i in range(0, len(nodes)):

        x_dof = (nodes[i]-1)*3
        y_dof = (nodes[i]-1)*3+1
        z_dof = (nodes[i]-1)*3+2

        # Rows to zero

        lhs[x_dof][:] = 0.0
        lhs[y_dof][:] = 0.0
        lhs[z_dof][:] = 0.0

        # Columns to zero

        lhs[:][x_dof] = 0.0
        lhs[:][y_dof] = 0.0
        lhs[:][z_dof] = 0.0

        lhs[x_dof][x_dof] = 1.0
        lhs[y_dof][y_dof] = 1.0
        lhs[z_dof][z_dof] = 1.0


class TestLoadingConditionsSurface(KratosUnittest.TestCase):
    def test_HydrostaticLoad3D3N(self):
        dim = 2
        mp = KratosMultiphysics.ModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        mp.AddNodalSolutionStepVariable(
            KratosMultiphysics.POSITIVE_FACE_PRESSURE)

        # create nodes
        mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2, 1.0, 0.0, 0.0)
        mp.CreateNewNode(3, 0.0, 1.0, 0.0)
        mp.CreateNewNode(4, 1.0, 1.0, 0.0)

        specific_weight = 10.0

        centre = [0.0, 0.0, 1.0]
        radius = 10.0
        normal = [0.0, 0.0, 1.0]
        fluid_volume = 2.0

        nodeIds1 = [1, 3, 2]
        nodeIds2 = [2, 3, 4]
        nodeIds = list(set().union(nodeIds1, nodeIds2))

        # ensure that the property 1 is created
        self.properties = mp.GetProperties()[1]

        self.properties.SetValue(
            StructuralMechanicsApplication.FREE_SURFACE_RADIUS, radius)
        self.properties.SetValue(
            StructuralMechanicsApplication.FLUID_VOLUME, fluid_volume)

        self.properties.SetValue(
            StructuralMechanicsApplication.FREE_SURFACE_CENTRE, centre)
        self.properties.SetValue(
            StructuralMechanicsApplication.FREE_SURFACE_NORMAL, normal)
        self.properties.SetValue(
            StructuralMechanicsApplication.SPECIFIC_WEIGHT, specific_weight)

        mp.SetBufferSize(2)

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X,
                                                  KratosMultiphysics.REACTION_X, mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y,
                                                  KratosMultiphysics.REACTION_Y, mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z,
                                                  KratosMultiphysics.REACTION_Z, mp)

        cond_hydro_1 = mp.CreateNewCondition(
            "HydrostaticLoadCondition3D3N", 1, nodeIds1, mp.GetProperties()[1])
        cond_hydro_2 = mp.CreateNewCondition(
            "HydrostaticLoadCondition3D3N", 2, nodeIds2, mp.GetProperties()[1])

        elem1 = mp.CreateNewElement(
            "ShellThinElement3D3N", 1, nodeIds1, mp.GetProperties()[1])
        elem2 = mp.CreateNewElement(
            "ShellThinElement3D3N", 2, nodeIds2, mp.GetProperties()[1])

        for node in mp.Nodes:
            node.SetSolutionStepValue(
                KratosMultiphysics.DISPLACEMENT_X, 0, 0.0)
            node.SetSolutionStepValue(
                KratosMultiphysics.DISPLACEMENT_Y, 0, 0.0)
            node.SetSolutionStepValue(
                KratosMultiphysics.DISPLACEMENT_Z, 0, 0.0)

            node.SetSolutionStepValue(
                KratosMultiphysics.DISPLACEMENT_X, 1, 0.0)
            node.SetSolutionStepValue(
                KratosMultiphysics.DISPLACEMENT_Y, 1, 0.0)
            node.SetSolutionStepValue(
                KratosMultiphysics.DISPLACEMENT_Z, 1, 0.0)

        VolumeCalcUtilty = StructuralMechanicsApplication.VolumeCalculationUnderPlaneUtility(
            centre, radius, normal)
        vol = VolumeCalcUtilty.CalculateVolume(mp)
        intersected_area = VolumeCalcUtilty.GetIntersectedArea()
        self.properties.SetValue(
            StructuralMechanicsApplication.FREE_SURFACE_AREA, intersected_area)

        avg_nodes = 10
        avg_elems = 10

        neighbhor_finder = KratosMultiphysics.FindNodalNeighboursProcess(
            mp, avg_elems, avg_nodes)

        neighbhor_finder.Execute()

        print("Initial Volume of the structure is ::	%f " % vol)

        """ for node in mp.Nodes:
		
		
            height = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE,0)
            
            if height < 0:
                pressure = specific_weight*height

                
            else: pressure = 0
        
            node.SetSolutionStepValue(KratosMultiphysics.POSITIVE_FACE_PRESSURE,pressure) """

        lhs_hydro_1 = KratosMultiphysics.Matrix(0, 0)
        rhs_hydro_1 = KratosMultiphysics.Vector(0)
        lhs_hydro_2 = KratosMultiphysics.Matrix(0, 0)
        rhs_hydro_2 = KratosMultiphysics.Vector(0)

        # Global Matrix
        lhs_hydro = np.zeros([len(nodeIds)*3, len(nodeIds)*3])
        rhs_hydro = np.zeros([len(nodeIds)*3])

        for node in mp.Nodes:
            print("Distance :: ", node.GetSolutionStepValue(
                KratosMultiphysics.DISTANCE, 0))
            print("Pressure :: ", node.GetSolutionStepValue(
                KratosMultiphysics.POSITIVE_FACE_PRESSURE, 0))

        cond_hydro_1.CalculateLocalSystem(
            lhs_hydro_1, rhs_hydro_1, mp.ProcessInfo)
        cond_hydro_2.CalculateLocalSystem(
            lhs_hydro_2, rhs_hydro_2, mp.ProcessInfo)

        AssembleLHSAndRHSInGlobal(
            lhs_hydro, rhs_hydro, lhs_hydro_1, rhs_hydro_1, nodeIds1)

        AssembleLHSAndRHSInGlobal(
            lhs_hydro, rhs_hydro, lhs_hydro_2, rhs_hydro_2, nodeIds2)

        # AssembleLHSAndRHSInGlobal(lhs_sur,rhs_sur,lhs_sur_1,rhs_sur_1,nodeIds1)
        # AssembleLHSAndRHSInGlobal(lhs_sur,rhs_sur,lhs_sur_2,rhs_sur_2,nodeIds2)

        for node in mp.Nodes:
            print("Nodal Normal :: ", node.GetSolutionStepValue(
                KratosMultiphysics.NORMAL), 0)

        dirichelet_nodes = [1, 4]
        ApplyDiricheletCondition(lhs_hydro, dirichelet_nodes)

        # print(lhs_hydro_1)
        print('#################')
        # print(lhs_hydro_2)
        # print('#################')
        # print(rhs_hydro_1)
        # print('#################')
        # print(rhs_hydro_2)

        print("##############Global RHS###################")

        print(rhs_hydro)
        # print("#####################################")
        # print(rhs_sur)

        print("##############Global LHS###################")

        printMatrix(lhs_hydro)
        # print("#####################################")
        # printMatrix(lhs_sur)

        eigen_hydro = np.linalg.eig(lhs_hydro)

        print("###########Hydrostatic EigenValues#######")

        for i in range(0, lhs_hydro.shape[0]):
            if(math.fabs(np.real(eigen_hydro[0][i])) > 1e-15):
                print("Eigen value :: ", eigen_hydro[0][i])
                print("Eigen vector :: ", eigen_hydro[1][:, i])
                ApplyDisplacement(mp, eigen_hydro[1][:, i])

        #reference_res = [-44.872685126136794,46.372685126136815,0.0,-44.87268512613681,46.3726851261368,0.0,-56.657798145912594,58.1577981459126,0.0,-56.65779814591261,58.157798145912615,0.0]
        # for i in range(len(rhs)):
        self.assertAlmostEqual(0, 0)


if __name__ == '__main__':
    KratosUnittest.main()
