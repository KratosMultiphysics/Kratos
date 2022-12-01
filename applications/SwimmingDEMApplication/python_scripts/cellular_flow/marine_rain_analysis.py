import KratosMultiphysics as Kratos
from KratosMultiphysics import Vector, Parameters
import KratosMultiphysics.SwimmingDEMApplication as SDEM
import math
import KratosMultiphysics.SwimmingDEMApplication.swimming_DEM_analysis as swimming_DEM_analysis
BaseAnalysis = swimming_DEM_analysis.SwimmingDEMAnalysis

class MarineRainAnalysis(BaseAnalysis):
    def __init__(self, varying_parameters = Parameters("{}")):
        BaseAnalysis.__init__(self, varying_parameters)
        self.project_parameters.do_search_dem_neighbours = False
        self.project_parameters.vorticity_calculation_type = 0
        self.project_parameters.do_impose_flow_from_field = True

        self.project_parameters.L = 0.1
        self.project_parameters.U = 0.3
        self.project_parameters.k = 2.72
        self.project_parameters.omega = math.pi

    def GetFieldUtility(self):
        self.flow_field = SDEM.CellularFlowField(self.project_parameters.L, self.project_parameters.U, self.project_parameters.k, self.project_parameters.omega)
        space_time_set = SDEM.SpaceTimeSet()
        self.field_utility = SDEM.FluidFieldUtility(space_time_set, self.flow_field, 1000.0, 1e-6)
        return self.field_utility

    def PerformZeroStepInitializations(self):
        import random
        #from random import randint
        L = self.project_parameters.L
        #U = self.project_parameters.U
        #k = self.project_parameters.k
        #omega = self.project_parameters.omega
        N_positions = 100
        L *= math.pi
        possible_xs = [2 * L * (i + 1) / N_positions for i in range(N_positions)]
        i_position = 0

        for node in self.spheres_model_part.Nodes:
            rand_x = 2 * L * random.random()
            rand_y = self.project_parameters.BoundingBoxMinY + (
                     self.project_parameters.BoundingBoxMaxY -
                     self.project_parameters.BoundingBoxMinY
                     ) * random.random()
            init_x = node.X
            init_y = node.Y
            #rand_x = possible_xs[randint(0, N_positions - 1)]
            #rand_y = possible_xs[randint(0, N_positions - 1)]
            #rand_x = possible_xs[min(i_position // N_positions, N_positions - 1)]
            #rand_y = possible_xs[min(i_position // N_positions, N_positions - 1)]
            rand_x = possible_xs[i_position // N_positions]
            rand_y = possible_xs[i_position % N_positions]
            node.X = rand_x
            node.Y = rand_y
            node.SetSolutionStepValue(Kratos.DISPLACEMENT_X, rand_x - init_x)
            node.SetSolutionStepValue(Kratos.DISPLACEMENT_Y, rand_y - init_y)
            i_position += 1

        for node in self.spheres_model_part.Nodes:
            vel= Vector(3)
            coor= Vector(3)
            coor[0]=node.X
            coor[1]=node.Y
            coor[2]=node.Z
            self.flow_field.Evaluate(0.0,coor,vel,0)
            node.SetSolutionStepValue(Kratos.VELOCITY_X, vel[0])
            node.SetSolutionStepValue(Kratos.VELOCITY_Y, vel[1])
            node.SetSolutionStepValue(Kratos.VELOCITY_Z, vel[2])
            node.SetSolutionStepValue(Kratos.VELOCITY_OLD_X, vel[0])
            node.SetSolutionStepValue(Kratos.VELOCITY_OLD_Y, vel[1])
            node.SetSolutionStepValue(Kratos.VELOCITY_OLD_Z, vel[2])
            node.SetSolutionStepValue(Kratos.SLIP_VELOCITY_X, vel[0])
            node.SetSolutionStepValue(Kratos.SLIP_VELOCITY_Y, vel[1])
            node.SetSolutionStepValue(Kratos.SLIP_VELOCITY_Z, vel[2])

    def FluidSolve(self, time = 'None', solve_system=True):
        pass
