from KratosMultiphysics import *
from KratosMultiphysics.SwimmingDEMApplication import *
from KratosMultiphysics.DEMApplication import *
import math
import swimming_DEM_analysis
BaseAnalysis = swimming_DEM_analysis.SwimmingDEMAnalysis

class MarineRainAnalysis(BaseAnalysis):
    def __init__(self, varying_parameters = Parameters("{}")):
        BaseAnalysis.__init__(self, varying_parameters)
        self.pp.CFD_DEM.do_search_neighbours = False
        self.pp.CFD_DEM.vorticity_calculation_type = 0
        self.pp.CFD_DEM.do_impose_flow_from_field = True
        self.pp.CFD_DEM.print_MATERIAL_ACCELERATION_option = False

        self.pp.CFD_DEM.L = 0.1
        self.pp.CFD_DEM.U = 0.3
        self.pp.CFD_DEM.k = 2.72
        self.pp.CFD_DEM.omega = math.pi

    def GetFieldUtility(self):
        self.flow_field = CellularFlowField(self.pp.CFD_DEM.L, self.pp.CFD_DEM.U, self.pp.CFD_DEM.k, self.pp.CFD_DEM.omega)
        space_time_set = SpaceTimeSet()
        self.field_utility = FluidFieldUtility(space_time_set, self.flow_field, 1000.0, 1e-6)
        return self.field_utility

    def PerformZeroStepInitializations(self):
        import random
        #from random import randint
        L = self.pp.CFD_DEM.L
        #U = self.pp.CFD_DEM.U
        #k = self.pp.CFD_DEM.k
        #omega = self.pp.CFD_DEM.omega
        N_positions = 100
        L *= math.pi
        possible_xs = [2 * L * (i + 1) / N_positions for i in range(N_positions)]
        i_position = 0

        for node in self.spheres_model_part.Nodes:
            rand_x = 2 * L * random.random()
            rand_y = self.pp.CFD_DEM.BoundingBoxMinY + (self.pp.CFD_DEM.BoundingBoxMaxY - self.pp.CFD_DEM.BoundingBoxMinY) * random.random()
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
            node.SetSolutionStepValue(DISPLACEMENT_X, rand_x - init_x)
            node.SetSolutionStepValue(DISPLACEMENT_Y, rand_y - init_y)
            i_position += 1

        for node in self.spheres_model_part.Nodes:
            vel= Vector(3)
            coor= Vector(3)
            coor[0]=node.X
            coor[1]=node.Y
            coor[2]=node.Z
            self.flow_field.Evaluate(0.0,coor,vel,0)
            node.SetSolutionStepValue(VELOCITY_X, vel[0])
            node.SetSolutionStepValue(VELOCITY_Y, vel[1])
            node.SetSolutionStepValue(VELOCITY_Z, vel[2])
            node.SetSolutionStepValue(VELOCITY_OLD_X, vel[0])
            node.SetSolutionStepValue(VELOCITY_OLD_Y, vel[1])
            node.SetSolutionStepValue(VELOCITY_OLD_Z, vel[2])
            node.SetSolutionStepValue(SLIP_VELOCITY_X, vel[0])
            node.SetSolutionStepValue(SLIP_VELOCITY_Y, vel[1])
            node.SetSolutionStepValue(SLIP_VELOCITY_Z, vel[2])

    def FluidSolve(self, time = 'None', solve_system=True):
        pass
