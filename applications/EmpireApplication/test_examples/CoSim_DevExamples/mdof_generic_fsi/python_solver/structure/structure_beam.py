from python_solver.structure.StructureMDoF import *
from python_solver.structure.StructuralProperties import *
from python_solver.element.beam import *
from python_solver.element.torsional_bar import *
from math import sin, cos, radians


class Structure():

    def __init__(self, _project_parameters):

        self.properties = StructuralProperties(_project_parameters)

        self.struct_X, self.struct_Y, self.struct_R = self.structure()
        self.solver_X, self.solver_Y, self.solver_R = self.solver()

        self.position = self.initial_position()
        self.old_results = self.predict_displacement()
        self.results = self.old_results

    def structure(self):

        # Restrained degrees of freedom / ex Cantelever
        restrained_dof_beam = [1, 0]
        restrained_dof_spring = [0]
        target_frequency_x = 0.23
        target_frequency_y = 0.2
        target_frequency_r = 0.4

        struct_X = Beam(self.properties,
                        restrained_dof_beam, target_frequency_x)
        struct_Y = Beam(self.properties,
                        restrained_dof_beam, target_frequency_y)
        struct_R = TorsionalBar(self.properties,
                                restrained_dof_spring, target_frequency_r)

        return [struct_X, struct_Y, struct_R]

    def solver(self):

        solver_X = StructureMDoF(
            self.properties.dt, self.struct_X.M, self.struct_X.K, self.struct_X.B, self.properties.rho_inf,
            self.properties.disp_X, self.properties.vel_X, self.properties.acc_X, self.properties.output_filename_X, self.properties.output_filename_Result + "X", self.struct_X.K_big, self.struct_X.M_big, self.struct_X.B_big)

        solver_Y = StructureMDoF(
            self.properties.dt, self.struct_Y.M, self.struct_Y.K, self.struct_Y.B, self.properties.rho_inf,
            self.properties.disp_Y, self.properties.vel_Y, self.properties.acc_Y, self.properties.output_filename_Y, self.properties.output_filename_Result + "Y", self.struct_Y.K_big, self.struct_Y.M_big, self.struct_Y.B_big)

        solver_R = StructureMDoF(
            self.properties.dt, self.struct_R.M, self.struct_R.K, self.struct_R.B, self.properties.rho_inf,
            self.properties.disp_R, self.properties.vel_R, self.properties.acc_R, self.properties.output_filename_R, self.properties.output_filename_Result + "R", self.struct_R.K_big, self.struct_R.M_big, self.struct_R.B_big)

        return [solver_X, solver_Y, solver_R]

    def predict_displacement(self):

        result_X = self.solver_X.predictDisplacement()
        result_Y = self.solver_Y.predictDisplacement()
        result_R = self.solver_R.predictDisplacement()

        disp_X, rot_X = result_X[::2], result_Y[1::2]
        disp_Y, rot_Y = result_Y[::2], result_X[1::2]
        disp_Z = np.zeros(len(result_R))

        return [disp_X, disp_Y, disp_Z, rot_X, rot_Y, result_R]

    def solve(self, forces):

        self.solver_X.solveStructure(forces[0])
        self.solver_Y.solveStructure(forces[1])
        self.solver_R.solveStructure(forces[2])

    def solve_static(self, forces):

        result_X = np.linalg.solve(self.struct_X.K, forces[0])
        result_Y = np.linalg.solve(self.struct_Y.K, forces[1])
        result_R = np.linalg.solve(self.struct_R.K, forces[2])

        disp_X, rot_X = result_X[::2], result_Y[1::2]
        disp_Y, rot_Y = result_Y[::2], result_X[1::2]
        disp_Z = np.zeros(len(result_R))

        results = disp_X, disp_Y, disp_Z, rot_X, rot_Y, result_R

        return results

    def get_displacement(self):

        result_X = self.solver_X.getDisplacement()
        result_Y = self.solver_Y.getDisplacement()
        result_R = self.solver_R.getDisplacement()

        disp_X, rot_X = result_X[::2], result_Y[1::2]
        disp_Y, rot_Y = result_Y[::2], result_X[1::2]
        disp_Z = np.zeros(len(result_R))

        self.results = disp_X, disp_Y, disp_Z, rot_X, rot_Y, result_R

    # def get_reaction(self):
    #     def get_values(solver, structure):

    #         filename_force = self.properties.output_filename_Result + str(structure) + "_force" + ".dat"
    #         filename_moment = self.properties.output_filename_Result + str(structure) + "_moment" + ".dat"

    #         support_output_force = open(filename_force, 'w')
    #         out = "#Results for group " + "\n"
    #         out += "#time    BASE FORCE REACTION \n"
    #         support_output_force.write(out)

    #         support_output_moment = open(filename_moment, 'w')
    #         out = "#Results for group " + "\n"
    #         out += "#time    BASE MOMENT REACTION \n"
    #         support_output_force.write(out)

    #         acc = solver.getAcceleration()
    #         vel = solver.getVelocity()
    #         dis = solver.getDisplacement()

    #         acc = np.insert(acc, 0, 0)
    #         if len(acc) != len(structure.M_big):
    #             acc = np.insert(acc, 0, 0)
    #         vel = np.insert(vel, 0, 0)
    #         if len(vel) != len(structure.B_big):
    #             vel = np.insert(vel, 0, 0)
    #         disp = np.insert(disp, 0, 0)
    #         if len(disp) != len(structure.K_big):
    #             disp = disp.insert(disp, 0, 0)

    #         reaction = np.dot(structure.M_big, acc) + np.dot(
    #             structure.B_big, vel) + np.dot(structure.K_big, disp)

    #         force, moment = reaction[::2], reaction[1::2]

    #         support_output_force.write(str(time) + " " + str(force[0]) + "\n")
    #         support_output_force.flush()

    #         support_output_moment.write(str(time) + " " + str(moment[0]) + "\n")
    #     self.support_output_moment.flush()

    def solve_eigen(self, mode, direction):
        mode = mode - 1
        if direction == "X":
            result_X = self.struct_X.eigen_value_load(mode)
            result_Y = np.zeros(len(result_X))
            result_R = np.zeros(len(self.struct_R.K))

        elif direction == "Y":
            result_Y = self.struct_Y.eigen_value_load(mode)
            result_X = np.zeros(len(result_Y))
            result_R = np.zeros(len(self.struct_R.K))

        elif direction == "R":
            result_R = self.struct_R.eigen_value_load(mode)
            result_X = np.zeros(len(self.struct_X.K))
            result_Y = np.zeros(len(self.struct_Y.K))

        disp_X, rot_X = result_X[::2], result_Y[1::2]
        disp_Y, rot_Y = result_Y[::2], result_X[1::2]
        disp_Z = np.zeros(len(result_R))

        self.results = disp_X, disp_Y, disp_Z, rot_X, rot_Y, result_R

    def update_relaxed_result(self, _new_result):

        self.results = _new_result
        self.old_results = _new_result

    def update_result(self):

        self.old_results = self.results

    def update_structure_time(self):

        updateX = self.solver_X.updateStructureTimeStep()
        updateY = self.solver_Y.updateStructureTimeStep()
        updateR = self.solver_R.updateStructureTimeStep()

        return updateX, updateY, updateR

    def print_support_output(self, time):

        outputX = self.solver_X.printSupportOutput(time)
        outputY = self.solver_Y.printSupportOutput(time)
        outputR = self.solver_R.printSupportOutput(time)

        return outputX, outputY, outputR

    def close_outpu(self):

        closeX = self.solver_X.support_output.close()
        closeY = self.solver_Y.support_output.close()
        closeR = self.solver_R.support_output.close()

        return closeX, closeY, closeR

    def initial_position(self):
        level_height = self.properties.height / self.properties.levels
        position = [[0, 0, ((i + 1) * level_height)]
                    for i in range(self.properties.levels)]

        return position

    def get_position(self):

        def transformation_matrix(result, level):

            alpha = result[2][level]
            beta = 0
            gamma = 0
            dispX = result[0][level]
            dispY = result[1][level]
            dispZ = 0

            # Transformation Matrix
            T = np.zeros((4, 4))
            T[0, 0] = cos(radians(alpha)) * cos(radians(beta))
            T[0, 1] = cos(radians(alpha)) * sin(radians(beta)) * sin(
                radians(gamma)) - sin(radians(alpha)) * cos(radians(gamma))
            T[0, 2] = cos(radians(alpha)) * sin(radians(beta)) * cos(
                radians(gamma)) + sin(radians(alpha)) * sin(radians(gamma))
            T[0, 3] = dispX
            T[1, 0] = sin(radians(alpha)) * cos(radians(beta))
            T[1, 1] = sin(radians(alpha)) * sin(radians(beta)) * sin(
                radians(gamma)) + cos(radians(alpha)) * cos(radians(gamma))
            T[1, 2] = sin(radians(alpha)) * sin(radians(beta)) * cos(
                radians(gamma)) - cos(radians(alpha)) * sin(radians(gamma))
            T[1, 3] = dispY
            T[2, 0] = -sin(radians(beta))
            T[2, 1] = cos(radians(beta)) * sin(radians(gamma))
            T[2, 2] = cos(radians(beta)) * cos(radians(gamma))
            T[2, 3] = dispZ
            T[3, 0] = 0
            T[3, 1] = 0
            T[3, 2] = 0
            T[3, 3] = 1

            return T

        for level in range(self.properties.levels):

            ref_position = [0, 0, 0]
            r_0 = [a - b for a, b in zip(self.position[level], ref_position)]
            r_0.append(1)
            T = transformation_matrix(self.results, level)
            r = np.dot(T, r_0)
            d = [a - b for a, b in zip(r, r_0)]
            print(d)

            self.position[level] = [
                a + b for a, b in zip(self.position[level], d)]

    def get_forces_back(self, time):

        fX = self.solver_X.getForcesBack(time)
        fY = self.solver_Y.getForcesBack(time)
        fR = self.solver_R.getForcesBack(time)

        return fX, fY, fR

