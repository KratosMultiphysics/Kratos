from numpy import linalg as la


class Convergence():

    def __init__(self, structure):

        self.abs_residual = structure.properties.fsi_abs_res
        self.relax_coef = structure.properties.fsi_relax_coef
        self.rel_residual = structure.properties.fsi_rel_res
        self.residual = None
        self.old_residual = None
        self.relaxed_solution = None

    def cal_residual(self, structure):

        self.old_residual = self.residual

        solution = structure.results
        old_solution = structure.old_results

        residual = [0] * len(solution)
        for index in range(len(solution)):
            residual[index] = solution[index] - old_solution[index]

        self.residual = residual

    def cal_relaxation(self, structure):

        print("RELAXATION FUNCTION CALLED:")
        print("OLD SOLUTION: ", structure.old_results[0])
        print("NEW SOLUTION: ", structure.results[0])
        old_solution = structure.results

        for i in range(len(old_solution)):
            for j in range(len(old_solution[i])):
                old_solution[i][j] = old_solution[i][
                    j] + self.relax_coef * self.residual[i][j]
        print("RELAXED SOLUTION")
        print(old_solution[0])

        self.relaxed_solution = old_solution

    def aitken_relaxation(self, iteration):

        max_initial_coeff = 0.125

        if iteration < 1:
            aitken_coef = min(max_initial_coeff, self.relax_coef)
        else:
            numerator = 0
            denominator = 0
            for i in range(len(self.residual)):
                for j in range(len(self.residual[i])):
                    numerator += self.old_residual[i][j] * (
                        self.residual[i][j] - self.old_residual[i][j])
                    denominator += pow(self.residual[i][j] - self.old_residual[i][j], 2)
            aitken_coef = - self.relax_coef * (numerator / denominator)

        self.relax_coef = aitken_coef
