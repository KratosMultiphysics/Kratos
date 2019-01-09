import KratosMultiphysics
import matplotlib.pyplot as plt
import numpy as np

class PlotManufacturedFieldUtility(object):
    def __init__(self, settings = KratosMultiphysics.Parameters("""{}""")):

        default_settings = KratosMultiphysics.Parameters("""{
            "benchmark_name"  : "polynomial_vortex",
            "Parameters"      : {}
        }""")

        settings.ValidateAndAssignDefaults(default_settings)

        # from polynomial_vortex import ManufacturedSolution
        benchmark_module = __import__(settings["benchmark_name"].GetString())
        self.benchmark = benchmark_module.CreateManufacturedSolution(settings["Parameters"])

        L = 1
        n_points = 20
        step = L / n_points
        x = np.arange(0, L, step)
        y = np.arange(0, L, step)
        self.X, self.Y = np.meshgrid(x, y)

    def PlotVelocity(self, time = 1.0):

        fig, ax = plt.subplots()
        q = ax.quiver(self.X, self.Y, self.benchmark.u1(self.X, self.Y, time), self.benchmark.u2(self.X, self.Y, time))
        ax.quiverkey(q, X=0.3, Y=1.1, U=1,label='Quiver key, length = 1', labelpos='E')
        plt.show()

    def PlotBodyForce(self, time = 1.0):

        fig, ax = plt.subplots()
        q = ax.quiver(self.X, self.Y, self.benchmark.body_force1(self.X, self.Y, time), self.benchmark.body_force2(self.X, self.Y, time))
        ax.quiverkey(q, X=0.3, Y=1.1, U=1,label='Quiver key, length = 1', labelpos='E')
        plt.show()

if __name__ == "__main__":
    PlotManufacturedFieldUtility().PlotVelocity()
