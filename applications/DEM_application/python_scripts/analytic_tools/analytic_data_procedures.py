import matplotlib.pyplot as plt

class WatcherAnalyzer:
    def __init__(self, analytic_face_watcher):
        self.face_watcher = analytic_face_watcher

    def MakeTotalFluxPlot(self):
        times = []
        fluxes = []
        self.face_watcher.GetTotalFlux(times, fluxes)
        self.accumulated_flux = [abs(sum(fluxes[:i])) for i in range(1, len(fluxes)+1)]
        plt.plot(times, self.accumulated_flux)
        plt.show()
