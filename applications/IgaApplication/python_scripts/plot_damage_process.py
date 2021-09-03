import KratosMultiphysics as KM
import KratosMultiphysics.IgaApplication as IGA
import KratosMultiphysics.StructuralMechanicsApplication as SM

import matplotlib.pyplot as plt

from matplotlib.pylab import *
from mpl_toolkits.axes_grid1 import host_subplot
import matplotlib.animation as animation
import numpy as np

def Factory(settings, model):
    if(type(settings) != KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return PlotDamageProcess(model, settings["Parameters"])

## All the processes python should be derived from "Process"
class PlotDamageProcess(KM.Process):
    def __init__(self, model, settings ):
        KM.Process.__init__(self)

        default_settings = KM.Parameters("""{
            "model_part_name" : "please_specify_model_part_name",
            "on_nodes" : true,
            "on_elements" : true
        }""")

        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = model[settings["model_part_name"].GetString()]

        self.on_nodes = settings["on_nodes"].GetBool()
        self.on_elements = settings["on_elements"].GetBool()

        self.fig, (self.ax1, self.ax2) = plt.subplots(nrows=2)

    def ExecuteFinalizeSolutionStep(self):
        time_step = self.model_part.ProcessInfo[KM.TIME]

        npts = 2000
        ngridx = 1000
        ngridy = 2000


        # -----------------------
        # Interpolation on a grid
        # -----------------------
        # A contour plot of irregularly spaced data coordinates
        # via interpolation on a grid.

        # Create grid values first.
        xi = np.linspace(0, 10, ngridx)
        yi = np.linspace(0, 10, ngridy)

        x = []
        y = []
        z = []

        damage_t = []
        damage_c = []

        for element in self.model_part.Elements:
            location = element.GetGeometry().Center()

            x.append(location[0])
            y.append(location[1])
            z.append(location[2])

            damage_t.append(element.CalculateOnIntegrationPoints(SM.DAMAGE_TENSION, self.model_part.ProcessInfo)[0])
            damage_c.append(element.CalculateOnIntegrationPoints(SM.DAMAGE_COMPRESSION, self.model_part.ProcessInfo)[0])

        import matplotlib.tri as tri
        # Perform linear interpolation of the data (x,y)
        # on a grid defined by (xi,yi)
        triang = tri.Triangulation(x, y)
        interpolator = tri.LinearTriInterpolator(triang, damage_c)
        Xi, Yi = np.meshgrid(xi, yi)
        zi = interpolator(Xi, Yi)

        # Note that scipy.interpolate provides means to interpolate data on a grid
        # as well. The following would be an alternative to the four lines above:
        #from scipy.interpolate import griddata
        #zi = griddata((x, y), z, (xi[None,:], yi[:,None]), method='linear')

        self.ax1.clear()
        self.ax1.contour(xi, yi, zi, levels=35, linewidths=0.0, colors='k')
        cntr1 = self.ax1.contourf(xi, yi, zi, levels=35, cmap="RdBu_r")

        #self.fig.colorbar(cntr1, ax=self.ax1)
        self.ax1.plot(x, y, 'ko', ms=3)
        self.ax1.set(xlim=(0, 1), ylim=(0, 1))
        self.ax1.set_title('grid and contour (%d points, %d grid points)' %
                      (npts, ngridx * ngridy))


        # ----------
        # Tricontour
        # ----------
        # Directly supply the unordered, irregularly spaced coordinates
        # to tricontour.
        self.ax2.clear()
        self.ax2.tricontour(x, y, z, levels=35, linewidths=0.5, colors='k')
        cntr2 = self.ax2.tricontourf(x, y, damage_t, levels=35, cmap="RdBu_r")

        #self.fig.colorbar(cntr2, ax=self.ax2)
        self.ax2.plot(x, y, 'ko', ms=3)
        self.ax2.set(xlim=(0, 1), ylim=(0, 1))
        self.ax2.set_title('tricontour (%d points)' % npts)

        #plt.subplots_adjust(hspace=1)
        #plt.show()

        plt.draw()
        #plt.savefig('reactions_' + str(time_step) + '.png')
        plt.pause(1e-17)
