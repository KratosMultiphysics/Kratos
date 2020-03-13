# Import of utilities
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.pyplot import xticks, yticks
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
import sys
import os  # to be able to run Linux terminal commands
from matplotlib import cm
import warnings
import pickle
import matplotlib.animation as animation

# Setting of plot parameters
matplotlib.rcParams['font.size'] = 10
matplotlib.rcParams['axes.linewidth'] = 1
matplotlib.rcParams['lines.linewidth'] = 2
matplotlib.rcParams['axes.labelweight'] = 'bold'
matplotlib.rcParams['xtick.major.size'] = 6
matplotlib.rcParams['xtick.minor.size'] = 2.4
matplotlib.rcParams['ytick.major.size'] = 6
matplotlib.rcParams['ytick.minor.size'] = 2.4
matplotlib.rcParams['savefig.dpi'] = 300

# Test_examples to plot
plot_04_tube = False
plot_04_tube_inert = False
plot_04_tube2D = False
plot_04_tube2D_inert = False
plot_tube2D_pipe_flow_abaqus = False
plot_tube2D_fluent_pipe_structure = False
plot_tube2D_fluent_abaqus = True
plot_tube2D_tango = True
plot_tube3D_fluent_abaqus = True
plot_tube3D_tango = True
plot_tube3D_tango_thick = True

# Initialization
fig, ax = plt.subplots()
def file(i):
    return f"/CSM/Area_TS{i}"
path = "../"

# Plot definitions
if plot_04_tube:
    name = "04_tube"
    tmp = path + name + file(0)
    A = np.array(pd.read_csv(tmp, sep="\s+"))
    x = x_04_tube = A[:, 0]
    line_04_tube, = ax.plot(x_04_tube, A[:, 1], label=name)

    def init_04_tube():  # only required for blitting to give a clean slate.
        line_04_tube.set_ydata([np.nan] * len(x_04_tube))
        return line_04_tube,

    def animate_04_tube(i):
        tmp = path + "04_tube" + file(i)
        A = np.array(pd.read_csv(tmp, sep="\s+"))
        line_04_tube.set_ydata(A[:, 1])  # update the data.
        return line_04_tube,

if plot_04_tube_inert:
    name = "04_tube_inert"
    tmp = path + name + file(0)
    A = np.array(pd.read_csv(tmp, sep="\s+"))
    x = x_04_tube_inert = A[:, 0]
    line_04_tube_inert, = ax.plot(x_04_tube_inert, A[:, 1], label=name)

    def init_04_tube_inert():  # only required for blitting to give a clean slate.
        line_04_tube_inert.set_ydata([np.nan] * len(x_04_tube_inert))
        return line_04_tube_inert,

    def animate_04_tube_inert(i):
        tmp = path + "04_tube_inert" + file(i)
        A = np.array(pd.read_csv(tmp, sep="\s+"))
        line_04_tube_inert.set_ydata(A[:, 1])  # update the data.
        return line_04_tube_inert,

if plot_04_tube2D:
    name = "04_tube2D"
    tmp = path + name + file(0)
    A = np.array(pd.read_csv(tmp, sep="\s+"))
    x = x_04_tube2D = A[:, 0]
    line_04_tube2D, = ax.plot(x_04_tube2D, A[:, 1], label="1D PipeFlow + PipeStructure (CoCoNuT)", color="tab:blue", linestyle="--")

    def init_04_tube2D():  # only required for blitting to give a clean slate.
        line_04_tube2D.set_ydata([np.nan] * len(x_04_tube2D))
        return line_04_tube2D,

    def animate_04_tube2D(i):
        tmp = path + "04_tube2D" + file(i)
        A = np.array(pd.read_csv(tmp, sep="\s+"))
        line_04_tube2D.set_ydata(A[:, 1])  # update the data.
        return line_04_tube2D,

if plot_04_tube2D_inert:
    name = "04_tube2D_inert"
    tmp = path + name + file(0)
    A = np.array(pd.read_csv(tmp, sep="\s+"))
    x = x_04_tube2D_inert = A[:, 0]
    line_04_tube2D_inert, = ax.plot(x_04_tube2D_inert, A[:, 1], label="1D PipeFlow + Inert PipeStructure (CoCoNuT)", color="tab:blue", linestyle="-")

    def init_04_tube2D_inert():  # only required for blitting to give a clean slate.
        line_04_tube2D_inert.set_ydata([np.nan] * len(x_04_tube2D_inert))
        return line_04_tube2D_inert,

    def animate_04_tube2D_inert(i):
        tmp = path + "04_tube2D_inert" + file(i)
        A = np.array(pd.read_csv(tmp, sep="\s+"))
        line_04_tube2D_inert.set_ydata(A[:, 1])  # update the data.
        return line_04_tube2D_inert,

if plot_tube2D_pipe_flow_abaqus:
    def file_tube2D_pipe_flow_abaqus(i):
        return f"/CSM/CSM_Time{i}Surface0Output.dat"
    name = "tube2D_pipe_flow_abaqus"
    tmp = path + name + f"/CSM/CSM_Time0Surface0Nodes.dat"
    A = np.array(pd.read_csv(tmp, sep="\s+"))
    x0 = A[:, 1]
    y0 = A[:, 0]
    sort_tube2D_pipe_flow_abaqus = np.argsort(x0)
    x = x0_tube2D_pipe_flow_abaqus = x0[sort_tube2D_pipe_flow_abaqus]
    y0_tube2D_pipe_flow_abaqus = y0[sort_tube2D_pipe_flow_abaqus]
    line_tube2D_pipe_flow_abaqus, = ax.plot(x0_tube2D_pipe_flow_abaqus, y0_tube2D_pipe_flow_abaqus ** 2 * np.pi, label="2D PipeFlow + Abaqus (CoCoNuT)", color="tab:orange", linestyle="--")

    def init_tube2D_pipe_flow_abaqus():  # only required for blitting to give a clean slate.
        line_tube2D_pipe_flow_abaqus.set_ydata([np.nan] * len(x0_tube2D_pipe_flow_abaqus))
        return line_tube2D_pipe_flow_abaqus,

    def animate_tube2D_pipe_flow_abaqus(i):
        tmp = path + "tube2D_pipe_flow_abaqus" + file_tube2D_pipe_flow_abaqus(i)
        A = np.array(pd.read_csv(tmp, sep="\s+"))[sort_tube2D_pipe_flow_abaqus]
        line_tube2D_pipe_flow_abaqus.set_ydata((A[:, 0] + y0_tube2D_pipe_flow_abaqus) ** 2 * np.pi)  # update the data.
        line_tube2D_pipe_flow_abaqus.set_xdata(A[:, 1] + x0_tube2D_pipe_flow_abaqus)  # update the data.
        return line_tube2D_pipe_flow_abaqus,

if plot_tube2D_fluent_pipe_structure:
    def file_fluent_pipe_structure(i):
        return f"/CSM/Area_TS{i}"
    name = "tube2D_fluent_pipe_structure"
    tmp = path + name + file(0)
    A = np.array(pd.read_csv(tmp, sep="\s+"))
    x = x_tube2D_fluent_pipe_structure = A[:, 0]
    line_tube2D_fluent_pipe_structure, = ax.plot(x_tube2D_fluent_pipe_structure, A[:, 1], label="2D Fluent + Inert PipeStructure (CoCoNuT)", color="tab:orange", linestyle=":")

    def init_tube2D_fluent_pipe_structure():  # only required for blitting to give a clean slate.
        line_tube2D_fluent_pipe_structure.set_ydata([np.nan] * len(x_tube2D_fluent_pipe_structure))
        return line_tube2D_fluent_pipe_structure,

    def animate_tube2D_fluent_pipe_structure(i):
        tmp = path + "tube2D_fluent_pipe_structure" + file_fluent_pipe_structure(i)
        A = np.array(pd.read_csv(tmp, sep="\s+"))
        line_tube2D_fluent_pipe_structure.set_ydata(A[:, 1])  # update the data.
        return line_tube2D_fluent_pipe_structure,

if plot_tube2D_fluent_abaqus:
    # 2D Fluent Abaqus CoCoNuT
    def file_tube2D_fluent_abaqus(i):
        return f"/CSM/CSM_Time{i}Surface0Output.dat"
    name = "tube2D_fluent_abaqus"
    tmp = path + name + f"/CSM/CSM_Time0Surface0Nodes.dat"
    A = np.array(pd.read_csv(tmp, sep="\s+"))
    x = x0 = A[:, 1]
    y0 = A[:, 0]
    sort_tube2D_fluent_abaqus = np.argsort(x0)
    x0_tube2D_fluent_abaqus = x0[sort_tube2D_fluent_abaqus]
    y0_tube2D_fluent_abaqus = y0[sort_tube2D_fluent_abaqus]
    line_tube2D_fluent_abaqus, = ax.plot(x0_tube2D_fluent_abaqus, y0_tube2D_fluent_abaqus ** 2 * np.pi, label="2D Fluent + Abaqus (CoCoNuT)", color="tab:orange", linestyle="-")


    def init_tube2D_fluent_abaqus():  # only required for blitting to give a clean slate.
        line_tube2D_fluent_abaqus.set_ydata([np.nan] * len(x0_tube2D_fluent_abaqus))
        return line_tube2D_fluent_abaqus,

    def animate_tube2D_fluent_abaqus(i):
        tmp = path + "tube2D_fluent_abaqus" + file_tube2D_fluent_abaqus(i)
        A = np.array(pd.read_csv(tmp, sep="\s+"))[sort_tube2D_fluent_abaqus]
        line_tube2D_fluent_abaqus.set_ydata(
            (A[:, 0] + y0_tube2D_fluent_abaqus) ** 2 * np.pi)  # update the data.
        line_tube2D_fluent_abaqus.set_xdata(A[:, 1] + x0_tube2D_fluent_abaqus)  # update the data.
        return line_tube2D_fluent_abaqus,

if plot_tube2D_tango:
    # 2D Fluent Abaqus Tango
    dir_Tango = '/cfdfile2/data/fm/nicolas/Tango/cases/Data_Tube2D/FluentSolver0/'
    tmp = dir_Tango + f'FSI1Refine0Time0Surface0NodesOrdered.dat'
    A_tango = np.array(pd.read_csv(tmp, header=None, skiprows=1, sep="\s+"))
    x_tango = A_tango[:, 1]
    line_tango, = ax.plot(x_tango, A_tango[:, 2] ** 2 * np.pi, label="2D Fluent + Abaqus (Tango)", color="tab:orange", marker='o', ms=2, linestyle="")

    def init_tube2D_tango():  # only required for blitting to give a clean slate.
        line_tango.set_ydata([np.nan] * len(x_tango))
        return line_tango,

    def animate_tube2D_tango(i):
        # print(i)
        tmp = dir_Tango + f'FSI1Refine0Time{i}Surface0NodesOrdered.dat'
        A_tango = np.array(pd.read_csv(tmp, header=None, skiprows=1, sep="\s+"))
        line_tango.set_xdata(A_tango[:, 0])
        line_tango.set_ydata(A_tango[:, 1] ** 2 * np.pi)  # update the data.
        return line_tango,

if plot_tube3D_fluent_abaqus:
    # 3D Fluent Abaqus CoCoNuT
    def file_tube3D_fluent_abaqus(i):
        return f"/CSM/CSM_Time{i}Surface0Output.dat"
    name = "tube3D_fluent_abaqus_finer"
    tmp = path + name + f"/CSM/CSM_Time0Surface0Nodes.dat"
    A = np.array(pd.read_csv(tmp, sep="\s+"))
    mask_3dcoco = (abs(A[:, 2]) < 1e-16) & (A[:, 1] > 0)
    A = A[mask_3dcoco]
    sort_tube3D_fluent_abaqus = np.argsort(A[:, 0])
    pos0_3dcoco = A[sort_tube3D_fluent_abaqus, :]
    x0_tube3D_fluent_abaqus = pos0_3dcoco[:, 0]
    y0_tube3D_fluent_abaqus = pos0_3dcoco[:, 1]
    line_tube3D_fluent_abaqus, = ax.plot(x0_tube3D_fluent_abaqus, y0_tube3D_fluent_abaqus ** 2 * np.pi, label="3D Fluent + Abaqus (CoCoNuT)", color="tab:green", linestyle="-")

    def init_tube3D_fluent_abaqus():  # only required for blitting to give a clean slate.
        line_tube3D_fluent_abaqus.set_ydata([np.nan] * len(x0_tube3D_fluent_abaqus))
        return line_tube3D_fluent_abaqus,

    def animate_tube3D_fluent_abaqus(i):
        tmp = path + "tube3D_fluent_abaqus_finer" + file_tube3D_fluent_abaqus(i)
        A = np.array(pd.read_csv(tmp, sep="\s+"))
        A = A[mask_3dcoco]
        A = A[sort_tube3D_fluent_abaqus, :] + pos0_3dcoco
        line_tube3D_fluent_abaqus.set_xdata(A[:, 0])  # update the data.
        line_tube3D_fluent_abaqus.set_ydata(A[:, 1] ** 2 * np.pi)  # update the data.
        return line_tube3D_fluent_abaqus,

if plot_tube3D_tango:
    # 3D Fluent Abaqus Tango
    dir_3dtango = '/cfdfile2/data/fm/nicolas/Tango/cases/Data_Tube3D_finer/AbaqusSolver0/'
    tmp = dir_3dtango + f'FSI1Refine0Time0Surface0Nodes.dat'
    A_3dtango = np.array(pd.read_csv(tmp, header=None, skiprows=1, sep="\s+"))
    mask_3dtango = (abs(A_3dtango[:, 2]) < 1e-16) & (A_3dtango[:, 1] > 0)
    A_3dtango = A_3dtango[mask_3dtango]
    sort_id_3dtango = np.argsort(A_3dtango[:, 0])
    pos0_3dtango = A_3dtango[sort_id_3dtango, :]
    x_3dtango = pos0_3dtango[:, 0]
    line_3dtango, = ax.plot(x_3dtango, A_3dtango[:, 2] ** 2 * np.pi, label="3D Fluent + Abaqus (Tango)", color="tab:green", marker='o', ms=2, linestyle="")

    def init_tube3D_tango():  # only required for blitting to give a clean slate.
        line_3dtango.set_ydata([np.nan] * len(x_3dtango))
        return line_3dtango,

    def animate_tube3D_tango(i):
        tmp = dir_3dtango + f'FSI1Refine0Time{i}Surface0Output.dat'
        A_3dtango = np.array(pd.read_csv(tmp, header=None, skiprows=1, sep="\s+"))
        A_3dtango = A_3dtango[mask_3dtango]
        A_3dtango = A_3dtango[sort_id_3dtango, :] + pos0_3dtango
        line_3dtango.set_xdata(A_3dtango[:, 0])
        line_3dtango.set_ydata(A_3dtango[:, 1] ** 2 * np.pi)  # update the data.
        return line_3dtango,

if plot_tube3D_tango_thick:
    # 3D Fluent Abaqus Tango
    dir_3dtango_thick = '/cfdfile2/data/fm/nicolas/Tango/cases/Data/AbaqusSolver0/'
    tmp = dir_3dtango_thick + f'FSI1Refine0Time0Surface0Nodes.dat'
    A_3dtango_thick = np.array(pd.read_csv(tmp, header=None, skiprows=1, sep="\s+"))
    mask_3dtango_thick = (abs(A_3dtango_thick[:, 2]) < 1e-16) & (A_3dtango_thick[:, 1] > 0)
    A_3dtango_thick = A_3dtango_thick[mask_3dtango_thick]
    sort_id_3dtango_thick = np.argsort(A_3dtango_thick[:, 0])
    pos0_3dtango_thick = A_3dtango_thick[sort_id_3dtango_thick, :]
    x_3dtango_thick = pos0_3dtango_thick[:, 0]
    line_3dtango_thick, = ax.plot(x_3dtango_thick, A_3dtango_thick[:, 2] ** 2 * np.pi, label="3D Fluent + Abaqus thick(Tango)", color="tab:green", marker='o', ms=2, linestyle="")

    def init_tube3D_tango_thick():  # only required for blitting to give a clean slate.
        line_3dtango_thick.set_ydata([np.nan] * len(x_3dtango_thick))
        return line_3dtango_thick,

    def animate_tube3D_tango_thick(i):
        tmp = dir_3dtango_thick + f'FSI1Refine0Time{i}Surface0Output.dat'
        A_3dtango_thick = np.array(pd.read_csv(tmp, header=None, skiprows=1, sep="\s+"))
        A_3dtango_thick = A_3dtango_thick[mask_3dtango_thick]
        A_3dtango_thick = A_3dtango_thick[sort_id_3dtango_thick, :] + pos0_3dtango_thick
        line_3dtango_thick.set_xdata(A_3dtango_thick[:, 0])
        line_3dtango_thick.set_ydata(A_3dtango_thick[:, 1] ** 2 * np.pi)  # update the data.
        return line_3dtango_thick,

# Setting scale
a0 = 0.005 ** 2 * np.pi #A[0, 1]
amp = 5e-06#3e-07
plt.ylim([a0 - amp, a0 + amp])

# Animation parameters
interv = 100  # Interval between frames in ms
T = range(1, 21)  # Number of frames
svc = 100

# # Animations
# if plot_04_tube:
#     ani_04_tube = animation.FuncAnimation(fig, animate_04_tube, init_func=init_04_tube, interval=interv, blit=False, save_count=svc, repeat=True, frames=T)
# if plot_04_tube_inert:
#     ani_04_tube_inert = animation.FuncAnimation(fig, animate_04_tube_inert, init_func=init_04_tube_inert, interval=interv, blit=False, save_count=svc, repeat=True, frames=T)
# if plot_04_tube2D:
#     ani_04_tube2D = animation.FuncAnimation(fig, animate_04_tube2D, init_func=init_04_tube2D, interval=interv, blit=False, save_count=svc, repeat=True, frames=T)
# if plot_04_tube2D_inert:
#     ani_04_tube2D_inert = animation.FuncAnimation(fig, animate_04_tube2D_inert, init_func=init_04_tube2D_inert, interval=interv, blit=False, save_count=svc, repeat=True, frames=T)
# if plot_tube2D_pipe_flow_abaqus:
#     ani_tube2D_pipe_flow_abaqus = animation.FuncAnimation(fig, animate_tube2D_pipe_flow_abaqus, init_func=init_tube2D_pipe_flow_abaqus, interval=interv, blit=False, save_count=svc, repeat=True, frames=T)
# if plot_tube2D_fluent_pipe_structure:
#     ani_tube2D_fluent_pipe_structure = animation.FuncAnimation(fig, animate_tube2D_fluent_pipe_structure, init_func=init_tube2D_fluent_pipe_structure, interval=interv, blit=False, save_count=svc, repeat=True, frames=T)
# if plot_tube2D_fluent_abaqus:
#     ani_tube2D_fluent_abaqus = animation.FuncAnimation(fig, animate_tube2D_fluent_abaqus, init_func=init_tube2D_fluent_abaqus, interval=interv, blit=False, save_count=svc, repeat=True, frames=T)
# if plot_tube2D_tango:
#     ani_tube2D_tango = animation.FuncAnimation(fig, animate_tube2D_tango, init_func=init_tube2D_tango, interval=interv, blit=False, save_count=svc, repeat=True, frames=T)
# if plot_tube3D_fluent_abaqus:
#     ani_tube3D_fluent_abaqus = animation.FuncAnimation(fig, animate_tube3D_fluent_abaqus, init_func=init_tube3D_fluent_abaqus, interval=interv, blit=False, save_count=svc, repeat=True, frames=T)
# if plot_tube3D_tango:
#     ani_tube3D_tango = animation.FuncAnimation(fig, animate_tube3D_tango, init_func=init_tube3D_tango, interval=interv, blit=False, save_count=svc, repeat=True, frames=T)


# Animations
def init():
    lines = ()
    if plot_04_tube:
        line = init_04_tube()
        lines += line
    if plot_04_tube_inert:
        line = init_04_tube_inert()
        lines += line
    if plot_04_tube2D:
        line = init_04_tube2D()
        lines += line
    if plot_04_tube2D_inert:
        line = init_04_tube2D_inert()
        lines += line
    if plot_tube2D_pipe_flow_abaqus:
        line = init_tube2D_pipe_flow_abaqus()
        lines += line
    if plot_tube2D_fluent_pipe_structure:
        line = init_tube2D_fluent_pipe_structure()
        lines += line
    if plot_tube2D_fluent_abaqus:
        line = init_tube2D_fluent_abaqus()
        lines += line
    if plot_tube2D_tango:
        line = init_tube2D_tango()
        lines += line
    if plot_tube3D_fluent_abaqus:
        line = init_tube3D_fluent_abaqus()
        lines += line
    if plot_tube3D_tango:
        line = init_tube3D_tango()
        lines += line
    if plot_tube3D_tango_thick:
        line = init_tube3D_tango_thick()
        lines += line
    return lines,


def animate(i):
    lines = ()
    if plot_04_tube:
        line = animate_04_tube(i)
        lines += line
    if plot_04_tube_inert:
        line = animate_04_tube_inert(i)
        lines += line
    if plot_04_tube2D:
        line = animate_04_tube2D(i)
        lines += line
    if plot_04_tube2D_inert:
        line = animate_04_tube2D_inert(i)
        lines += line
    if plot_tube2D_pipe_flow_abaqus:
        line = animate_tube2D_pipe_flow_abaqus(i)
        lines += line
    if plot_tube2D_fluent_pipe_structure:
        line = animate_tube2D_fluent_pipe_structure(i)
        lines += line
    if plot_tube2D_fluent_abaqus:
        line = animate_tube2D_fluent_abaqus(i)
        lines += line
    if plot_tube2D_tango:
        line = animate_tube2D_tango(i)
        lines += line
    if plot_tube3D_fluent_abaqus:
        line = animate_tube3D_fluent_abaqus(i)
        lines += line
    if plot_tube3D_tango:
        line = animate_tube3D_tango(i)
        lines += line
    if plot_tube3D_tango_thick:
        line = animate_tube3D_tango_thick(i)
        lines += line
    return lines,

ani = animation.FuncAnimation(fig, animate, init_func=init, interval=interv, blit=False, save_count=svc, repeat=True, frames=T)


plt.plot(x, a0 * np.ones_like(x), color='k')
plt.title("Animation of the area of each crosssection")
plt.ylabel("area (m²)")
plt.xlabel("z-coordinate (m)")
plt.legend(loc='lower center')

# Set up formatting for the movie files
### To get an mp4-file
print(plt.rcParams['animation.ffmpeg_path'])
plt.rcParams['animation.ffmpeg_path'] = u'/apps/SL6.3/FFmpeg/3.1.3/bin/ffmpeg'
print(plt.rcParams['animation.ffmpeg_path'])
writer = animation.FFMpegFileWriter(codec='mpeg1video', metadata=dict(artist='NicolasDelaissé'), fps=24, bitrate=2000)

# if bool_save:
#     line_ani.save(save_name, writer=writer)
#
# Writer = animation.writers['ffmpeg']
# writer = Writer(fps=15, metadata=dict(artist='NicolasDelaissé'), bitrate=1800)

# ani.save('finer2.mp4', writer=writer)

plt.show()

