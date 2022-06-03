

# from math import min
import numpy as np
from pyevtk.hl import imageToVTK 
from imageio import imwrite
from .ErrorMessages import *


def exportVTK(FileName, cellData):
    shape = list(cellData.values())[0].shape
    ndim  = len(shape)
    N = min(shape)
    spacing = (1./N, 1./N, 1./N)

    if ndim==3:
        imageToVTK(FileName, cellData = cellData, spacing = spacing)

    elif ndim==2:
        cellData2D = {}
        for key in cellData.keys(): cellData2D[key] = np.expand_dims(cellData[key], axis=2)
        imageToVTK(FileName, cellData = cellData2D, spacing = spacing)

    else:
        msgDimError(ndim)

 

#=====================================================================================================
#=====================================================================================================
#
#                                           EXPORT
#
#=====================================================================================================
#=====================================================================================================


def save_png(phase, filename):
    imwrite(filename, (1-phase))

def save_vtk(phase, filename):
    exportVTK(filename, cellData = {'phase' : phase})