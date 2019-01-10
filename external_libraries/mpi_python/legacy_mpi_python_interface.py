from __future__ import print_function, absolute_import, division

import KratosMultiphysics
import KratosMultiphysics.mpi as mpi

class LegacyMPIPythonInterface(object):

    def __init__(self):
        self.world = KratosMultiphysics.ParallelEnvironment.GetDataCommunicator("World")

