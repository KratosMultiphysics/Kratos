from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
import KratosMultiphysics as Kratos
from KratosMultiphysics import MultiFileFlag
from KratosMultiphysics import GiDPostMode
from KratosMultiphysics import Logger
from KratosMultiphysics import DEMApplication
from KratosMultiphysics import DemStructuresCouplingApplication
from KratosMultiphysics import StructuralMechanicsApplication
from KratosMultiphysics import gid_output # this is deprecated

class FemDemCoupledGiDOutput(gid_output.GiDOutput):

    # ------------------------------------------------------------------------------------------- #
    def __init__(self,
                file_name,
                vol_output,
                post_mode,
                multifile,
                deformed_mesh,
                write_conditions,
                solid_model_part,
                fluid_model_part,
                balls_model_part,
                clusters_model_part,
                rigid_faces_model_part,
                mixed_solid_fluid_model_part,
                mixed_solid_balls_model_part,
                mixed_solid_balls_fluid_model_part):

        gid_output.GiDOutput.__init__(self,
                                  file_name,
                                  vol_output,
                                  post_mode,
                                  multifile,
                                  deformed_mesh,
                                  write_conditions)

        self.GiDMultiFileFlag = multifile
        self.outerlistfilename = os.path.join("..", self.listfilename)

        self.solid_model_part                   = solid_model_part
        self.balls_model_part                   = balls_model_part
        self.clusters_model_part                = clusters_model_part
        self.rigid_faces_model_part             = rigid_faces_model_part
        self.mixed_model_part                   = mixed_model_part
        self.mixed_solid_fluid_model_part       = mixed_solid_fluid_model_part
        self.mixed_solid_balls_model_part       = mixed_solid_balls_model_part
        self.mixed_solid_balls_fluid_model_part = mixed_solid_balls_fluid_model_part

        self.solid_nodal_results                   = []
        self.dem_nodal_results                     = []
        self.clusters_nodal_results                = []
        self.rigid_faces_nodal_results             = []
        self.mixed_solid_fluid_nodal_results       = []
        self.mixed_solid_balls_nodal_results       = []
        self.mixed_solid_balls_fluid_nodal_results = []
        self.structures_gauss_points_results       = []