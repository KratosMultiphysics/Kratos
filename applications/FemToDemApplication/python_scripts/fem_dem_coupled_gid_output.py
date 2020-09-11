from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
import KratosMultiphysics as Kratos
import KratosMultiphysics.FemToDemApplication as FEMDEM
from KratosMultiphysics import MultiFileFlag
from KratosMultiphysics import GiDPostMode
from KratosMultiphysics import Logger
from KratosMultiphysics import DEMApplication
from KratosMultiphysics import DemStructuresCouplingApplication
from KratosMultiphysics import gid_output # this is deprecated

class FemDemCoupledGiDOutput(gid_output.GiDOutput):

    """The constructor of the FemDemCoupledGiDOutput-Object.
    This object prints the results from the PFEM, FEMDEM and DEM 
    apps in the same post file 
    """
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
        self.fluid_model_part                   = fluid_model_part

        self.balls_model_part                   = balls_model_part
        self.clusters_model_part                = clusters_model_part
        self.rigid_faces_model_part             = rigid_faces_model_part

        self.mixed_solid_fluid_model_part       = mixed_solid_fluid_model_part
        self.mixed_solid_balls_model_part       = mixed_solid_balls_model_part
        self.mixed_solid_balls_fluid_model_part = mixed_solid_balls_fluid_model_part

        self.solid_nodal_results                   = []
        self.fluid_nodal_results                   = []
        self.dem_nodal_results                     = []
        self.clusters_nodal_results                = []
        self.rigid_faces_nodal_results             = []

        self.mixed_solid_fluid_nodal_results       = []
        self.mixed_solid_balls_nodal_results       = []
        self.mixed_solid_balls_fluid_nodal_results = []

        self.solid_gauss_points_results            = []


    """ initialize_dem_fem_results
    This method initializes the post coupled post process
    """
    def initialize_dem_fem_results(self,
                                   solid_nodal_results,
                                   fluid_nodal_results,
                                   dem_nodal_results,
                                   clusters_nodal_results,
                                   rigid_faces_nodal_results,
                                   mixed_solid_fluid_nodal_results,
                                   mixed_solid_balls_nodal_results,
                                   mixed_solid_balls_fluid_nodal_results,
                                   solid_gauss_points_results):

        self.solid_nodal_results       = solid_nodal_results
        self.fluid_nodal_results       = fluid_nodal_results
        self.dem_nodal_results         = dem_nodal_results
        self.clusters_nodal_results    = clusters_nodal_results
        self.rigid_faces_nodal_results = rigid_faces_nodal_results

        self.mixed_solid_fluid_nodal_results       = mixed_solid_fluid_nodal_results
        self.mixed_solid_balls_nodal_results       = mixed_solid_balls_nodal_results
        self.mixed_solid_balls_fluid_nodal_results = mixed_solid_balls_fluid_nodal_results

        self.solid_gauss_points_results = solid_gauss_points_results

        if self.multi_file == MultiFileFlag.SingleFile:
            print("Singlefile option is not available for the DEM-Structures Coupling application!")
            mesh_name = 0.0
            self.io.InitializeMesh(mesh_name)
            self.io.WriteSphereMesh(self.balls_model_part.GetMesh())
            self.io.WriteMesh(self.clusters_model_part.GetMesh())
            self.io.WriteMesh(self.rigid_faces_model_part.GetMesh())
            # self.io.WriteMesh(self.mixed_solid_fluid_model_part.GetMesh())
            # self.io.WriteMesh(self.mixed_solid_balls_model_part.GetMesh())
            # self.io.WriteMesh(self.mixed_solid_balls_fluid_model_part.GetMesh())
            self.io.WriteMesh(self.fluid_model_part.GetMesh())
            self.io.WriteMesh(self.solid_model_part.GetMesh())
            self.io.WriteMesh(self.mixed_solid_balls_fluid_model_part.GetMesh())
            self.io.FinalizeMesh()
            self.io.InitializeResults(mesh_name, self.mixed_solid_balls_model_part.GetMesh())

        # Initialize list file
        with open(self.listfilename, "w") as listfile:
            if self.multi_file == MultiFileFlag.MultipleFiles:
                listfile.write("Multiple\n")
            elif self.multi_file == MultiFileFlag.SingleFile:
                listfile.write("Single\n")
        if self.multi_file == MultiFileFlag.SingleFile:
            if self.post_mode == GiDPostMode.GiD_PostBinary:
                self.write_step_to_list()
            else:
                self.write_step_to_list(0)


    """ write_step_to_list
    """
    def write_step_to_list(self, step_label):
    
        if self.post_mode   == GiDPostMode.GiD_PostBinary:
            ext = ".post.bin"
        elif self.post_mode == GiDPostMode.GiD_PostAscii:
            ext = ".post.res"
        elif self.post_mode == GiDPostMode.GiD_PostAsciiZipped:
            ext = ".post.res"

        with open(self.listfilename, "a") as listfile:
            listfile.write("Multiple\n")
            listfile.write(self.filename+"_"+"%.12g"%step_label+ext+"\n")


    """ write_step_to_outer_list
    """
    def write_step_to_outer_list(self, step_label):

        if self.post_mode   == GiDPostMode.GiD_PostBinary:
            ext = ".post.bin"
        elif self.post_mode == GiDPostMode.GiD_PostAscii:
            ext = ".post.res"
        elif self.post_mode == GiDPostMode.GiD_PostAsciiZipped:
            ext = ".post.res"

        with open(self.outerlistfilename, "a") as listfile:
            listfile.write("Multiple\n")
            folder_name = self.filename + "_Post_Files"
            full_string_to_write = os.path.join(folder_name, self.filename + "_" + "%.12g"%step_label + ext)
            listfile.write(full_string_to_write + "\n")


    """ Writeresults
    """
    def Writeresults(self, time):

        # We reorder the Id of the model parts
        femdem_util = FEMDEM.FEMDEMCouplingUtilities()
        reorder_util_elem = FEMDEM.RenumberingNodesUtility(self.solid_model_part, self.fluid_model_part)
        reorder_util_elem.RenumberElements()


        Logger.PrintInfo("","")
        Logger.PrintInfo("","*****************  PRINTING RESULTS FOR GID  *************************")
        Logger.PrintInfo("","")
        Logger.Flush()

        number_pfem_nodes = femdem_util.GetNumberOfNodes(self.fluid_model_part)
        for node in self.balls_model_part.Nodes:
            node.Id = node.Id + number_pfem_nodes

        if self.GiDMultiFileFlag == "Multiples":
            self.mixed_solid_fluid_model_part.Elements.clear()
            self.mixed_solid_fluid_model_part.Nodes.clear()

            self.mixed_solid_balls_model_part.Elements.clear()
            self.mixed_solid_balls_model_part.Nodes.clear()

            self.mixed_solid_balls_fluid_model_part.Elements.clear()
            self.mixed_solid_balls_fluid_model_part.Nodes.clear()

            # Now we fill the mixed MDPA in order to print
            post_utils = DEMApplication.PostUtilities()
            post_utils.AddModelPartToModelPart(self.mixed_solid_fluid_model_part, self.fluid_model_part)

            post_utils.AddModelPartToModelPart(self.mixed_solid_balls_model_part, self.balls_model_part)
            post_utils.AddModelPartToModelPart(self.mixed_solid_balls_model_part, self.rigid_faces_model_part)

            post_utils.AddModelPartToModelPart(self.mixed_solid_balls_fluid_model_part, self.balls_model_part)
            post_utils.AddModelPartToModelPart(self.mixed_solid_balls_fluid_model_part, self.rigid_faces_model_part)
            post_utils.AddModelPartToModelPart(self.mixed_solid_balls_fluid_model_part, self.solid_model_part)
            post_utils.AddModelPartToModelPart(self.mixed_solid_balls_fluid_model_part, self.fluid_model_part)

            FEMDEM.FEMDEMCouplingUtilities().RemoveDuplicates(self.mixed_solid_fluid_model_part)

        self.write_dem_fem_results(time)
        reorder_util_elem.UndoRenumberElements()


    """ write_dem_fem_results
    """
    def write_dem_fem_results(self, label):

        if self.multi_file == MultiFileFlag.MultipleFiles:
            self.io.InitializeMesh(label)
            self.io.WriteSphereMesh(self.balls_model_part.GetMesh())
            self.io.WriteMesh(self.rigid_faces_model_part.GetMesh())
            self.io.WriteMesh(self.mixed_solid_balls_fluid_model_part.GetMesh())
            self.io.FinalizeMesh()
            self.io.InitializeResults(label, self.mixed_solid_balls_fluid_model_part.GetMesh())

        for var in  self.solid_nodal_results:
            kratos_variable = Kratos.KratosGlobals.GetVariable(var)
            self._write_nodal_results(label, self.solid_model_part, kratos_variable)

        for var in  self.fluid_nodal_results:
            kratos_variable = Kratos.KratosGlobals.GetVariable(var)
            self._write_nodal_results(label, self.fluid_model_part, kratos_variable)

        for var in self.dem_nodal_results:
            kratos_variable = Kratos.KratosGlobals.GetVariable(var)
            self._write_nodal_results(label, self.balls_model_part, kratos_variable)

        for var in self.clusters_nodal_results:
            kratos_variable = Kratos.KratosGlobals.GetVariable(var)
            self._write_nodal_results(label, self.clusters_model_part, kratos_variable)

        for var in self.rigid_faces_nodal_results:
            kratos_variable = Kratos.KratosGlobals.GetVariable(var)
            self._write_nodal_results(label, self.rigid_faces_model_part, kratos_variable)

        for var in self.mixed_solid_balls_fluid_nodal_results:
            kratos_variable = Kratos.KratosGlobals.GetVariable(var)
            self._write_nodal_results(label, self.mixed_solid_balls_fluid_model_part, kratos_variable)

        for var in self.mixed_solid_balls_nodal_results:
            kratos_variable = Kratos.KratosGlobals.GetVariable(var)
            self._write_nodal_results(label, self.mixed_solid_balls_model_part, kratos_variable)

        for var in self.mixed_solid_fluid_nodal_results:
            kratos_variable = Kratos.KratosGlobals.GetVariable(var)
            self._write_nodal_results(label, self.mixed_solid_fluid_model_part, kratos_variable)

        for var in self.solid_gauss_points_results:
            kratos_variable = Kratos.KratosGlobals.GetVariable(var)
            self._write_gp_results(label, self.solid_model_part, kratos_variable)

        if self.multi_file == MultiFileFlag.MultipleFiles:
            self._finalize_results()

            with open(self.listfilename, "a") as listfile:
                self.write_step_to_list(label)

            self.write_step_to_outer_list(label)


