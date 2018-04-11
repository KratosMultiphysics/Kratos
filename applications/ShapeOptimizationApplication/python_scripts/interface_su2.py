# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
#
# ==============================================================================
# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Import Kratos core and apps
from KratosMultiphysics import *

# Import additional apps
import SU2
from itertools import islice
import shutil
from contextlib import contextmanager
import sys, os

# Functions
@contextmanager
def suppress_stdout():
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:
            yield
        finally:
            sys.stdout = old_stdout

# ==============================================================================
class InterfaceSU2():

    # --------------------------------------------------------------------------
    def __init__(self,interface_parameters ):
        default_parameters = Parameters("""
        {
            "su2_related":
            {
                "config_file"                 : "None",
                "number_of_cores"             : 1,
                "number_of_zones"             : 1,
                "gradient_method"             : "DISCRETE_ADJOINT",
                "write_frequency_adjoint_run" : 1000,
                "mesh_file"                   : "None",
                "mesh_motion_file"            : "mesh_motion.dat",
                "design_surface_tag"          : "None"
            },
            "kratos_related":
            {
                "mdpa_file"      : "converted_su2_mesh.mdpa",
                "write_elements" : false
            },
            "echo_level" : 1
        }""")
        interface_parameters.RecursivelyValidateAndAssignDefaults(default_parameters)
        self.interface_parameters = interface_parameters

        self.su2_line_id = 3
        self.su2_tria_id = 5
        self.su2_quad_id = 9
        self.su2_tet_id = 10
        self.su2_hex_id = 12
        self.su2_prism_id = 13
        self.su2_pyramid_id = 14

        self.su2_mesh_data  = {}
        self.node_id_su2_to_kratos = {}
        self.node_id_kratos_to_su2 = {}
        self.element_id_su2_to_kratos = {}
        self.list_of_optimization_conditions = []

        self.__ReadSU2Mesh()
        self.__TranslateSU2IdsToKratosIds()

    # --------------------------------------------------------------------------
    def WriteSU2MeshAsMDPA(self):
        if self.interface_parameters["echo_level"].GetInt()==1:
            print("> Start writing mdpa file...")

        self.new_file = open(self.interface_parameters["kratos_related"]["mdpa_file"].GetString(), 'w')

        self.__WriteHeader()
        self.__WriteNodes()
        if self.interface_parameters["kratos_related"]["write_elements"].GetBool():
            self.__WriteElements()
        if self.interface_parameters["su2_related"]["design_surface_tag"].GetString() != "None":
            self.__WriteShapeOptimizationConditions()
            self.__WriteDesignSurfaceAsSubModelPart()
        self.__WriteSubModelPartsWithoutDesignSurface()

        self.new_file.close()

        if self.interface_parameters["echo_level"].GetInt()==1:
            print("> Finished writing mdpa file!\n")

    # --------------------------------------------------------------------------
    def WriteNodesAsSU2MeshMotionFile(self,kratos_node_set,target_directory = "."):
        self.new_file = open(os.path.join(target_directory,self.interface_parameters["su2_related"]["mesh_motion_file"].GetString()), 'w')

        if self.su2_mesh_data["NDIME"] == 2:
            for kratos_node in kratos_node_set:
                su2_point_id = self.node_id_kratos_to_su2[kratos_node.Id]
                self.new_file.write(str(su2_point_id))
                self.new_file.write("\t")
                self.new_file.write(str("%.12f"%(kratos_node.X0)))
                self.new_file.write("\t")
                self.new_file.write(str("%.12f"%(kratos_node.Y0)))
                self.new_file.write("\n")
        else:
            for kratos_node in kratos_node_set:
                su2_point_id = self.node_id_kratos_to_su2[kratos_node.Id]
                self.new_file.write(str(su2_point_id))
                self.new_file.write("\t")
                self.new_file.write(str("%.12f"%(kratos_node.X0)))
                self.new_file.write("\t")
                self.new_file.write(str("%.12f"%(kratos_node.Y0)))
                self.new_file.write("\t")
                self.new_file.write(str("%.12f"%(kratos_node.Z0)))
                self.new_file.write("\n")

        self.new_file.close()

    # --------------------------------------------------------------------------
    def InitializeNewSU2Project(self):
        if self.interface_parameters["echo_level"].GetInt()==1:
            print("\n> Initializing SU2 project...")

        if os.path.isdir("DESIGNS"):
            shutil.rmtree("DESIGNS")

        su2_config = SU2.io.Config(self.interface_parameters["su2_related"]["config_file"].GetString())
        su2_config.NUMBER_PART = self.interface_parameters["su2_related"]["number_of_cores"].GetInt()
        su2_config.NZONES = self.interface_parameters["su2_related"]["number_of_zones"].GetInt()
        su2_config.GRADIENT_METHOD= self.interface_parameters["su2_related"]["gradient_method"].GetString()

        with suppress_stdout():
            state = SU2.io.State()
            state.find_files(su2_config)

            self.project = SU2.opt.Project(su2_config,state)
            self.project.config["MOTION_FILENAME"] = self.interface_parameters["su2_related"]["mesh_motion_file"].GetString()
            self.project.state.FILES["MOTION_FILE"] = self.interface_parameters["su2_related"]["mesh_motion_file"].GetString()

        if self.interface_parameters["echo_level"].GetInt()==1:
            print("> Finished initializing SU2 project!\n")

    # --------------------------------------------------------------------------
    def ComputeValues(self, list_of_response_ids, update_mesh, design_number):
        if self.interface_parameters["echo_level"].GetInt()==2:
            return self.project.f (list_of_response_ids, update_mesh, design_number)
        else:
            with suppress_stdout():
                return self.project.f (list_of_response_ids, update_mesh, design_number)

    # --------------------------------------------------------------------------
    def ComputeGradient(self, response_id, update_mesh, design_number):
        direct_solution_frequency = self.project.config["WRT_SOL_FREQ"]
        self.project.config["WRT_SOL_FREQ"] = self.interface_parameters["su2_related"]["write_frequency_adjoint_run"].GetInt()

        if self.interface_parameters["echo_level"].GetInt()==2:
            su2_gradient = self.project.df(response_id, update_mesh, design_number)[0]
        else:
            with suppress_stdout():
                su2_gradient = self.project.df(response_id, update_mesh, design_number)[0]

        kratos_gradient = {}
        for su2_node_id in su2_gradient.keys():
            kratos_gradient[self.node_id_su2_to_kratos[su2_node_id]] = su2_gradient[su2_node_id]

        self.project.config["WRT_SOL_FREQ"] = direct_solution_frequency

        return kratos_gradient

    # --------------------------------------------------------------------------
    def __ReadSU2Mesh(self):

        if self.interface_parameters["echo_level"].GetInt()==1:
            print("\n> Start reading SU2 mesh file...")

        ''' imports mesh and builds python dictionary structure
            su2_mesh_data                               mesh data dictionary
            su2_mesh_data["NDIME"]                      number of dimensions
            su2_mesh_data["NELEM"]                      number of elements
            su2_mesh_data["ELEM"]                       element array [ type, nodes, su_id ]
            su2_mesh_data["NPOIN"]                      number of points
            su2_mesh_data["POIN"]                       dict with point infos [su2_id: [x,y,z], su2_id: [x,y,z]...]
            su2_mesh_data["NMARK"]                      number of markers
            su2_mesh_data["MARKS"]                      marker data dictionary
            su2_mesh_data["MARKS"]["tag_name"]          marker data for 'tag_name'
            su2_mesh_data["MARKS"]["tag_name"]["NELEM"] number of elements
            su2_mesh_data["MARKS"]["tag_name"]["ELEM"]  element array [type,nodes]
            su2_mesh_data["MARKS"]["tag_name"]["POIN_IDS"]  stores su2_ids of points on marker
        '''

        # initialize variables
        self.su2_mesh_data["MARKS"]  = {}

        # open meshfile
        f_tb_read = open(self.interface_parameters["su2_related"]["mesh_file"].GetString(),'r')

        # readline helper functin
        def mesh_readlines(n_lines=1):
            fileslice = islice(f_tb_read,n_lines)
            return list(fileslice)

        # scan file until end of file
        keepon = True
        while keepon:

            # read line
            line = mesh_readlines()

            # stop if line is empty
            if not line:
                keepon = False
                break

            # fix white space
            line = line[0]
            line = line.replace('\t',' ')
            line = line.replace('\n',' ')

            # skip comments
            if line[0] == "%":
                pass

            # number of dimensions
            elif "NDIME=" in line:
                # save to SU2_MESH data
                self.su2_mesh_data["NDIME"] = int( line.split("=")[1].strip() )

            # elements
            elif "NELEM=" in line:

                # number of elements
                nelem = int( line.split("=")[1].strip() )
                # save to SU2_MESH data
                self.su2_mesh_data["NELEM"] = nelem

                # only read nelem lines
                fileslice = islice(f_tb_read,nelem)

                # the data pattern
                pattern = tuple( [int] + [int]*9 )

                # scan next lines for element data
                elem = [
                    [ t(s) for t,s in zip(pattern,line.split()) ]
                    for line in fileslice
                ]

                # save to SU2_MESH data
                self.su2_mesh_data["ELEM"] = elem

            # points
            elif "NPOIN=" in line:

                # number of points
                npoin = int( line.split("=")[1].strip().split(' ')[0] )
                # save to SU2_MESH data
                self.su2_mesh_data["NPOIN"] = npoin

                # only read npoin lines
                fileslice = islice(f_tb_read,npoin)

                # the data pattern
                pattern = None
                if self.su2_mesh_data["NDIME"] == 2:
                    pattern = tuple( [float]*2 + [int] )
                else:
                    pattern = tuple( [float]*3 + [int] )

                # scan next lines for element data
                poin = [
                    [ t(s) for t,s in zip(pattern,line.split()) ]
                    for line in fileslice
                ]

                # save to SU2_MESH data
                su2_point_ids = [row[-1] for row in poin]
                su2_point_coordinates = [row[0:-1] for row in poin]
                self.su2_mesh_data["POIN"] = dict(zip(su2_point_ids,su2_point_coordinates))

            # number of markers
            elif "NMARK=" in line:
                nmark = int( line.split("=")[1].strip() )
                # save to SU2_MESH data
                self.su2_mesh_data["NMARK"] = nmark

            # a marker
            elif "MARKER_TAG=" in line:
                # marker tag
                thistag = line.split("=")[1].strip()

                # start SU2_MARK dictionary
                thismark = {}

                # read number of marker elements
                line = mesh_readlines()[0]
                if not "MARKER_ELEMS=" in line:
                    raise Exception("Marker Specification Error")

                # convert string to int int
                thisnelem = int( line.split("=")[1].strip() )

                # save to SU2_MARK data
                thismark["NELEM"] = thisnelem

                # only read thisnelem lines
                fileslice = islice(f_tb_read,thisnelem)

                # the data pattern
                pattern = tuple( [int] + [int]*9 )

                # scan next lines for element data
                markelem = [
                    [ t(s) for t,s in zip(pattern,line.split()) ]
                    for line in fileslice
                ]

                # save to SU2_MARK data
                thismark["ELEM"] = markelem

                # Assign to marker their corresponding SU2 node IDs
                thismark["POIN_IDS"] = []
                for elem in thismark["ELEM"]:
                    # loop over all points
                    for itr in range(1,len(elem)):
                        su2_point_id = int("%i " % elem[itr])
                        if su2_point_id not in thismark["POIN_IDS"]:
                            thismark["POIN_IDS"].append(su2_point_id)

                # add to marker list
                self.su2_mesh_data["MARKS"][thistag] = thismark

        # Close read file
        f_tb_read.close()

        if self.interface_parameters["echo_level"].GetInt()==1:
            print("> Finished reading SU2 mesh file!\n")

    # --------------------------------------------------------------------------
    def __TranslateSU2IdsToKratosIds(self):
        # Translate node ids
        for su2_id in self.su2_mesh_data["POIN"]:
            if su2_id == 0:
                kratos_id = self.su2_mesh_data["NPOIN"]
                self.node_id_su2_to_kratos[su2_id] = kratos_id
                self.node_id_kratos_to_su2[kratos_id] = su2_id
            else:
                kratos_id = su2_id
                self.node_id_su2_to_kratos[su2_id] = kratos_id
                self.node_id_kratos_to_su2[kratos_id] = su2_id

        # Translate element ids
        for elem in self.su2_mesh_data["ELEM"]:
            su2_id = elem[-1]

            if su2_id == 0:
                kratos_id = self.su2_mesh_data["NELEM"]
                self.element_id_su2_to_kratos[su2_id] = kratos_id
            else:
                kratos_id = elem[-1]
                self.element_id_su2_to_kratos[su2_id] = kratos_id

    # --------------------------------------------------------------------------
    def __WriteHeader(self):
        self.new_file.write("Begin Properties 1\n")
        self.new_file.write("End Properties\n\n")

    # --------------------------------------------------------------------------
    def __WriteNodes(self):
        self.new_file.write("Begin Nodes\n")

        if self.su2_mesh_data["NDIME"] == 2:
            for su2_id, coordinates in sorted(self.su2_mesh_data["POIN"].items()):
                self.new_file.write("\t")
                self.new_file.write(str(self.node_id_su2_to_kratos[su2_id]))
                self.new_file.write("\t")
                self.new_file.write(str("%.12f"%(coordinates[0])))
                self.new_file.write("\t")
                self.new_file.write(str("%.12f"%(coordinates[1])))
                self.new_file.write("\t")
                self.new_file.write(str(0.0))
                self.new_file.write("\n")
        else:
            for su2_id, coordinates in sorted(self.su2_mesh_data["POIN"].items()):
                self.new_file.write("\t")
                self.new_file.write(str(self.node_id_su2_to_kratos[su2_id]))
                self.new_file.write("\t")
                self.new_file.write(str("%.12f"%(coordinates[0])))
                self.new_file.write("\t")
                self.new_file.write(str("%.12f"%(coordinates[1])))
                self.new_file.write("\t")
                self.new_file.write(str("%.12f"%(coordinates[2])))
                self.new_file.write("\n")

        self.new_file.write("End Nodes\n\n")

    # --------------------------------------------------------------------------
    def __WriteElements(self):

        # Identify all indices of the specific element types in list of elements
        line_entries = []
        tria_entries = []
        quad_entries = []
        tet_entries = []
        hex_entries = []
        prism_entries = []
        pyramid_entries = []
        for index,entry in enumerate(self.su2_mesh_data["ELEM"]):
            if entry[0] == self.su2_line_id:
                line_entries.append(index)
            if entry[0] == self.su2_tria_id:
                tria_entries.append(index)
            if entry[0] == self.su2_quad_id:
                quad_entries.append(index)
            if entry[0] == self.su2_tet_id:
                tet_entries.append(index)
            if entry[0] == self.su2_hex_id:
                hex_entries.append(index)
            if entry[0] == self.su2_prism_id:
                prism_entries.append(index)
            if entry[0] == self.su2_pyramid_id:
                pyramid_entries.append(index)

        # write triangular elements
        if tria_entries:
            self.new_file.write("Begin Elements GenericElement2D3N\n")
            for list_index in tria_entries:
                entry = self.su2_mesh_data["ELEM"][list_index]
                self.new_file.write("\t")
                self.new_file.write(str(self.element_id_su2_to_kratos[entry[-1]]))
                self.new_file.write("\t")
                self.new_file.write("1")
                self.new_file.write("\t")
                self.new_file.write(str(self.node_id_su2_to_kratos[entry[1]]))
                self.new_file.write("\t")
                self.new_file.write(str(self.node_id_su2_to_kratos[entry[2]]))
                self.new_file.write("\t")
                self.new_file.write(str(self.node_id_su2_to_kratos[entry[3]]))
                self.new_file.write("\n")
            self.new_file.write("End Elements\n\n")

        # write tet elements (not that pyramides are written as two tets)
        if tet_entries:
            self.new_file.write("Begin Elements GenericElement3D4N\n")
            for list_index in tet_entries:
                entry = self.su2_mesh_data["ELEM"][list_index]
                self.new_file.write("\t")
                self.new_file.write(str(self.element_id_su2_to_kratos[entry[-1]]))
                self.new_file.write("\t")
                self.new_file.write("1")
                self.new_file.write("\t")
                self.new_file.write(str(self.node_id_su2_to_kratos[entry[1]]))
                self.new_file.write("\t")
                self.new_file.write(str(self.node_id_su2_to_kratos[entry[2]]))
                self.new_file.write("\t")
                self.new_file.write(str(self.node_id_su2_to_kratos[entry[3]]))
                self.new_file.write("\t")
                self.new_file.write(str(self.node_id_su2_to_kratos[entry[4]]))
                self.new_file.write("\n")
            self.new_file.write("End Elements\n\n")

        # write pyramids as two tet elements
        if pyramid_entries:
            print("Waring!!!!! Pyramid elements included but not supported by Kratos. Hence, all pyramids are translated to two tets.")
            self.new_file.write("Begin Elements GenericElement3D4N\n")
            num_pyramids = 0
            for list_index in pyramid_entries:
                entry = self.su2_mesh_data["ELEM"][list_index]

                num_pyramids = num_pyramids + 1
                elemID_tet1 = self.element_id_su2_to_kratos[entry[-1]]
                elemID_tet2 = self.su2_mesh_data["NELEM"] + num_pyramids

                self.new_file.write("\t")
                self.new_file.write(str(elemID_tet1))
                self.new_file.write("\t")
                self.new_file.write("1")
                self.new_file.write("\t")
                self.new_file.write(str(self.node_id_su2_to_kratos[entry[2]]))
                self.new_file.write("\t")
                self.new_file.write(str(self.node_id_su2_to_kratos[entry[3]]))
                self.new_file.write("\t")
                self.new_file.write(str(self.node_id_su2_to_kratos[entry[1]]))
                self.new_file.write("\t")
                self.new_file.write(str(self.node_id_su2_to_kratos[entry[5]]))
                self.new_file.write("\n")

                self.new_file.write("\t")
                self.new_file.write(str(elemID_tet2))
                self.new_file.write("\t")
                self.new_file.write("1")
                self.new_file.write("\t")
                self.new_file.write(str(self.node_id_su2_to_kratos[entry[4]]))
                self.new_file.write("\t")
                self.new_file.write(str(self.node_id_su2_to_kratos[entry[1]]))
                self.new_file.write("\t")
                self.new_file.write(str(self.node_id_su2_to_kratos[entry[3]]))
                self.new_file.write("\t")
                self.new_file.write(str(self.node_id_su2_to_kratos[entry[5]]))
                self.new_file.write("\n")

            self.new_file.write("End Elements\n\n")

        # write prism elements
        if prism_entries:
            self.new_file.write("Begin Elements GenericElement3D6N\n")
            for list_index in prism_entries:
                entry = self.su2_mesh_data["ELEM"][list_index]
                self.new_file.write("\t")
                self.new_file.write(str(self.element_id_su2_to_kratos[entry[-1]]))
                self.new_file.write("\t")
                self.new_file.write("1")
                self.new_file.write("\t")
                self.new_file.write(str(self.node_id_su2_to_kratos[entry[4]]))
                self.new_file.write("\t")
                self.new_file.write(str(self.node_id_su2_to_kratos[entry[5]]))
                self.new_file.write("\t")
                self.new_file.write(str(self.node_id_su2_to_kratos[entry[6]]))
                self.new_file.write("\t")
                self.new_file.write(str(self.node_id_su2_to_kratos[entry[1]]))
                self.new_file.write("\t")
                self.new_file.write(str(self.node_id_su2_to_kratos[entry[2]]))
                self.new_file.write("\t")
                self.new_file.write(str(self.node_id_su2_to_kratos[entry[3]]))
                self.new_file.write("\n")
            self.new_file.write("End Elements\n\n")

        # write hex elements
        if hex_entries:
            self.new_file.write("Begin Elements GenericElement3D8N\n")
            for list_index in prism_entries:
                entry = self.su2_mesh_data["ELEM"][list_index]
                self.new_file.write("\t")
                self.new_file.write(str(self.element_id_su2_to_kratos[entry[-1]]))
                self.new_file.write("\t")
                self.new_file.write("1")
                self.new_file.write("\t")
                self.new_file.write(str(self.node_id_su2_to_kratos[entry[1]]))
                self.new_file.write("\t")
                self.new_file.write(str(self.node_id_su2_to_kratos[entry[2]]))
                self.new_file.write("\t")
                self.new_file.write(str(self.node_id_su2_to_kratos[entry[3]]))
                self.new_file.write("\t")
                self.new_file.write(str(self.node_id_su2_to_kratos[entry[4]]))
                self.new_file.write("\t")
                self.new_file.write(str(self.node_id_su2_to_kratos[entry[5]]))
                self.new_file.write("\t")
                self.new_file.write(str(self.node_id_su2_to_kratos[entry[6]]))
                self.new_file.write("\t")
                self.new_file.write(str(self.node_id_su2_to_kratos[entry[7]]))
                self.new_file.write("\t")
                self.new_file.write(str(self.node_id_su2_to_kratos[entry[8]]))
                self.new_file.write("\n")
            self.new_file.write("End Elements\n\n")

    # --------------------------------------------------------------------------
    def __WriteShapeOptimizationConditions(self):
        design_surface_tag = self.interface_parameters["su2_related"]["design_surface_tag"].GetString()

        # Identify all indices of the specific element types in list of elements on design surface
        line_entries = []
        tria_entries = []
        quad_entries = []
        for index,entry in enumerate(self.su2_mesh_data["MARKS"][design_surface_tag]["ELEM"]):
            if entry[0] == self.su2_line_id:
                line_entries.append(index)
            if entry[0] == self.su2_tria_id:
                tria_entries.append(index)
            if entry[0] == self.su2_quad_id:
                quad_entries.append(index)

        condition_id = 0

        # write line elements on design surface as line conditions
        if line_entries:
            self.new_file.write("Begin Conditions ShapeOptimizationCondition2D2N\n")
            for list_index in line_entries:
                entry = self.su2_mesh_data["MARKS"][design_surface_tag]["ELEM"][list_index]

                condition_id = condition_id+1
                self.list_of_optimization_conditions.append(condition_id)

                self.new_file.write("\t")
                self.new_file.write(str(condition_id))
                self.new_file.write("\t")
                self.new_file.write("1")
                self.new_file.write("\t")
                self.new_file.write(str(self.node_id_su2_to_kratos[entry[1]]))
                self.new_file.write("\t")
                self.new_file.write(str(self.node_id_su2_to_kratos[entry[2]]))
                self.new_file.write("\n")
            self.new_file.write("End Conditions\n\n")

        # write triangular elements on design surface as triangle conditions
        if tria_entries:
            self.new_file.write("Begin Conditions ShapeOptimizationCondition3D3N\n")
            for list_index in tria_entries:
                entry = self.su2_mesh_data["MARKS"][design_surface_tag]["ELEM"][list_index]

                condition_id = condition_id+1
                self.list_of_optimization_conditions.append(condition_id)

                self.new_file.write("\t")
                self.new_file.write(str(condition_id))
                self.new_file.write("\t")
                self.new_file.write("1")
                self.new_file.write("\t")
                self.new_file.write(str(self.node_id_su2_to_kratos[entry[1]]))
                self.new_file.write("\t")
                self.new_file.write(str(self.node_id_su2_to_kratos[entry[2]]))
                self.new_file.write("\t")
                self.new_file.write(str(self.node_id_su2_to_kratos[entry[3]]))
                self.new_file.write("\n")
            self.new_file.write("End Conditions\n\n")

        # write quad elements on design surface as quad conditions
        if quad_entries:
            self.new_file.write("Begin Conditions ShapeOptimizationCondition3D4N\n")
            for list_index in quad_entries:
                entry = self.su2_mesh_data["MARKS"][design_surface_tag]["ELEM"][list_index]

                condition_id = condition_id+1
                self.list_of_optimization_conditions.append(condition_id)

                self.new_file.write("\t")
                self.new_file.write(str(condition_id))
                self.new_file.write("\t")
                self.new_file.write("1")
                self.new_file.write("\t")
                self.new_file.write(str(self.node_id_su2_to_kratos[entry[1]]))
                self.new_file.write("\t")
                self.new_file.write(str(self.node_id_su2_to_kratos[entry[2]]))
                self.new_file.write("\t")
                self.new_file.write(str(self.node_id_su2_to_kratos[entry[3]]))
                self.new_file.write("\t")
                self.new_file.write(str(self.node_id_su2_to_kratos[entry[4]]))
                self.new_file.write("\n")
            self.new_file.write("End Conditions\n\n")

    # --------------------------------------------------------------------------
    def __WriteSubModelPartsWithoutDesignSurface(self):
        for mark_tag in self.su2_mesh_data["MARKS"]:

            if mark_tag == self.interface_parameters["su2_related"]["design_surface_tag"].GetString():
                continue

            self.new_file.write("Begin SubModelPart " + mark_tag + "\n")
            self.new_file.write("\tBegin SubModelPartNodes\n")
            for su2_point_id in self.su2_mesh_data["MARKS"][mark_tag]["POIN_IDS"]:
                self.new_file.write("\t\t")
                self.new_file.write(str(self.node_id_su2_to_kratos[su2_point_id]))
                self.new_file.write("\n")
            self.new_file.write("\tEnd SubModelPartNodes\n")
            self.new_file.write("End SubModelPart\n\n")

    # --------------------------------------------------------------------------
    def __WriteDesignSurfaceAsSubModelPart(self):
        design_surface_tag = self.interface_parameters["su2_related"]["design_surface_tag"].GetString()

        self.new_file.write("Begin SubModelPart " + design_surface_tag + "\n")
        self.new_file.write("\tBegin SubModelPartNodes\n")
        for su2_point_id in self.su2_mesh_data["MARKS"][design_surface_tag]["POIN_IDS"]:
            self.new_file.write("\t\t")
            self.new_file.write(str(self.node_id_su2_to_kratos[su2_point_id]))
            self.new_file.write("\n")
        self.new_file.write("\tEnd SubModelPartNodes\n")
        self.new_file.write("\tBegin SubModelPartConditions\n")
        for kratos_condition_id in self.list_of_optimization_conditions:
            self.new_file.write("\t\t")
            self.new_file.write(str(kratos_condition_id))
            self.new_file.write("\n")
        self.new_file.write("\tEnd SubModelPartConditions\n")
        self.new_file.write("End SubModelPart\n\n")

# ==============================================================================