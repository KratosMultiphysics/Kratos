from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
if not "DO_NOT_PARTITION_DOMAIN" in os.environ:
    from KratosMultiphysics.MetisApplication import *
    from KratosMultiphysics.MPISearchApplication import *
    import DEM_material_test_script_mpi as DEM_material_test_script
else:
    import DEM_material_test_script

from KratosMultiphysics.mpi import *
import DEM_procedures

from glob import glob

class MdpaCreator(DEM_procedures.MdpaCreator):

    def __init__(self, path, DEM_parameters):
        super(MdpaCreator,self).__init__(path, DEM_parameters)


class GranulometryUtils(DEM_procedures.GranulometryUtils):

    def __init__(self, domain_volume, model_part):
        super(GranulometryUtils,self).__init__(domain_volume, model_part)


class PostUtils(DEM_procedures.PostUtils):

    def __init__(self, DEM_parameters, balls_model_part):
        super(PostUtils,self).__init__(DEM_parameters, balls_model_part)


class Procedures(DEM_procedures.Procedures):

    def __init__(self, DEM_parameters):
        super(Procedures,self).__init__(DEM_parameters)

    def Barrier(self):
        mpi.world.barrier()

    def AddMpiVariables(self, model_part):
        model_part.AddNodalSolutionStepVariable(PARTITION_INDEX)
        model_part.AddNodalSolutionStepVariable(PARTITION_MASK)

    def CreateDirectories(self, main_path, problem_name):

        root             = main_path + '/' + problem_name
        post_path        = root + '_Post_Files'
        data_and_results = root + '_Results_and_Data'
        graphs_path      = root + '_Graphs'
        MPI_results      = root + '_MPI_results'

        if mpi.rank == 0:
            for directory in [post_path, data_and_results, graphs_path, MPI_results]:
                if not os.path.isdir(directory):
                    os.makedirs(str(directory))

        self.Barrier()

        return [post_path, data_and_results, graphs_path, MPI_results]

    def PreProcessModel(self, DEM_parameters):
        if (mpi.rank == 0):
            print("Creating MPIer...")
            #MPIClassObject = MPIer.MPIerClass(str(DEM_parameters["problem_name"].GetString()) + "DEM.mdpa")
            print("done.")
        self.Barrier() #TODO: maybe not necessary (debugging)

    def FindMaxNodeIdInModelPart(self, model_part):

        node_max = super(Procedures,self).FindMaxNodeIdInModelPart(model_part)
        node_max_gath = mpi.allgather(mpi.world,node_max)
        total_max = reduce(lambda x,y: max(x,y), node_max_gath)
        return total_max

    def DeleteFiles(self):
        if (mpi.rank == 0):
            files_to_delete_list = glob('*.time')
            for to_erase_file in files_to_delete_list:
                os.remove(to_erase_file)

    def KRATOSprint(self, message):
        if (mpi.rank == 0):
            Logger.Print(*args, label="DEM")
            Logger.Flush()


class DEMFEMProcedures(DEM_procedures.DEMFEMProcedures):
    def PrintGraph(self, time):
        if (mpi.rank == 0):
            super(DEMFEMProcedures,self).PrintGraph(time)

    def FinalizeGraphs(self,RigidFace_model_part):
        if (mpi.rank == 0):
            super(DEMFEMProcedures,self).FinalizeGraphs(RigidFace_model_part)

    def FinalizeBallsGraphs(self,spheres_model_part):
        if (mpi.rank == 0):
            super(DEMFEMProcedures,self).FinalizeBallsGraphs(spheres_model_part)


class Report(DEM_procedures.Report):

    def __init__(self):
        super(Report,self).__init__()


class MaterialTest(DEM_procedures.MaterialTest):

    def __init__(self):
        super(MaterialTest,self).__init__()

    # Important: This has to be defined here as the imports from
    # the superclase and the derived clase are different
    def Initialize(self, DEM_parameters, procedures, solver, graphs_path, post_path, balls_model_part, rigid_face_model_part):
        self.TestType = DEM_parameters["TestType"].GetString()

        if (self.TestType != "None"):
            self.script = DEM_material_test_script.MaterialTest(DEM_parameters, procedures, solver, graphs_path, post_path, balls_model_part, rigid_face_model_part)
            self.script.Initialize()


class MultifileList(object):

    def __init__(self,name,step):
        self.index = 0
        self.step = step
        self.name = name
        self.which_folder = which_folder
        if which_folder == "inner":
            absolute_path_to_file = os.path.join(post_path, "_list_" + self.name + "_" + str(mpi.rank) + "_" + str(step) + ".post.lst")
        else:
            absolute_path_to_file = os.path.join(post_path, self.name + "_" + str(mpi.rank) + ".post.lst")
        #THIS BREAKS THE AUTOMATIC OPENING OF POSTPROCESS FILES IN GID WHEN SWITCHING TO POST

        self.file = open(absolute_path_to_file, "w")


class DEMIo(DEM_procedures.DEMIo):

    def __init__(self, DEM_parameters, post_path, all_model_parts):
        super(DEMIo,self).__init__(DEM_parameters, post_path, all_model_parts)

    def AddMpiVariables(self):
        self.spheres_variables.append(PARTITION_INDEX)

    def SetOutputName(self,name):
        self.gid_io.ChangeOutputName(name + "_" + str(mpi.rank))

    def GetMultiFileListName(self, name):
        return name + "_" + str(mpi.rank)


class ParallelUtils(DEM_procedures.ParallelUtils):

    def __init__(self):
        super(ParallelUtils,self).__init__()
        self.mpi_utilities = MpiUtilities()

    def Repart(self, balls_model_part):
        self.mpi_utilities.Repart(balls_model_part, 0, 1)

    def CalculateModelNewIds(self, balls_model_part):
        self.mpi_utilities.CalculateModelNewIds(balls_model_part, 40000)

    def PerformInitialPartition(self, model_part_io_spheres):
        domain_size = 3

        #print("(" + str(mpi.rank) + "," + str(mpi.size) + ")" + "before performing the division")
        number_of_partitions = mpi.size

        if mpi.rank == 0:
            #print("(" + str(mpi.rank) + "," + str(mpi.size) + ")" + "start partition process")
            partitioner = MetisDivideNodalInputToPartitionsProcess(model_part_io_spheres, number_of_partitions, domain_size);
            partitioner.Execute()

        #print("(" + str(mpi.rank) + "," + str(mpi.size) + ")" + "division performed")
        mpi.world.barrier()
        #return model_part_io_spheres

    def SetCommunicator(self, spheres_model_part, model_part_io_spheres, spheres_mp_filename):

        MPICommSetup = SetMPICommunicatorProcess(spheres_model_part)
        MPICommSetup.Execute()

        print("(" + str(mpi.rank) + "," + str(mpi.size) + ")" + "Communicator Set")
        print("(" + str(mpi.rank) + "," + str(mpi.size) + ")" + "Reading: "+spheres_mp_filename+"_"+str(mpi.rank))

        my_input_filename = spheres_mp_filename + "_" + str(mpi.rank)
        model_part_io_spheres = ModelPartIO(my_input_filename)

        return [model_part_io_spheres, spheres_model_part, MPICommSetup]

    def GetSearchStrategy(self, solver, model_part):
        return MPI_DEMSearch(model_part.GetCommunicator())

class SetOfModelParts(DEM_procedures.SetOfModelParts):
    pass

class DEMEnergyCalculator(DEM_procedures.DEMEnergyCalculator):
    pass
