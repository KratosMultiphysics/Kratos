from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.MetisApplication import *
from KratosMultiphysics.MPISearchApplication import *
from KratosMultiphysics.mpi import *

import DEM_procedures

import DEM_material_test_script_mpi as DEM_material_test_script

import MPIer


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

    def AddMpiVariables(self, model_part):
        model_part.AddNodalSolutionStepVariable(PARTITION_INDEX)
        model_part.AddNodalSolutionStepVariable(PARTITION_MASK)

    def CreateDirectories(self, main_path, problem_name):
        root             = main_path + '/' + problem_name

        post_path        = root + '_Post_Files'
        list_path        = root + '_Post_Lists'
        data_and_results = root + '_Results_and_Data'
        graphs_path      = root + '_Graphs'
        MPI_results      = root + '_MPI_results'

        if mpi.rank == 0:
            for directory in [post_path, list_path, data_and_results, graphs_path, MPI_results]:
                if not os.path.isdir(directory):
                    os.makedirs(str(directory))

        mpi.world.barrier()

        return [post_path,list_path,data_and_results,graphs_path,MPI_results]

    def PreProcessModel(self, DEM_parameters):
        if (mpi.rank == 0):
            MPIClassObject = MPIer.MPIerClass(str(DEM_parameters.problem_name) + "DEM.mdpa")
            
    def FindMaxNodeIdInModelPart(self, model_part):

        node_max = super(Procedures,self).FindMaxNodeIdInModelPart(model_part)
        node_max_gath = mpi.allgather(mpi.world,node_max)
        total_max = reduce(lambda x,y: max(x,y), node_max_gath)
        return total_max
        

    def KRATOSprint(self, message):
        if (mpi.rank == 0):
            print(message)
            sys.stdout.flush()


class DEMFEMProcedures(DEM_procedures.DEMFEMProcedures):
    def __init__(self, DEM_parameters, graphs_path, balls_model_part, RigidFace_model_part):
        super(DEMFEMProcedures,self).__init__(DEM_parameters,graphs_path,balls_model_part,RigidFace_model_part)

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
        self.TestType = DEM_parameters.TestType

        if (self.TestType != "None"):
            self.script = DEM_material_test_script.MaterialTest(DEM_parameters, procedures, solver, graphs_path, post_path, balls_model_part, rigid_face_model_part)
            self.script.Initialize()


class MultifileList(object):

    def __init__(self,name,step):
        self.index = 0
        self.step = step
        self.name = name
        self.file = open(self.name+"_"+str(mpi.rank)+"_"+str(step)+".post.lst","w")
        

class DEMIo(DEM_procedures.DEMIo):

    def __init__(self):
        super(DEMIo,self).__init__()

    def EnableMpiVariables(self):
        self.spheres_variables.append(PARTITION_INDEX)

    def SetOutputName(self,name):
        self.gid_io.ChangeOutputName(name + "_" + str(mpi.rank))

    def PrintMultifileLists(self,time, post_path):
        for mfilelist in self.multifilelists:
            
            if mfilelist.index == mfilelist.step:
                
                if (self.encoding == GiDPostMode.GiD_PostBinary):
                    mfilelist.file.write(post_path+"/"+mfilelist.name+"_"+str(mpi.rank)+"_"+"%.12g"%time+".post.bin\n")
                else:
                    mfilelist.file.write(post_path+"/"+mfilelist.name+"_"+str(mpi.rank)+"_"+"%.12g"%time+".post.msh\n")
                    mfilelist.file.write(post_path+"/"+mfilelist.name+"_"+str(mpi.rank)+"_"+"%.12g"%time+".post.res\n")
                mfilelist.file.flush()
                mfilelist.index = 0
            mfilelist.index += 1


class ParallelUtils(DEM_procedures.ParallelUtils):

    def __init__(self):
        super(ParallelUtils,self).__init__()
        self.mpi_utilities = MpiUtilities()

    def Repart(self, balls_model_part):
        self.mpi_utilities.Repart(balls_model_part, 0, 1)

    def CalculateModelNewIds(self, balls_model_part):
        self.mpi_utilities.CalculateModelNewIds(balls_model_part, 40000)

    def PerformInitialPartition(self, model_part, model_part_io_solid, input_file_name):
        domain_size = 3

        print("(" + str(mpi.rank) + "," + str(mpi.size) + ")" + "before performing the division")
        number_of_partitions = mpi.size
        
        if mpi.rank == 0:
            print("(" + str(mpi.rank) + "," + str(mpi.size) + ")" + "start partition process")
            partitioner = MetisDivideNodalInputToPartitionsProcess(model_part_io_solid, number_of_partitions, domain_size);
            partitioner.Execute()

        print("(" + str(mpi.rank) + "," + str(mpi.size) + ")" + "division performed")
        mpi.world.barrier()

        MPICommSetup = SetMPICommunicatorProcess(model_part)
        MPICommSetup.Execute()

        print("(" + str(mpi.rank) + "," + str(mpi.size) + ")" + "Communicator Set")
        print("(" + str(mpi.rank) + "," + str(mpi.size) + ")" + "Reading: "+input_file_name+"_"+str(mpi.rank))

        my_input_filename = input_file_name + "_" + str(mpi.rank)
        model_part_io_solid = ModelPartIO(my_input_filename, True)

        return [model_part_io_solid, model_part, MPICommSetup]

    def GetSearchStrategy(self, solver, model_part):
        return MPI_DEMSearch(model_part.GetCommunicator())
