from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics

# MPI
import KratosMultiphysics.mpi as mpi
import KratosMultiphysics.TrilinosApplication as TrilinosApplication
import KratosMultiphysics.MetisApplication as MetisApplication

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

def PartitionCreator(model_part, settings, verbosity, sync_conditions):
    # Creating the partitions
    domain_size = model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
            
    input_file_name = settings["model_import_settings"]["input_filename"].GetString()
    
    number_of_partitions = mpi.mpi.size  # we set it equal to the number of processors
    
    if mpi.mpi.size > 1:
        model_part_io = KratosMultiphysics.ModelPartIO(input_file_name)
        if mpi.mpi.rank == 0:

            partitioner = MetisApplication.MetisDivideHeterogeneousInputProcess( model_part_io, number_of_partitions, domain_size, verbosity, sync_conditions)
            partitioner.Execute()

        mpi.mpi.world.barrier()

        settings["model_import_settings"]["input_filename"].SetString(input_file_name + "_" + str(mpi.mpi.rank))

        MPICommSetup = MetisApplication.SetMPICommunicatorProcess(model_part)
        MPICommSetup.Execute()

        Comm = TrilinosApplication.CreateCommunicator()
    else:
        raise NameError('Your number of mpi.size is 1')
    
    print("MPI partitions created")
        
    return [model_part, settings, Comm, domain_size]