import KratosMultiphysics
import KratosMultiphysics.kratos_utilities as kratos_utils

import os

metis_available = False
try: # Importing MPI before mapping makes the logo only print once
    import KratosMultiphysics.MetisApplication as KratosMetis
    metis_available = True
except ImportError:
    print("KratosMPI or MetisApplication is ot available, many MPI tests will be skipped!")
    warn_msg = "KratosMPI or MetisApplication is ot available, many MPI tests will be skipped!\n"
    Logger.PrintWarning("\nMappingApplication", warn_msg)

def CreatePartitionedMdpaFiles2D(number_of_partitions, mdpa_file_name):
    # Preparations
    folder_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "partitioned_mdpa_files")
    kratos_utils.DeleteDirectoryIfExisting(folder_path)

    if kratos_utils.IsRankZero():
        os.mkdir(folder_path)

        # Create a ModelPart
        node_1 = KratosMultiphysics.Node(1, 0.0, 0.0, 0.0)
        node_2 = KratosMultiphysics.Node(2, 0.0, 1.0, 0.0)
        node_3 = KratosMultiphysics.Node(3, 1.0, 1.0, 0.0)
        node_4 = KratosMultiphysics.Node(4, 1.0, 0.0, 0.0)

        geometry = KratosMultiphysics.Quadrilateral2D4(node_1, node_2, node_3, node_4)

        model_part = KratosMultiphysics.ModelPart("generated")

        settings = KratosMultiphysics.Parameters("""{
            "number_of_divisions" : 0,
            "element_name"        : "Element2D3N",
            "create_skin_sub_model_part": false
        }""")

        settings["number_of_divisions"].SetInt(number_of_partitions*3)

        mesh_generator = KratosMultiphysics.StructuredMeshGeneratorProcess(geometry,model_part, settings)
        mesh_generator.Execute()

        name_out_file = os.path.join(folder_path, mdpa_file_name)
        print(name_out_file)
        file = open(name_out_file + ".mdpa","w")
        file.close()
        KratosMultiphysics.ModelPartIO(name_out_file, KratosMultiphysics.IO.WRITE).WriteModelPart(model_part)

    if metis_available and number_of_partitions > 1 and kratos_utils.IsRankZero():
        print("Performing Partitioning")

        # Original .mdpa file reading
        model_part_io = KratosMultiphysics.ReorderConsecutiveModelPartIO(name_out_file)

        # Partition of the original .mdpa file
        domain_size = 2 # self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        verbosity = 0 #self.settings["echo_level"].GetInt()
        sync_conditions = True # Make sure that the condition goes to the same partition as the element is a face of
        partitioner = KratosMetis.MetisDivideHeterogeneousInputProcess(model_part_io,
                                                                    number_of_partitions,
                                                                    domain_size,
                                                                    verbosity,
                                                                    sync_conditions)
        partitioner.Execute()

        KratosMultiphysics.Logger.PrintInfo("::[TrilinosImportModelPartUtility]::", "Metis divide finished.")
    elif not metis_available and kratos_utils.IsRankZero():
        print("No partitioning performed")

if __name__ == '__main__':
    try: # Importing MPI before mapping makes the logo only print once
        import KratosMultiphysics.mpi as KratosMPI
        number_of_partitions = KratosMPI.mpi.size
    except:
        number_of_partitions = 1
    CreatePartitionedMdpaFiles(number_of_partitions)