from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
try:
    import KratosMultiphysics.mpi as KratosMPI
    import KratosMultiphysics.MetisApplication as KratosMetis
    import KratosMultiphysics.TrilinosApplication as KratosTrilinos
    print("MPI imported, rank =", KratosMPI.mpi.rank)
    my_rank = KratosMPI.mpi.rank
    number_of_processors = KratosMPI.mpi.size
except:
    print("MPI not imported!")
    my_rank = 0   
    number_of_processors = 1
import KratosMultiphysics.MappingApplication as KratosMapping

def create_model_part(name, geom_range, number_nodes_per_rank, number_of_processors):
    model_part = KratosMultiphysics.ModelPart(name)

    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)

    model_part.AddProperties(KratosMultiphysics.Properties(1))
    props = model_part.GetProperties()[1]

    # number_nodes_per_rank = 5
    num_nodes_global = number_nodes_per_rank * number_of_processors
    # geom_range = [0,10]
    d_range = (geom_range[1] - geom_range[0]) / number_of_processors
    dx = d_range / (number_nodes_per_rank-1)
    local_min = geom_range[0] + d_range * my_rank
    create_line_conditions = True

    node_start_id = my_rank * number_nodes_per_rank
    for i in range(number_nodes_per_rank):
        node_id = node_start_id + i + 1
        x_coord = local_min + dx*i
        # print(my_rank, node_id, x_coord)
        model_part.CreateNewNode(node_id, x_coord, 0.0, 0.0)
            
    p_i = 1
    for node in model_part.Nodes:
        node.SetSolutionStepValue(KratosMultiphysics.PARTITION_INDEX, my_rank)
        node.SetSolutionStepValue(KratosMultiphysics.PRESSURE, p_i)
        p_i += 1


    # TODO: Create Elements
    # model_part.CreateNewElement("Element2D4N", i+1,[int(elements[4*i+0]),int(elements[4*i+1]),int(elements[4*i+2]),int(elements[4*i+3])], prp)

    # Jordi I think this is not needed, since the "ParallelFillCommunicator" also sets an mpi-communicator
    # It is also done in the trilinos import modelpart utility
    # MPICommSetupOrigin = KratosMetis.SetMPICommunicatorProcess(model_part_origin)
    # MPICommSetupOrigin.Execute()
    if number_of_processors > 1:
        ParallelFillCommOrigin = KratosTrilinos.ParallelFillCommunicator(model_part.GetRootModelPart())
        ParallelFillCommOrigin.Execute()

    return model_part



def print_custom(string2print):
    rank = 0
    try:
        rank = KratosMPI.mpi.rank
    except:
        pass
    if rank == 0:
        print(string2print)

def barrier_custom():
    try:
        KratosMPI.mpi.world.barrier()
    except:
        pass


        # model_part = ModelPart(model_part_name)
        # for variable in variable_list:
        #     model_part.AddNodalSolutionStepVariable(variable)

        # if (number_of_partitions > 1):
        #     if (mpi.size > 1):
        #         if (mpi.rank == 0):
        #             model_part_io = ReorderConsecutiveModelPartIO(model_part_input_file)

        #             partitioner = MetisDivideHeterogeneousInputProcess(
        #                 model_part_io,
        #                 number_of_partitions,
        #                 size_domain,
        #                 0, # verbosity, set to 1 for more detailed output
        #                 True)

        #             partitioner.Execute()

        #         mpi.world.barrier()
        #         model_part_input_file = model_part_input_file + "_" + str(mpi.rank)

        # model_part_io = ModelPartIO(model_part_input_file)
        # model_part_io.ReadModelPart(model_part)

        # if (number_of_partitions > 1):
        #     MPICommSetup = SetMPICommunicatorProcess(model_part)
        #     MPICommSetup.Execute()

        #     ParallelFillComm = ParallelFillCommunicator(model_part.GetRootModelPart())
        #     ParallelFillComm.Execute()

        # model_part.ProcessInfo.SetValue(DOMAIN_SIZE, size_domain)
        # model_part.SetBufferSize(1)

        # return model_part

geom_range_origin = [0,10]
geom_range_destination = [-10,10]
number_nodes_per_rank_origin = 3
number_nodes_per_rank_destination = 5
print_custom("\n##### Creating the ModelParts #####\n\n")
model_part_origin = create_model_part("origin", geom_range_origin, number_nodes_per_rank_origin, number_of_processors)
model_part_destination = create_model_part("destination", geom_range_destination, number_nodes_per_rank_destination, number_of_processors)

# model_part_origin.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
# model_part_destination.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
# model_part_origin.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)
# model_part_destination.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)

# model_part_origin.AddProperties(KratosMultiphysics.Properties(1))
# model_part_destination.AddProperties(KratosMultiphysics.Properties(1))
# props_origin = model_part_origin.GetProperties()[1]
# props_destination = model_part_destination.GetProperties()[1]




# num_nodes_global = number_nodes_per_rank * number_of_processors
# d_range = (geom_range[1] - geom_range[0]) / number_of_processors
# dx = d_range / (number_nodes_per_rank-1)
# local_min = d_range * my_rank
# create_line_conditions = True




# node_start_id = my_rank * number_nodes_per_rank
# for i in range(number_nodes_per_rank):
#     node_id = node_start_id + i + 1
#     x_coord = local_min + dx*i
#     print(my_rank, node_id, x_coord)
#     model_part_origin.CreateNewNode(node_id, x_coord, 0.0, 0.0)
#     model_part_destination.CreateNewNode(node_id, x_coord, 0.0, 0.0)
        
# for node in model_part_origin.Nodes:
#     node.SetSolutionStepValue(KratosMultiphysics.PARTITION_INDEX, my_rank)
# for node in model_part_destination.Nodes:
#     node.SetSolutionStepValue(KratosMultiphysics.PARTITION_INDEX, my_rank)


# # TODO: Create Elements
# # model_part.CreateNewElement("Element2D4N", i+1,[int(elements[4*i+0]),int(elements[4*i+1]),int(elements[4*i+2]),int(elements[4*i+3])], prp)

# # Jordi I think this is not needed, since the "ParallelFillCommunicator" also sets an mpi-communicator
# # It is also done in the trilinos import modelpart utility
# # MPICommSetupOrigin = KratosMetis.SetMPICommunicatorProcess(model_part_origin)
# # MPICommSetupOrigin.Execute()
# ParallelFillCommOrigin = KratosTrilinos.ParallelFillCommunicator(model_part_origin.GetRootModelPart())
# ParallelFillCommOrigin.Execute()

# # MPICommSetupDestination = KratosMetis.SetMPICommunicatorProcess(model_part_destination)
# # MPICommSetupDestination.Execute()
# ParallelFillCommDestination = KratosTrilinos.ParallelFillCommunicator(model_part_destination.GetRootModelPart())
# ParallelFillCommDestination.Execute()



mortar_mapper_settings = KratosMultiphysics.Parameters("""
        {
            "mapper_type": "Mortar",
            "interface_submodel_part_origin": "interface_chimera_background",
            "interface_submodel_part_destination": "Inlet2D_inlet"
        }
""")

nearest_element_mapper_settings = KratosMultiphysics.Parameters("""
        {
            "mapper_type": "NearestElement",
            "interface_submodel_part_origin": "interface_chimera_background",
            "interface_submodel_part_destination": "Inlet2D_inlet"
        }
""")

nearest_neighbor_matrix_mapper_settings = KratosMultiphysics.Parameters("""
        {
            "mapper_type": "NearestNeighborMatrixBased"
        }
""")
print_custom("\n##### Creating the mappers with the MapperFactory #####\n\n")

# mortar_mapper = KratosMapping.MapperFactoryNew.CreateMapper(model_part_origin, model_part_destination, mortar_mapper_settings)
# mortar_mapper.Map(KratosMultiphysics.PRESSURE, KratosMultiphysics.PRESSURE, KratosMapping.Mapper.ADD_VALUES)
# mortar_mapper.UpdateInterface()
# print_custom("\n\n")
# barrier_custom()

# nearest_element_mapper = KratosMapping.MapperFactoryNew.CreateMapper(model_part_origin, model_part_destination, nearest_element_mapper_settings)
# nearest_element_mapper.Map(KratosMultiphysics.PRESSURE, KratosMultiphysics.PRESSURE, KratosMapping.Mapper.ADD_VALUES)
# nearest_element_mapper.UpdateInterface()
# print_custom("\n\n")
# barrier_custom()

nearest_neighbor_matrix_mapper = KratosMapping.MapperFactoryNew.CreateMapper(model_part_origin, model_part_destination, nearest_neighbor_matrix_mapper_settings)
nearest_neighbor_matrix_mapper.Map(KratosMultiphysics.PRESSURE, KratosMultiphysics.PRESSURE, KratosMapping.Mapper.CONSERVATIVE)
# nearest_neighbor_matrix_mapper.InverseMap(KratosMultiphysics.PRESSURE, KratosMultiphysics.PRESSURE)
# nearest_neighbor_matrix_mapper.UpdateInterface()
# print_custom("\n\n")
# barrier_custom()

'''
print_custom("##### Creating the mappers DIRECTLY #####\n\n")
# 30.08.2017
# This does not work in a parallel compilation if the number of processors is 1 or executed in serial, the wrong space is used...
# Or can the Trilinos space also be used in serial / 1 process ?
mortar_mapper_2 = KratosMapping.MortarMapper(model_part_origin, model_part_destination, mortar_mapper_settings)
mortar_mapper_2.Map(KratosMultiphysics.PRESSURE, KratosMultiphysics.PRESSURE, KratosMapping.Mapper.ADD_VALUES)
mortar_mapper_2.UpdateInterface()
print_custom("\n\n")
barrier_custom()

# nearest_element_mapper_2 = KratosMapping.NearestElementMapper(model_part_origin, model_part_destination, nearest_element_mapper_settings)
# nearest_element_mapper_2.Map(KratosMultiphysics.PRESSURE, KratosMultiphysics.PRESSURE, KratosMapping.Mapper.ADD_VALUES)
# nearest_element_mapper_2.UpdateInterface()
# print_custom("\n\n")
# barrier_custom()

nearest_neighbor_matrix_mapper_2 = KratosMapping.NearestNeighborMapperMatrix(model_part_origin, model_part_destination, nearest_neighbor_matrix_mapper_settings)
nearest_neighbor_matrix_mapper_2.Map(KratosMultiphysics.PRESSURE, KratosMultiphysics.PRESSURE, KratosMapping.Mapper.ADD_VALUES)
nearest_neighbor_matrix_mapper_2.UpdateInterface()
print_custom("\n\n")
barrier_custom()
'''