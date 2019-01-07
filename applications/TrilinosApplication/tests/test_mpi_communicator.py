from __future__ import print_function, absolute_import, division

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics
import KratosMultiphysics.mpi as KratosMPI
import KratosMultiphysics.MetisApplication as KratosMetis
import KratosMultiphysics.TrilinosApplication as KratosTrilinos


def GetFilePath(fileName):
    return os.path.dirname(os.path.realpath(__file__)) + "/" + fileName


class TestMPICommunicator(KratosUnittest.TestCase):

    def _read_model_part_mpi(self,main_model_part):
        
        if(KratosMPI.mpi.size == 1):
            self.skipTest("Test can be run only using more than one mpi process")
        
        ## Add variables to the model part
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)  
     
        ## Serial partition of the original .mdpa file
        input_filename = "test_mpi_communicator"   
        if KratosMPI.mpi.rank == 0 :

            # Original .mdpa file reading
            model_part_io = KratosMultiphysics.ModelPartIO(input_filename)

            # Partition of the original .mdpa file
            number_of_partitions = KratosMPI.mpi.size # Number of partitions equals the number of processors
            domain_size = main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
            verbosity = 0
            sync_conditions = True # Make sure that the condition goes to the same partition as the element is a face of
            
            partitioner = KratosMetis.MetisDivideHeterogeneousInputProcess(model_part_io, number_of_partitions , domain_size, verbosity, sync_conditions)
            partitioner.Execute() 

            print("Metis divide finished.")

        KratosMPI.mpi.world.barrier()
    
        ## Read the partitioned .mdpa files
        mpi_input_filename = input_filename + "_" + str(KratosMPI.mpi.rank)
        model_part_io = KratosMultiphysics.ModelPartIO(mpi_input_filename)
        model_part_io.ReadModelPart(main_model_part)
    
        ## Construct and execute the MPICommunicator
        KratosMetis.SetMPICommunicatorProcess(main_model_part).Execute()
        
        ## Construct and execute the Parallel fill communicator
        ParallelFillCommunicator = KratosTrilinos.ParallelFillCommunicator(main_model_part.GetRootModelPart())
        ParallelFillCommunicator.Execute()
      
        ## Check submodelpart of each main_model_part of each processor
        self.assertTrue(main_model_part.HasSubModelPart("Skin"))
        skin_sub_model_part = main_model_part.GetSubModelPart("Skin")
        
        
    def test_assemble_variable_in_model_part(self):
        current_model = KratosMultiphysics.Model()
        main_model_part = current_model.CreateModelPart("MainModelPart")
        main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)
                
        self._read_model_part_mpi(main_model_part)
        
        ## Initialize DENSITY variable
        for node in main_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.DENSITY, 0, 0.0)
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 0, 0.0)
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 0, 0.0)
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z, 0, 0.0)
            
        ## Fill the DENSITY value in each processor
        for elem in main_model_part.Elements:
            for node in elem.GetNodes():
                d = node.GetSolutionStepValue(KratosMultiphysics.DENSITY)
                d += 1
                node.SetSolutionStepValue(KratosMultiphysics.DENSITY,0,d)
                node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X,0,d)
                node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y,0,d)
                node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z,0,d)
        
        ## Communicate between processors
        main_model_part.GetCommunicator().AssembleCurrentData(KratosMultiphysics.DENSITY)
        main_model_part.GetCommunicator().AssembleCurrentData(KratosMultiphysics.DISPLACEMENT)

        ## Expected nodal values after communication
        expected_values = {1 : 2,
                           2 : 3,
                           3 : 1,
                           4 : 3,
                           5 : 6,
                           6 : 3,
                           7 : 1,
                           8 : 3,
                           9 : 2}
        
        ## Check the obtained values
        for node in main_model_part.Nodes:
            self.assertEqual(node.GetSolutionStepValue(KratosMultiphysics.DENSITY), expected_values[node.Id])
            self.assertEqual(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X), expected_values[node.Id])
            self.assertEqual(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y), expected_values[node.Id])
            self.assertEqual(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z), expected_values[node.Id])
            
                
    def test_assemble_variable_in_sub_model_part(self):
        current_model = KratosMultiphysics.Model()
        main_model_part = current_model.CreateModelPart("MainModelPart")
        main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)
         
        self._read_model_part_mpi(main_model_part)
       
        submodelpart = main_model_part.GetSubModelPart("Skin")      
        
        KratosMPI.mpi.world.barrier()
        
        # Check partitioning
        #~ for i in range(KratosMPI.mpi.size):
            #~ if(KratosMPI.mpi.rank == i):
                #~ print(" ")
                #~ print("Submodelpart elements of processor id = ", i)
                #~ for elem in submodelpart.Elements:
                    #~ print(elem.Id)
                #~ print("Submodelpart conditions of processor id = ", i)
                #~ for elem in submodelpart.Conditions:
                    #~ print(elem.Id)
                #~ print("Submodelpart local nodes of processor id = ", i)
                #~ for node in submodelpart.GetCommunicator().LocalMesh().Nodes:
                    #~ print(node.Id)
                #~ print("Submodelpart ghost nodes of processor id = ", i)
                #~ for node in submodelpart.GetCommunicator().GhostMesh().Nodes:
                    #~ print(node.Id)
                #~ print("Submodelpart interface nodes of processor id = ", i)
                #~ for node in submodelpart.GetCommunicator().InterfaceMesh().Nodes:
                    #~ print(node.Id)
                #~ print(" ")
            #~ KratosMPI.mpi.world.barrier()
        
        for node in submodelpart.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.DENSITY, 0, 0.0)
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 0, 0.0)
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 0, 0.0)
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z, 0, 0.0)
            
        for condition in submodelpart.Conditions:
            for node in condition.GetNodes():
                d = node.GetSolutionStepValue(KratosMultiphysics.DENSITY)
                d += 1
                node.SetSolutionStepValue(KratosMultiphysics.DENSITY,0,d)
                node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 0, d)
                node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 0, d)
                node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z, 0, d)
        
        submodelpart.GetCommunicator().AssembleCurrentData(KratosMultiphysics.DENSITY)
        submodelpart.GetCommunicator().AssembleCurrentData(KratosMultiphysics.DISPLACEMENT)
                
        for node in submodelpart.Nodes:
            self.assertEqual(node.GetSolutionStepValue(KratosMultiphysics.DENSITY), 2.0)
            self.assertEqual(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X), 2.0)
            self.assertEqual(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y), 2.0)
            self.assertEqual(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z), 2.0)



    #def test_model_part_io_properties_block(self):
    #    model_part = ModelPart("Main")
    #    model_part_io = ModelPartIO("test_model_part_io")
    #    model_part_io.ReadProperties(model_part.Properties)

if __name__ == '__main__':
    KratosUnittest.main()
