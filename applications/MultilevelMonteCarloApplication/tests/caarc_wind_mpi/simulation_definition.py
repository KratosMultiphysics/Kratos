# Import Python libraries
import os
import numpy as np
import pickle

# Import Kratos
import KratosMultiphysics
import KratosMultiphysics.mpi
import KratosMultiphysics.MultilevelMonteCarloApplication
import KratosMultiphysics.MappingApplication
from applications.MultilevelMonteCarloApplication.tests.caarc_wind_mpi.FluidDynamicsAnalysisMC import FluidDynamicsAnalysisMC
from KratosMultiphysics.FluidDynamicsApplication import check_and_prepare_model_process_fluid

# Avoid printing of Kratos informations
# KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)


class SimulationScenario(FluidDynamicsAnalysisMC):
    def __init__(self,input_model,input_parameters,sample):
        super().__init__(input_model,input_parameters)
        self.sample = sample
        self.mapping = False
        self.interest_model_part = "FluidModelPart.NoSlip3D_structure"
        self.number_instances_time_power_sums = 0
        self.IsVelocityFieldPerturbed = False
        self.filename = "filename"

    def ModifyInitialProperties(self):
        """
        function changing print to file settings
        input:  self: an instance of the class
        """
        super().ModifyInitialProperties()
        random_inlet = max(0.1,self.sample[1])
        self.project_parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["value"][0].SetDouble(random_inlet)
        if (self.project_parameters["problem_data"]["perturbation"]["type"].GetString() == "correlated"):
            self.project_parameters["processes"]["initial_conditions_process_list"][0]["Parameters"]["seed"].SetInt(self.sample[0])

    def ComputeNeighbourElements(self):
        """
        function computing neighbour elements, required by our boundary conditions
        input:  self: an instance of the class
        """
        tmoc = KratosMultiphysics.TetrahedralMeshOrientationCheck
        throw_errors = False
        flags = (tmoc.COMPUTE_NODAL_NORMALS).AsFalse() | (tmoc.COMPUTE_CONDITION_NORMALS).AsFalse()
        flags |= tmoc.ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS
        KratosMultiphysics.TetrahedralMeshOrientationCheck(self._GetSolver().main_model_part.GetSubModelPart("fluid_computational_model_part"),throw_errors, flags).Execute()

    def Initialize(self):
        """
        function initializing moment estimator array
        input:  self: an instance of the class
        """
        super().Initialize()
        # compute neighbour elements required for current boundary conditions and not automatically run due to remeshing
        self.ComputeNeighbourElements()
        # initialize moment estimator array for each qoi to build time power sums
        self.moment_estimator_array = [[[0.0],[0.0]] for _ in range (0,1)] # +1 is for drag force x
        if (self.mapping is True):
            power_sums_parameters = KratosMultiphysics.Parameters("""{
                "reference_variable_name": "PRESSURE"
                }""")
            self.power_sums_process_mapping = KratosMultiphysics.MultilevelMonteCarloApplication.PowerSumsStatistics(self.mapping_reference_model.GetModelPart(self.interest_model_part),power_sums_parameters)
            self.power_sums_process_mapping.ExecuteInitialize()
        else:
            power_sums_parameters = KratosMultiphysics.Parameters("""{
                "reference_variable_name": "PRESSURE"
                }""")
            self.power_sums_process = KratosMultiphysics.MultilevelMonteCarloApplication.PowerSumsStatistics(self.model.GetModelPart(self.interest_model_part),power_sums_parameters)
            self.power_sums_process.ExecuteInitialize()

    def FinalizeSolutionStep(self):
        """
        function applying mapping if required and updating moment estimator array
        input:  self: an instance of the class
        """
        super().FinalizeSolutionStep()
        # run if current index is index of interest
        if (self.is_current_index_maximum_index is True):
            # avoid burn-in time
            if (self.model.GetModelPart(self.interest_model_part).ProcessInfo.GetPreviousTimeStepInfo().GetValue(KratosMultiphysics.TIME) >= \
                self.project_parameters["problem_data"]["burnin_time"].GetDouble()):
                # update number of contributions to time power sums
                self.number_instances_time_power_sums = self.number_instances_time_power_sums + 1
                # update power sums of drag force x
                self.moment_estimator_array[0][0][0] = self.moment_estimator_array[0][0][0] + self.current_drag_force_x
                self.moment_estimator_array[0][1][0] = self.moment_estimator_array[0][1][0] + self.current_drag_force_x**2
                if (self.mapping is True):
                    # call parallel fill communicator
                    ParallelFillCommunicator = KratosMultiphysics.mpi.ParallelFillCommunicator(self.mapping_reference_model.GetModelPart("FluidModelPart"))
                    ParallelFillCommunicator.Execute()
                    # mapping from current model part of interest to reference model part the pressure
                    mapping_parameters = KratosMultiphysics.Parameters("""{
                        "mapper_type": "nearest_element",
                        "interface_submodel_part_origin": "NoSlip3D_structure",
                        "interface_submodel_part_destination": "NoSlip3D_structure",
                        "echo_level" : 3
                        }""")
                    mapper = KratosMultiphysics.MappingApplication.MapperFactory.CreateMPIMapper(self._GetSolver().main_model_part,self.mapping_reference_model.GetModelPart("FluidModelPart"),mapping_parameters)
                    mapper.Map(KratosMultiphysics.PRESSURE,KratosMultiphysics.PRESSURE)
                    # pressure field power sums
                    self.power_sums_process_mapping.ExecuteFinalizeSolutionStep()
                else:
                    # update pressure field power sums
                    self.power_sums_process.ExecuteFinalizeSolutionStep()
        else:
            pass

    def EvaluateQuantityOfInterest(self):
        """
        function evaluating the QoI of the problem: lift coefficient
        input:  self: an instance of the class
        """

        communicator = KratosMultiphysics.DataCommunicator.GetDefault()
        rank = communicator.Rank()

        # run if current index is index of interest
        if (self.is_current_index_maximum_index is True):
            qoi_list = []

            # append time average drag force
            qoi_list.append(self.mean_drag_force_x)

            # append time average pressure
            if (self.mapping is not True):
                model_part_of_interest = self.model.GetModelPart(self.interest_model_part)
            elif (self.mapping is True):
                model_part_of_interest = self.mapping_reference_model.GetModelPart(self.interest_model_part)
            pressure_list = \
                [node.GetValue(KratosMultiphysics.CONTACT_PRESSURE) for node in model_part_of_interest.Nodes if node.GetSolutionStepValue(KratosMultiphysics.PARTITION_INDEX)==rank]
            pressure_list_new = []

            qoi_averaged_pressure = []
            if rank==0:
                pressure_list_new = communicator.GathervDoubles(pressure_list, 0)
            else:
                communicator.GathervDoubles(pressure_list, 0)
            if rank==0:
                for list_pressure in pressure_list_new:
                    qoi_averaged_pressure.extend(list_pressure)
            qoi_list.append(qoi_averaged_pressure)

            # append number of contributions to the power sums list
            self.moment_estimator_array[0].append(self.number_instances_time_power_sums) # drag force x
            # append drag force x time series power sums
            qoi_list.append(self.moment_estimator_array[0]) # drag force x

            # append pressure time series power sums
            pressure_power_sums_list = []
            for node in model_part_of_interest.Nodes:
                if node.GetSolutionStepValue(KratosMultiphysics.PARTITION_INDEX) == rank:
                    S1 = node.GetValue(KratosMultiphysics.MultilevelMonteCarloApplication.POWER_SUM_1)
                    S2 = node.GetValue(KratosMultiphysics.MultilevelMonteCarloApplication.POWER_SUM_2)
                    M = self.number_instances_time_power_sums
                    power_sums = [[S1],[S2],M]
                    pressure_power_sums_list.append(power_sums)

            string_to_send = pickle.dumps(pressure_power_sums_list, 2).decode("latin1")
            pressure_power_sums_list_new = []
            if rank==0:
                pressure_power_sums_list_new.extend(pressure_power_sums_list)
                for rank_recv in range(1, communicator.Size()):
                    pickled_list = communicator.SendRecvString(string_to_send, rank_recv, rank_recv)
                    list_to_extend = pickle.loads(pickled_list.encode("latin1"))
                    pressure_power_sums_list_new.extend(list_to_extend)
            else:
                pickled_list = communicator.SendRecvString(string_to_send, 0, 0)
            qoi_list.append(pressure_power_sums_list_new)

            # qoi_list is complete only in rank 0
            # we synchronize with other ranks
            string_to_send_qoi = pickle.dumps(qoi_list, 2).decode("latin1")
            if rank==0:
                for rank_recv in range(1, communicator.Size()):
                    pickled_list_qoi = communicator.SendRecvString(string_to_send_qoi, rank_recv, rank_recv)
            else:
                pickled_list_qoi = communicator.SendRecvString(string_to_send_qoi, 0, 0)
                qoi_list = pickle.loads(pickled_list_qoi.encode("latin1"))

        else:
            qoi_list = None
        return qoi_list

    def MappingAndEvaluateQuantityOfInterest(self):
        """
        function mapping the weighted pressure on reference model and calling evaluation of quantit of interest
        input:  self: an instance of the class
        """
        # call parallel fill communicator
        ParallelFillCommunicator = KratosMultiphysics.mpi.ParallelFillCommunicator(self.mapping_reference_model.GetModelPart("FluidModelPart"))
        ParallelFillCommunicator.Execute()
        # map from current model part of interest to reference model part
        mapping_parameters = KratosMultiphysics.Parameters("""{
            "mapper_type": "nearest_element",
            "interface_submodel_part_origin": "NoSlip3D_structure",
            "interface_submodel_part_destination": "NoSlip3D_structure",
            "echo_level" : 3
            }""")
        mapper = KratosMultiphysics.MappingApplication.MapperFactory.CreateMPIMapper(self._GetSolver().main_model_part,self.mapping_reference_model.GetModelPart("FluidModelPart"),mapping_parameters)
        mapper.Map(KratosMultiphysics.CONTACT_PRESSURE, \
            KratosMultiphysics.CONTACT_PRESSURE,        \
            KratosMultiphysics.MappingApplication.Mapper.FROM_NON_HISTORICAL |     \
            KratosMultiphysics.MappingApplication.Mapper.TO_NON_HISTORICAL)
        # evaluate qoi
        qoi_list = self.EvaluateQuantityOfInterest()
        return qoi_list
