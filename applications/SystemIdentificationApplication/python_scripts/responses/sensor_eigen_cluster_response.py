from typing import Optional
import numpy as np
import scipy.sparse.linalg as spla
import scipy.linalg as sla

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.SystemIdentificationApplication as KratosSI
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.responses.response_function import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartOperation
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import GetMaskStatusControllers

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> ResponseFunction:
    if not parameters.Has("name"):
        raise RuntimeError(f"SensorEigenClusterResponse instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"SensorEigenClusterResponse instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return SensorEigenClusterResponse(parameters["name"].GetString(), model, parameters["settings"], optimization_problem)


class SensorEigenClusterResponse(ResponseFunction):
    def __init__(self, name: str, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        super().__init__(name)

        default_settings = Kratos.Parameters("""{
            "sensor_group_name"      : "",
            "sensor_mask_name"       : "",
            "number_of_eigen_values" : 5,
            "boltzmann_operator_beta": 30.0,
            "store_eigen_vectors"    : true,
            "echo_level"             : 0
        }""")

        parameters.ValidateAndAssignDefaults(default_settings)


        self.model = model
        self.sensor_group_name = parameters["sensor_group_name"].GetString()
        self.sensor_mask_name = parameters["sensor_mask_name"].GetString()
        self.number_of_eigen_values = parameters["number_of_eigen_values"].GetInt()
        self.store_eigen_vectors = parameters["store_eigen_vectors"].GetBool()
        self.echo_level = parameters["echo_level"].GetInt()

        self.model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"response_{self.GetName()}", [self.sensor_group_name], False)
        self.model_part: Optional[Kratos.ModelPart] = None
        self.optimization_problem = optimization_problem
        self.mask_status: 'Optional[KratosSI.SensorMaskStatus]' = None

        self.boltzmann_operator = KratosOA.BoltzmannOperator(parameters["boltzmann_operator_beta"].GetDouble())

    def GetImplementedPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return [KratosSI.SENSOR_STATUS]

    def Initialize(self) -> None:
        self.model_part = self.model_part_operation.GetModelPart()
        sensor_group_data = ComponentDataView(self.sensor_group_name, self.optimization_problem)

        for controller in GetMaskStatusControllers(sensor_group_data, self.sensor_mask_name):
            if isinstance(controller, KratosSI.SensorMaskStatus):
                self.mask_status = controller

        if self.mask_status == None:
            raise RuntimeError(f"SensorMaskStatus controller not found for the sensor_mask_name = \"{self.sensor_mask_name}\" in the sensor_group = \"{self.sensor_group_name}\".")

        if self.mask_status.GetSensorModelPart().NumberOfNodes() == 0:
            raise RuntimeError(f"No sensors found in the sensor group \"{self.sensor_group_name}\".")

        self.G_T = np.array(self.mask_status.GetMasks()) # (n,m) matrix where n is number of elements, m is number of masks
        self.domain_size_ta = Kratos.TensorAdaptors.GeometryMetricsTensorAdaptor(self.mask_status.GetMaskContainer(), Kratos.TensorAdaptors.GeometryMetricsTensorAdaptor.Metric.DomainSize)
        self.domain_size_ta.CollectData()
        self.total_domain_size = np.sum(self.domain_size_ta.data)

    def Check(self) -> None:
        pass

    def Finalize(self) -> None:
        pass

    def GetInfluencingModelPart(self) -> Kratos.ModelPart:
        if self.model_part is None:
            raise RuntimeError("Please call SensorEigenClusterResponse::Initialize first.")
        return self.model_part

    def CalculateValue(self) -> float:
        sensor_group_data = ComponentDataView(self.sensor_group_name, self.optimization_problem)
        ta = Kratos.TensorAdaptors.VariableTensorAdaptor(self.mask_status.GetMaskContainer(), KratosSI.SENSOR_STATUS)

        G_T = np.array(self.mask_status.GetMaskStatuses()) # (n,m) matrix where n is number of elements, m is number of masks
        self.A = G_T @ G_T.T
        self.eigen_values, self.eigen_vectors = spla.eigsh(self.A, k=self.number_of_eigen_values, which='LM')

        self.np_cluster_size_ratios = np.empty(self.number_of_eigen_values, dtype=np.float64)
        self.nd_data_cluster_size_ratios = Kratos.DoubleNDData(self.np_cluster_size_ratios, copy=False)
        for i in range(self.number_of_eigen_values):
            self.np_cluster_size_ratios[i] = np.sum(self.domain_size_ta.data[:] * self.eigen_vectors[:, i] ** 2) ** 2 / np.sum(self.domain_size_ta.data[:] * self.eigen_vectors[:, i] ** 4) / self.total_domain_size

            if self.store_eigen_vectors:
                ta.data[:] = self.eigen_vectors[:, i]
                sensor_group_data.GetUnBufferedData().SetValue(f"Eigen_Vector_{i}", ta.Clone(), overwrite=True)

        self.boltzmann_operator.Update(self.nd_data_cluster_size_ratios)

        if (self.echo_level > 0):
            Kratos.Logger.PrintInfo(self.__class__.__name__, f"Current Eigen values and respective cluster sizes:\n\t" + "\n\t".join([f"{i}: {v1:0.6e} - {v2:0.6e}" for i, (v1, v2) in enumerate(zip(self.eigen_values, self.np_cluster_size_ratios))]))

        return self.boltzmann_operator.CalculateValue()

    def CalculateGradient(self, physical_variable_combined_tensor_adaptor: 'dict[SupportedSensitivityFieldVariableTypes, Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor]') -> None:
        for physical_variable, combined_ta in physical_variable_combined_tensor_adaptor.items():
            if physical_variable != KratosSI.SENSOR_STATUS:
                raise RuntimeError(f"Unsupported variable = {physical_variable.Name()}.")

            if len(combined_ta.GetTensorAdaptors()) != 1:
                raise RuntimeError(f"Currently only supports sensitivities for one model part.")

            if combined_ta.GetTensorAdaptors()[0].GetContainer() != self.mask_status.GetSensorModelPart().Nodes:
                raise RuntimeError("Tensor adaptor container and mask status container mismatch.")

            sensor_status = np.array(self.mask_status.GetSensorStatus())
            alpha = 2.0 * sensor_status

            d_J_d_omega = self.boltzmann_operator.CalculateGradient().to_numpy()

            # d(omega_i) / d (v_a[c])
            data = combined_ta.GetTensorAdaptors()[0].data[:]
            data[:] = 0.0
            for i in range(self.number_of_eigen_values):
                lam_i = self.eigen_values[i]
                v_i = self.eigen_vectors[:, i]

                # Build and solve Nelson system directly for d(v_i)/d(theta) for all sensors at once.
                K_mod = self.A - lam_i * np.eye(self.A.shape[0], dtype=np.float64)
                max_idx = np.argmax(np.abs(v_i))
                K_mod[max_idx, :] = 0.0
                K_mod[:, max_idx] = 0.0
                K_mod[max_idx, max_idx] = 1.0
                lu_and_piv = sla.lu_factor(K_mod, check_finite=False)

                projection = self.G_T.T @ v_i
                d_lam_d_theta = alpha * projection ** 2

                rhs_all = np.multiply.outer(v_i, d_lam_d_theta)
                rhs_all -= self.G_T * (alpha * projection)[None, :]
                rhs_all[max_idx, :] = 0.0

                w_all = sla.lu_solve(lu_and_piv, rhs_all, check_finite=False)
                coeffs = v_i @ w_all
                d_v_i_d_theta = w_all - np.outer(v_i, coeffs)

                a = np.sum(self.domain_size_ta.data[:] * self.eigen_vectors[:, i] ** 4)
                b = np.sum(self.domain_size_ta.data[:] * self.eigen_vectors[:, i] ** 2)
                d_omega_d_v_a = 4.0 * (a * b * self.domain_size_ta.data[:] * self.eigen_vectors[:, i] - b ** 2 * (self.domain_size_ta.data[:] * self.eigen_vectors[:, i] ** 3)) / (self.total_domain_size * a ** 2)

                # d(omega_i) / d(theta) = d(omega_i) / d (v_a[c]) * d(v_a[c]) / d (theta)
                d_omega_i_d_theta = d_omega_d_v_a @ d_v_i_d_theta
                data[:] += d_omega_i_d_theta * d_J_d_omega[i]

            Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(combined_ta, perform_collect_data_recursively=False, copy=False).CollectData()

    def __str__(self) -> str:
        return f"Response [type = {self.__class__.__name__}, name = {self.GetName()}, model part name = {self.model_part.FullName()}]"