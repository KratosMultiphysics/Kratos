from typing import Optional
import numpy as np
import scipy.sparse.linalg as spla
import scipy.linalg as sla

import KratosMultiphysics as Kratos
import KratosMultiphysics.SystemIdentificationApplication as KratosSI
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.responses.response_function import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartOperation
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import GetMaskStatusControllers

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> ResponseFunction:
    if not parameters.Has("name"):
        raise RuntimeError(f"SensorMACResponse instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"SensorMACResponse instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return SensorMACResponse(parameters["name"].GetString(), model, parameters["settings"], optimization_problem)


class SensorMACResponse(ResponseFunction):
    def __init__(self, name: str, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        super().__init__(name)

        default_settings = Kratos.Parameters("""{
            "sensor_group_name"      : "",
            "sensor_mask_name"       : "",
            "number_of_eigen_values" : 5,
            "penalty_factor"         : 1.0,
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
        self.penalty_factor = parameters["penalty_factor"].GetDouble()

        self.model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"response_{self.GetName()}", [self.sensor_group_name], False)
        self.model_part: Optional[Kratos.ModelPart] = None
        self.optimization_problem = optimization_problem
        self.mask_status: 'Optional[KratosSI.SensorMaskStatus]' = None

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
        self.A = self.G_T @ self.G_T.T
        self.reference_eigen_values, self.reference_eigen_vectors = spla.eigsh(self.A, k=self.number_of_eigen_values, which='LM')

    def Check(self) -> None:
        pass

    def Finalize(self) -> None:
        pass

    def GetInfluencingModelPart(self) -> Kratos.ModelPart:
        if self.model_part is None:
            raise RuntimeError("Please call SensorMACResponse::Initialize first.")
        return self.model_part

    def CalculateValue(self) -> float:
        theta = np.array(self.mask_status.GetSensorStatus())
        G_T = np.array(self.mask_status.GetMaskStatuses()) # (n,m) matrix where n is number of elements, m is number of masks
        self.A = G_T @ G_T.T
        self.eigen_values, self.eigen_vectors = spla.eigsh(self.A, k=self.number_of_eigen_values, which='LM')

        # calculate the MAC matrix
        self.mac = np.empty((self.number_of_eigen_values, self.number_of_eigen_values))
        for i in range(self.number_of_eigen_values):
            for j in range(self.number_of_eigen_values):
                self.mac[i, j] = (self.reference_eigen_vectors[:, i].dot(self.eigen_vectors[:, j])) ** 2

        return 0.25 * np.linalg.norm(self.mac - np.identity(self.number_of_eigen_values)) ** 2 + self.penalty_factor * theta.dot(1 - theta)

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

            # d(omega_i) / d (v_a[c])
            data = combined_ta.GetTensorAdaptors()[0].data[:]
            data[:] = 0.0

            signed_mac = np.empty((self.number_of_eigen_values, self.number_of_eigen_values))
            for i in range(self.number_of_eigen_values):
                for j in range(self.number_of_eigen_values):
                    signed_mac[i, j] = (self.reference_eigen_vectors[:, i].dot(self.eigen_vectors[:, j]))

            identity = np.identity(self.number_of_eigen_values)

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

                # columns (sensor theta); rows (eigen vector component sensitivitis)
                d_v_i_d_theta = w_all - np.outer(v_i, coeffs)

                for j in range(self.number_of_eigen_values):
                    data[:] += (self.mac[i,j] - identity[i, j]) * (self.reference_eigen_vectors[:, j].dot(self.eigen_vectors[:, i]) * (self.reference_eigen_vectors[:, j].dot(d_v_i_d_theta)))

            data[:] += self.penalty_factor * ( 1 - 2.0 * sensor_status)
            Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(combined_ta, perform_collect_data_recursively=False, copy=False).CollectData()

    def __str__(self) -> str:
        return f"Response [type = {self.__class__.__name__}, name = {self.GetName()}, model part name = {self.model_part.FullName()}]"