from typing import Optional
import csv
from pathlib import Path
import numpy as np

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.SystemIdentificationApplication as KratosDT
import KratosMultiphysics.StructuralMechanicsApplication as KratosSM
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.responses.response_function import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartOperation
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy_decorator import ExecutionPolicyDecorator
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.SystemIdentificationApplication.sensor_sensitivity_solvers.system_identification_static_analysis import SystemIdentificationStaticAnalysis
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import GetSensors
from KratosMultiphysics.SystemIdentificationApplication.responses.damage_detection_response import DamageDetectionResponse

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> ResponseFunction:
    if not parameters.Has("name"):
        raise RuntimeError(f"AverageDamageDetectionResponse instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"AverageDamageDetectionResponse instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return AverageDamageDetectionResponse(parameters["name"].GetString(), model, parameters["settings"], optimization_problem)


class AverageDamageDetectionResponse(DamageDetectionResponse):
    def CalculateGradient(self, physical_variable_collective_expressions: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.CollectiveExpression]') -> None:
        # make everything zeros
        for physical_variable, collective_expression in physical_variable_collective_expressions.items():
            for container_expression in collective_expression.GetContainerExpressions():
                Kratos.Expression.LiteralExpressionIO.SetDataToZero(container_expression, physical_variable)

        # now compute sensitivities for each test scenario
        for exec_policy, sensor_measurement_data_file_name, test_case_weight in self.list_of_test_analysis_data:
            # read and replace the measurement data for each test scenario
            self.__SetSensorMeasuredValue(sensor_measurement_data_file_name)

            # run a single adjoint for each test scenario
            self.adjoint_analysis._GetSolver().GetComputingModelPart().ProcessInfo[KratosDT.TEST_ANALYSIS_NAME] = exec_policy.GetName()
            self.adjoint_analysis._GetSolver().GetComputingModelPart().ProcessInfo[Kratos.STEP] = self.optimization_problem.GetStep()
            self.adjoint_analysis.CalculateGradient(self.damage_response_function)

            for physical_variable, collective_expression in physical_variable_collective_expressions.items():
                sensitivity_variable = Kratos.KratosGlobals.GetVariable(Kratos.SensitivityUtilities.GetSensitivityVariableName(physical_variable))
                for container_expression in collective_expression.GetContainerExpressions():
                    current_gradient = container_expression.Clone()
                    self.adjoint_analysis.GetGradient(sensitivity_variable, current_gradient)
                    print("current var is ", sensitivity_variable)
                    print("current grad is ", current_gradient)
                    print("current grad is ", current_gradient.Evaluate())

                    patch_type = "multi"
                    vector_variables_independant = False

                    if patch_type == "single":
                        
                        if current_gradient.GetItemComponentCount() > 1 and vector_variables_independant:
                            print(current_gradient.GetItemComponentCount())
                            average = np.average(current_gradient.Evaluate().copy(), axis=0)
                            current_gradient_np = current_gradient.Evaluate().copy()
                            
                            rows, columns = current_gradient_np.shape
                            print("rows ", rows)
                            print("cols ", columns)
                            #print("sum ", Kratos.Expression.Utils.Sum(current_gradient.Clone()) )
                            #print("len ", len(current_gradient.Evaluate()))
                            #average = Kratos.Expression.Utils.Sum(current_gradient.Clone()) / len(current_gradient.Evaluate())
                            print("average ", average)
                            # average_np = np.zeros((rows, columns))
                            # for i in range(rows):
                            #     for j in range(columns):
                            #         average_np[i,j] = average[j]

                            #print("average_np ", average_np.shape)
                            mp = current_gradient.GetModelPart()
                            average_exp = Kratos.Expression.ElementExpression(mp)
                            print(average_exp)
                            #Kratos.Expression.CArrayExpressionIO.Read(average_exp, average_np) 
                            Kratos.Expression.LiteralExpressionIO.SetData(average_exp, average)
                            print("average_exp ", average_exp)
                            container_expression.SetExpression((container_expression.GetExpression() - average_exp.GetExpression() * test_case_weight))
                            print("updated grad is ", container_expression.Evaluate())
                            container_expression.SetExpression(Kratos.Expression.Utils.Collapse(container_expression).GetExpression())
                        else:
                            print(current_gradient.GetItemComponentCount())
                            print("sum ", Kratos.Expression.Utils.Sum(current_gradient.Clone()) )
                            print("len ", len(current_gradient.Evaluate()))
                            if sensitivity_variable == KratosSM.PRE_STRESS_SENSITIVITY:
                                print("avging prestress sens")
                                average = Kratos.Expression.Utils.Sum(current_gradient.Clone()) / (2*len(current_gradient.Evaluate()))
                            else:
                                average = Kratos.Expression.Utils.Sum(current_gradient.Clone()) / len(current_gradient.Evaluate())
                            print("average ", average)
                            container_expression.SetExpression((container_expression.GetExpression() - average * test_case_weight))
                            container_expression.SetExpression(Kratos.Expression.Utils.Collapse(container_expression).GetExpression())

                    elif patch_type == "multi":
                        if current_gradient.GetModelPart().GetSubModelPartNames():
                            sub_mp_names = current_gradient.GetModelPart().GetSubModelPartNames()
                        else: 
                            raise RuntimeError("Patch type 'multiple' is only supported for model parts with sub model parts")
                        
                        if current_gradient.GetItemComponentCount() > 1:   # PRESTRESS_VECTOR case where each component is independant
                            dummy_variable = Kratos.VELOCITY    
                        else:
                            dummy_variable = Kratos.PRESSURE

                        Kratos.Expression.VariableExpressionIO.Write(current_gradient, dummy_variable)
                            
                        for sub_mp_name in sub_mp_names:
                            sub_mp = current_gradient.GetModelPart().GetSubModelPart(sub_mp_name)    
                            exp = Kratos.Expression.ElementExpression(sub_mp)        
                            Kratos.Expression.VariableExpressionIO.Read(exp, dummy_variable)

                            if current_gradient.GetItemComponentCount() > 1 and vector_variables_independant:   # PRESTRESS_VECTOR case where each component is independant
                                avg = np.average(exp.Evaluate(), axis=0)
                                Kratos.Expression.LiteralExpressionIO.SetData(exp, avg)
                            else:
                                if current_gradient.GetItemComponentCount() > 1: # PRESTRESS_VECTOR case where components are not independant
                                    avg = np.sum(exp.Evaluate()) / (2*len(exp.Evaluate()))
                                    Kratos.Expression.LiteralExpressionIO.SetData(exp, [avg, avg, 0.0])
                                else: # scalar variable case
                                    avg = np.sum(exp.Evaluate()) / len(exp.Evaluate())
                                    Kratos.Expression.LiteralExpressionIO.SetData(exp, avg)

                            Kratos.Expression.VariableExpressionIO.Write(exp, dummy_variable)  

                        # now read back the updated gradient
                        mp_expression = Kratos.Expression.ElementExpression(current_gradient.GetModelPart())
                        Kratos.Expression.VariableExpressionIO.Read(mp_expression, dummy_variable)
                        container_expression.SetExpression((container_expression.GetExpression() - mp_expression.GetExpression() * test_case_weight))
                        container_expression.SetExpression(Kratos.Expression.Utils.Collapse(container_expression).GetExpression())
                        
                        # if current_gradient.GetItemComponentCount() > 1 and vector_variables_independant:   # PRESTRESS_VECTOR case where each component is independant
                        #     Kratos.Expression.VariableExpressionIO.Write(current_gradient, Kratos.VELOCITY, False)
                            
                        #     for sub_mp_name in sub_mp_names:
                        #         sub_mp = current_gradient.GetModelPart().GetSubModelPart(sub_mp_name)    
                        #         exp = Kratos.Expression.ElementExpression(sub_mp)        
                        #         Kratos.Expression.VariableExpressionIO.Read(exp, Kratos.VELOCITY, False)                
                        #         avg = np.average(exp.Evaluate(), axis=0)
                        #         Kratos.Expression.LiteralExpressionIO.SetData(exp, avg)
                        #         Kratos.Expression.VariableExpressionIO.Write(exp, Kratos.VELOCITY, False)                
                        # else:
                        #     if current_gradient.GetItemComponentCount() > 1: # PRESTRESS_VECTOR case where components are not independant
                        #         Kratos.Expression.VariableExpressionIO.Write(current_gradient, Kratos.VELOCITY, False)
                        #         for sub_mp_name in sub_mp_names:
                        #             sub_mp = current_gradient.GetModelPart().GetSubModelPart(sub_mp_name)    
                        #             exp = Kratos.Expression.ElementExpression(sub_mp)        
                        #             Kratos.Expression.VariableExpressionIO.Read(exp, Kratos.VELOCITY, False)                
                        #             avg = np.sum(exp.Evaluate()) / (2*len(exp.Evaluate()))
                        #             Kratos.Expression.LiteralExpressionIO.SetData(exp, Kratos.Array3(avg, avg, 0.0))
                        #             Kratos.Expression.VariableExpressionIO.Write(exp, Kratos.VELOCITY, False)
                        #     else: # scalar variable case
                        #         Kratos.Expression.VariableExpressionIO.Write(current_gradient, Kratos.WATER_PRESSURE, False)
                        #         for sub_mp_name in sub_mp_names:
                        #             sub_mp = current_gradient.GetModelPart().GetSubModelPart(sub_mp_name)    
                        #             exp = Kratos.Expression.ElementExpression(sub_mp)        
                        #             Kratos.Expression.VariableExpressionIO.Read(exp, Kratos.WATER_PRESSURE, False)                
                        #             avg = np.sum(exp.Evaluate()) / len(exp.Evaluate())
                        #             Kratos.Expression.LiteralExpressionIO.SetData(exp, avg)
                        #             Kratos.Expression.VariableExpressionIO.Write(exp, Kratos.WATER_PRESSURE, False)

                        


                        # if current_gradient.GetItemComponentCount() > 1 and vector_variables_independant:
                        #     if current_gradient.GetModelPart().GetSubModelPartNames():
                        #         sub_mp_names = current_gradient.GetModelPart().GetSubModelPartNames()
                        #         for sub_mp_name in sub_mp_names:
                        #             sub_mp = current_gradient.GetModelPart().GetSubModelPart(sub_mp_name)

                        #         mp_elements = current_gradient.GetModelPart().Elements

                        #         current_gradient_mp = current_gradient.Evaluate().copy()
                        #         averaged_current_gradient_mp = current_gradient.Evaluate().copy()*0.0

                        #         for sub_mp_name in sub_mp_names:
                        #             sub_mp = current_gradient.GetModelPart().GetSubModelPart(sub_mp_name)
                        #             sub_mp_elements = sub_mp.Elements

                        #             if current_gradient.GetItemComponentCount() > 1:
                        #                 sub_mp_current_gradient_np = np.zeros((len(sub_mp_elements), current_gradient.GetItemComponentCount()))
                        #             else: 
                        #                 sub_mp_current_gradient_np = np.zeros(len(sub_mp_elements))

                        #             for j, sub_elem in enumerate(sub_mp_elements):
                        #                 # find the index of this element in the main model part
                        #                 for i, elem in enumerate(mp_elements):
                        #                     if sub_elem.Id == elem.Id:
                        #                         if current_gradient.GetItemComponentCount() > 1:
                        #                             sub_mp_current_gradient_np[j,:] = current_gradient_mp[i,:]
                        #                         else:
                        #                             sub_mp_current_gradient_np[j] = current_gradient_mp[i]
                        #                         break

                        #             # now compute the average for this sub model part
                        #             if current_gradient.GetItemComponentCount() > 1:
                        #                 average = np.average(sub_mp_current_gradient_np.copy(), axis=0)
                        #             else:
                        #                 average = np.average(sub_mp_current_gradient_np.copy())

                        #             print("sub mp name ", sub_mp_name)
                        #             print("average ", average)

                        #             # now update the gradient for this sub model part
                        #             for j, sub_elem in enumerate(sub_mp_elements):
                        #                 # find the index of this element in the main model part
                        #                 for i, elem in enumerate(mp_elements):
                        #                     if sub_elem.Id == elem.Id:
                        #                         if current_gradient.GetItemComponentCount() > 1:
                        #                             averaged_current_gradient_mp[i,:] =  average[:] 
                        #                         else:
                        #                             averaged_current_gradient_mp[i] =  average
                        #                         break

                        #         print("updated grad is ", averaged_current_gradient_mp)

                        #         # finally set the updated gradient
                        #         averaged_current_gradient_mp_exp = Kratos.Expression.ElementExpression(current_gradient.GetModelPart())
                        #         Kratos.Expression.CArrayExpressionIO.Read(averaged_current_gradient_mp_exp, averaged_current_gradient_mp)
                        #         container_expression.SetExpression((container_expression.GetExpression() - averaged_current_gradient_mp_exp.GetExpression() * test_case_weight))
                        #         container_expression.SetExpression(Kratos.Expression.Utils.Collapse(container_expression).GetExpression())

                        #     else:
                        #         raise RuntimeError("Patch type 'multiple' is only supported for model parts with sub model parts")
                        # else:
                            # if current_gradient.GetModelPart().GetSubModelPartNames():
                            #     mp_elements = current_gradient.GetModelPart().Elements
                            #     current_gradient_mp = current_gradient.Evaluate().copy()
                            #     averaged_current_gradient_mp = current_gradient.Evaluate().copy()*0.0

                            #     sub_mp_names = current_gradient.GetModelPart().GetSubModelPartNames()
                            #     print("sub mp names ", sub_mp_names)
                            #     #raise RuntimeError(2222)
                            #     for sub_mp_name in sub_mp_names:
                            #         sub_mp = current_gradient.GetModelPart().GetSubModelPart(sub_mp_name)
                            #         sub_mp_elements = sub_mp.Elements

                            #         if current_gradient.GetItemComponentCount() > 1:
                            #             sub_mp_current_gradient_np = np.zeros((len(sub_mp_elements), current_gradient.GetItemComponentCount()))
                            #         else: 
                            #             sub_mp_current_gradient_np = np.zeros(len(sub_mp_elements))

                            #         for j, sub_elem in enumerate(sub_mp_elements):
                            #             # find the index of this element in the main model part
                            #             for i, elem in enumerate(mp_elements):
                            #                 if sub_elem.Id == elem.Id:
                            #                     if current_gradient.GetItemComponentCount() > 1:
                            #                         sub_mp_current_gradient_np[j,:] = current_gradient_mp[i,:]
                            #                     else:
                            #                         sub_mp_current_gradient_np[j] = current_gradient_mp[i]
                            #                     break

                            #         # now compute the average for this sub model part
                            #         if current_gradient.GetItemComponentCount() > 1:
                            #             average = np.sum(sub_mp_current_gradient_np.copy()) / (2*len(sub_mp_current_gradient_np))
                            #         else:
                            #             average = np.sum(sub_mp_current_gradient_np.copy()) / len(sub_mp_current_gradient_np)

                            #         print("sub mp name ", sub_mp_name)
                            #         print("average ", average)

                            #         # now update the gradient for this sub model part
                            #         for j, sub_elem in enumerate(sub_mp_elements):
                            #             # find the index of this element in the main model part
                            #             for i, elem in enumerate(mp_elements):
                            #                 if sub_elem.Id == elem.Id:
                            #                     if current_gradient.GetItemComponentCount() > 1:
                            #                         averaged_current_gradient_mp[i,:] =  average 
                            #                     else:
                            #                         averaged_current_gradient_mp[i] = average 
                            #                     break

                                    

                            #     print("updated grad is ", averaged_current_gradient_mp)

                            #     # finally set the updated gradient
                            #     averaged_current_gradient_mp_exp = Kratos.Expression.ElementExpression(current_gradient.GetModelPart())
                            #     Kratos.Expression.CArrayExpressionIO.Read(averaged_current_gradient_mp_exp, averaged_current_gradient_mp)
                            #     container_expression.SetExpression((container_expression.GetExpression() - averaged_current_gradient_mp_exp.GetExpression() * test_case_weight))
                            #     container_expression.SetExpression(Kratos.Expression.Utils.Collapse(container_expression).GetExpression())

                            # else:
                            #     raise RuntimeError("Patch type 'multiple' is only supported for model parts with sub model parts")

            #raise RuntimeError(10000)

    def __GetSensor(self, sensor_name: str) -> KratosDT.Sensors.Sensor:
        return self.sensor_name_dict[sensor_name]

    def __GetHeaderIndices(self, csv_stream: csv.reader) -> 'tuple[int, int]':
        headers = [s.strip() for s in next(csv_stream)]
        name_index = headers.index("name")
        value_index = headers.index("value")
        return name_index, value_index
    
    def __SetSensorMeasuredValue(self, sensor_measurement_data_file_name: str) -> None:
        with open(sensor_measurement_data_file_name, "r") as csv_measurement_file:
            csv_measurement_stream = csv.reader(csv_measurement_file, delimiter=",")
            measured_name_index, measured_value_index = self.__GetHeaderIndices(csv_measurement_stream)

            for measured_row in csv_measurement_stream:
                measured_sensor_name = measured_row[measured_name_index].strip()
                measured_value = float(measured_row[measured_value_index])
                self.__GetSensor(measured_sensor_name).GetNode().SetValue(KratosDT.SENSOR_MEASURED_VALUE, measured_value)

    def _GetResponsePrefix(self) -> str:
        return "AverageDamageDetectionResponse"
    
    def __str__(self) -> str:
        return f"Response [type = {self._GetResponsePrefix()}, name = {self.GetName()}, model part name = {self.model_part.FullName()}]"