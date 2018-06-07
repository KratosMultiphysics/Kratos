from __future__ import absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics as Kratos
import KratosMultiphysics.FluidDynamicsApplication
import KratosMultiphysics.AdjointFluidApplication
from fluid_dynamics_analysis import FluidDynamicsAnalysis
from adjoint_fluid_analysis import AdjointFluidAnalysis


def Run(optimization_model_part, response_data):
    primal_parameter_file_name = "ProjectParametersPrimal.json"
    adjoint_parameter_file_name = "ProjectParametersAdjoint.json"

    parameters = Kratos.Parameters(r'''{}''')

    with open(primal_parameter_file_name,'r') as primal_parameter_file:
        parameters.AddValue("primal_settings",Kratos.Parameters(primal_parameter_file.read()))

    with open(adjoint_parameter_file_name,'r') as adjoint_parameter_file:
        parameters.AddValue("adjoint_settings", Kratos.Parameters(adjoint_parameter_file.read()))
    
    model = Kratos.Model()

    print("> Running primal simulation")
    primal_simulation = FluidDynamicsAnalysis(model,parameters["primal_settings"])
    primal_simulation.Run()

    print("> Running adjoint simulation")
    adjoint_model = Kratos.Model()
    adjoint_simulation = AdjointFluidAnalysis(adjoint_model,parameters["adjoint_settings"])
    adjoint_simulation.Run()

    problem_data = parameters["adjoint_settings"]["problem_data"]
    end_time = problem_data["start_step"].GetDouble()
    nsteps = problem_data["nsteps"].GetInt()
    time_step = parameters["adjoint_settings"]["solver_settings"]["time_stepping"]["time_step"].GetDouble()

    start_time = end_time + time_step*float(nsteps)

    response_data["response_value"] = __calculate_drag("NoSlip2D_Cylinder.drag", start_time, end_time)
    response_data["response_gradients"] = __get_response_gradients(adjoint_simulation)

def __get_response_gradients(adjoint_simulation):
    gradient = {}
    for node in adjoint_simulation.main_model_part.Nodes:
        gradient[node.Id] = node.GetSolutionStepValue(KratosMultiphysics.SHAPE_SENSITIVITY)
    return gradient

def __calculate_drag(filename, start_time, end_time):
    with open(filename, "r") as file_input:
        lines = file_input.readlines()
    file_input.close()
    
    found_headers = False
    total_drag = 0.0
    current_time = 0.0
    previous_time = 0.0
    
    for line in lines:
        _data = line.strip().split()
        if not found_headers:
            if len(line)>4:
                if line[:4]=="Time":
                    found_headers = True
                    continue
        
        if line.strip()=="":
            continue
        
        if found_headers:
            _time = float(_data[0])
            _fx = float(_data[1])
            _fy = float(_data[2])
            _fz = float(_data[3])
            if (_time > end_time):
                break
            
            if (_time >= start_time):
                previous_time = current_time
                current_time = _time
                total_drag += _fx

    total_time = (current_time-start_time)
    delta_t = current_time - previous_time
    
    total_drag /= (total_time/delta_t)
    
    return total_drag  
    
    
if __name__ == "__main__":
    response_values = {}
    Run(None, response_values) 

    print(response_values)           

    
