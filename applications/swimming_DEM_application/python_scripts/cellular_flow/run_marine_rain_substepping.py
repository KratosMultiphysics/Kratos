import KratosSwimmingDEM as script
varying_parameters = dict()
basset_type = 2
# Nq_values = [1, 2, 4, 8, 16, 32, 64, 128]
Nq_values = [10]
Nq_values.reverse()
t_w = 1.0
m = 10
basset_force_integration_type = 2

def RunCase(varying_parameters):
    import sys
    import ProjectParameters as pp
    import DEM_explicit_solver_var as DEM_parameters
    import marine_rain_algorithm
    with marine_rain_algorithm.Algorithm(varying_parameters) as algorithm:
        try:
            test = script.Solution(algorithm)
            test.Run()
        except:
            error = sys.exc_info()
            print(error)
    del pp
    del DEM_parameters
    del marine_rain_algorithm
    del sys

varying_parameters['time_window'] = t_w
varying_parameters['basset_force_type'] = basset_type
varying_parameters['number_of_exponentials'] = m
varying_parameters['basset_force_integration_type'] = basset_force_integration_type

for Nq in Nq_values:
    varying_parameters['time_steps_per_quadrature_step'] = Nq
    RunCase(varying_parameters)
