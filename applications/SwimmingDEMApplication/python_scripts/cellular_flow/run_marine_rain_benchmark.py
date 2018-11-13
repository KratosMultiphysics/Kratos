import KratosSwimmingDEM as script
varying_parameters = dict()
types = {4}
m_values = {1, 2, 4, 5}
t_w = 0.1

def RunCase(varying_parameters):
    import sys
    import ProjectParameters as pp
    import DEM_explicit_solver_var as DEM_parameters
    import marine_rain_analysis
    with marine_rain_analysis.Algorithm(varying_parameters) as algorithm:
        try:
            test = script.Solution(algorithm)
            test.Run()
        except:
            error = sys.exc_info()
            print(error)
    del pp
    del DEM_parameters
    del marine_rain_analysis
    del sys

for basset_type in types:
    varying_parameters['time_window'] = t_w
    varying_parameters['basset_force_type'] = basset_type
    if basset_type == 4:
        for m in m_values:
            varying_parameters['number_of_exponentials'] = m
            RunCase(varying_parameters)
    else:
        RunCase(varying_parameters)
