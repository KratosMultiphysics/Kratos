import KratosSwimmingDEM as script
import sys
import ProjectParameters as pp
import DEM_explicit_solver_var as DEM_parameters
varying_parameters = dict()

irregular_mesh_sizes = set() #{0.1, 0.2, 0.4}
regular_mesh_n_points = {10}
derivatives_types = {6}
combinations_that_failed = []
for size in irregular_mesh_sizes.union(regular_mesh_n_points):
    varying_parameters['size_parameter'] = size
    for derivatives_type in derivatives_types:
        varying_parameters['material_acceleration_calculation_type'] = derivatives_type
        varying_parameters['laplacian_calculation_type'] = derivatives_type
        import pouliot_benchmark_2D_analysis
        with pouliot_benchmark_2D_analysis.Algorithm(varying_parameters) as algorithm:
            try:
                test = script.Solution(algorithm, varying_parameters)
                test.Run()
            except:
                error = sys.exc_info()
                combinations_that_failed.append({'size':size, 'type':derivatives_type})

print()
print('****************************************')

if len(combinations_that_failed):
    print('The following combinations produced an error:')
    print()
    for combination in combinations_that_failed:
        print(combination)
        print(error)
else:
    print('All combinations run without errors')
print('****************************************')
