domain_size = 2

bulk_modulus = -5.00000e+03
density = 2.50000e+03
viscosity = 6.50000e-01
blow_pressure = 2.0000e+05
#Dt = 2.00000e-03
#the below one is for the measured temperature and thus very high viscosity
Dt=0.005
#max_time = 0.4500000e+00
max_time = 0.6000000e+00
output_step = 2.00000e-02
alpha_shape = 2.00000e+00
erase_nodes = 1.00000e+00
adaptive_refinement = 0.00000e+00
delete_nodes_close_to_wall = 1.00000e+00
bounding_box_corner1_x = -1.00000e-04
bounding_box_corner1_y = -1.00000e-04
bounding_box_corner1_z = -1.00000e+00
bounding_box_corner2_x = 1.00000e+00
bounding_box_corner2_y = 1.00000e+00
bounding_box_corner2_z = 1.00000e+00
SolverType = "Quasi_Inc_Constant_Pressure"
# Declare Python Variables

problem_name = r'geom_corrected_geom_2_STAGES'
problem_path = r'/home/pavel/Kratos/EXAMPLES/BOTTLE_FORMING/geom_corrected_geom_2_STAGES.gid'
kratos_path =  r'/home/pavel'

