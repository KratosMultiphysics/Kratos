model_part_file_name = 'circle_plate.mdpa'
ref_model_part_file_name = 'ref_circle_plate.mdpa'

import sys
plate_radius = float(sys.argv[1])

# read ref model part
with open(ref_model_part_file_name, 'r') as model_part_file:
    lines = model_part_file.readlines()

# scale node coordinates by radius
node_lines = lines[lines.index('Begin Nodes\n')+1:lines.index('End Nodes\n')]
for line in node_lines:
    components = line.split()
    node_id = int(components[0])
    x = float(components[1])
    y = float(components[2])
    z = float(components[3])
    new_line = '{:5d}'.format(node_id) + ' ' \
               + '{:19.10f}'.format(x * plate_radius) + ' ' \
               + '{:19.10f}'.format(y * plate_radius) + ' ' \
               + '{:19.10f}'.format(z * plate_radius) + '\n'
    lines[lines.index(line)] = new_line

with open(model_part_file_name, 'w') as model_part_file:
    model_part_file.writelines(lines)
