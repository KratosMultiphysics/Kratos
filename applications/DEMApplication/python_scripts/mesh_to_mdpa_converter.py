# MESH TO COHESIVE MDPA CONVERTER
import sys

mesh_namefile = sys.argv[1] + '.msh'
mdpa_namefile = sys.argv[2] + '.mdpa'

SpheresMesh = open(mesh_namefile, 'r')
SpheresMdpa = open(mdpa_namefile, 'w')

top = 0.0
bottom = 0.0
radius = 0.0
Triaxial = False

# Test type
if Triaxial:
    top = 0.002771549
    bottom = 0.0
    radius = 0.0006907713
else:
    top = 0.0008956682
    bottom = 0.0
    radius = 0.001222598

node_section_started = False
node_section_finished = False
element_section_started = False
node_list = []
coord_x_list = []
coord_y_list = []
coord_z_list = []
element_list = []
radius_list = []
props_list = []
zeros_list = []

for Line in SpheresMesh:

    if Line.startswith('Coordinates'):
        node_section_started = True
        continue
    
    if node_section_started and not node_section_finished:
        if Line.startswith('End Coordinates'):
            node_section_finished = True
            continue
        data = Line.split(" ")
        data[0] = int(data[0])
        data[1] = float(data[1])
        data[2] = float(data[2])
        data[3] = float(data[3])
        
        if Triaxial:
            if data[2] > top:
                continue
            if data[2] < bottom:
                continue
            if (data[1] * data[1] + data[3] * data[3] > radius * radius):
                continue
        else:
            if data[3] > top:
                continue
            if data[3] < bottom:
                continue
            if (data[1] * data[1] + data[2] * data[2] > radius * radius):
                continue
                
        node_list.append(data[0])
        element_list.append(data[0])
        coord_x_list.append(data[1])
        coord_y_list.append(data[2])
        coord_z_list.append(data[3])
        zeros_list.append('0')

    if Line.startswith('Elements'):
        element_section_started = True
        continue
    
    if element_section_started: # and not element_section_finished:
        if Line.startswith('End Elements'):
            break
        data = Line.split(" ")
        data[0] = int(data[0])
        data[2] = float(data[2])
        data[3] = int(data[3])
        if data[0] in element_list:
            radius_list.append(data[2])
            props_list.append(data[3])

zeros_list = [int(i) for i in zeros_list]

SpheresMesh.close()

SpheresMdpa.write('''Begin ModelPartData
//  VARIABLE_NAME value
End ModelPartData\n
Begin Properties 1
PARTICLE_DENSITY 2650.0
YOUNG_MODULUS 8.03e7
POISSON_RATIO 0.25
FRICTION 0.6
PARTICLE_COHESION 0.0
COEFFICIENT_OF_RESTITUTION 0.01
PARTICLE_MATERIAL 1
ROLLING_FRICTION 0.025
ROLLING_FRICTION_WITH_WALLS 0.0
DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME DEM_KDEM_with_damage_parallel_bond_Hertz //_Hertz
INTERNAL_COHESION 13.72
INTERNAL_FRICTION_ANGLE 24.1
DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME DEM_D_Linear_HighStiffness //_Coulomb //DEM_D_Linear_HighStiffness
CONTACT_TAU_ZERO 10e4
CONTACT_SIGMA_MIN 13.72e6
CONTACT_INTERNAL_FRICC 0.0
ROTATIONAL_MOMENT_COEFFICIENT 0.01
KDEM_STANDARD_DEVIATION_TAU_ZERO 3.43
KDEM_STANDARD_DEVIATION_FRICTION 0.2235
SHEAR_ENERGY_COEF 4.0
LOOSE_MATERIAL_YOUNG_MODULUS 8e7
FRACTURE_ENERGY 0.0
End Properties\n
Begin Nodes\n''')

for i in range(len(node_list)):
    SpheresMdpa.write("%i %12.8f %12.8f %12.8f\n" % (node_list[i], coord_x_list[i], coord_y_list[i], coord_z_list[i]))
    
SpheresMdpa.write('''End Nodes\n
Begin Elements SphericContinuumParticle3D\n''')

for i in range(len(node_list)):
    SpheresMdpa.write("%i %i %i\n" % (element_list[i], props_list[i], node_list[i]))
    
SpheresMdpa.write('''End Elements\n
Begin NodalData RADIUS\n''')

for i in range(len(node_list)):
    SpheresMdpa.write("%i %i %12.8f\n" % (node_list[i], zeros_list[i], radius_list[i]))
    
SpheresMdpa.write('''End NodalData\n
Begin NodalData COHESIVE_GROUP\n''')
    
for i in range(len(node_list)):
    SpheresMdpa.write("%i %i %i\n" % (node_list[i], zeros_list[i], props_list[i]))
    
SpheresMdpa.write('''End NodalData\n
Begin NodalData SKIN_SPHERE
End NodalData\n
Begin SubModelPart PartsCont_dem // Group dem // Subtree PartsCont
    Begin SubModelPartNodes\n''')

for i in range(len(node_list)):
    SpheresMdpa.write("%i\n" % (node_list[i]))
    
SpheresMdpa.write('''End SubModelPartNodes
    Begin SubModelPartElements\n''')

for i in range(len(node_list)):
    SpheresMdpa.write("%i\n" % (element_list[i]))

SpheresMdpa.write('''End SubModelPartElements
    Begin SubModelPartConditions
    End SubModelPartConditions
End SubModelPart\n''')

SpheresMdpa.close()
