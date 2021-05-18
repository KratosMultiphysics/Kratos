import sys
import math

Mailfile = sys.argv[1]

in_file = open(Mailfile, 'r')

number_of_sectors = 12
span = 360.0 / number_of_sectors
accumulation = [0] * number_of_sectors
fraction = [0] * number_of_sectors
fraction_exp_300 = [0.0, 0.005645342312008972, 0.030788274905921897, 0.2692645408331682, 0.3890479963028982, 0.19824427279329238, 0.06181699346405223, 0.006438106555753553, 0.0002079619726678139, 0.0, 0.0, 0.0]
fraction_exp_500 = [0.0, 0.00490797546012256, 0.009815950920245342, 0.04539877300613493, 0.40981595092024536, 0.44417177914110423, 0.08036809815950907, 0.00490797546012256, 0.0, 0.0, 0.0, 0.0]
fraction_exp_650 = [0.0, 0.0, 0.006341436575547199, 0.047897039000100894, 0.47890931765871014, 0.40405553364939345, 0.015909436868987936, 0.0, 0.0, 0.0, 0.0, 0.0]
total_accumulation = 0

for line in in_file.readlines():
    numbers = line.split()
    x_coor = float(numbers[1])
    y_coor = float(numbers[2])
    
    angle_in_radians = math.atan2(y_coor,x_coor)
    angle_in_degrees = angle_in_radians * 180 / math.pi

    #CORRECTION
    angle_in_degrees *= -1.0 #Due to the orientation chosen by Liederkerke (clockwise in his paper)
    
    if angle_in_degrees < 0.0:
        angle_in_degrees += 360.0        
    
    for sector in range(0,number_of_sectors):
        if angle_in_degrees <= (sector +1) * span:
            final_sector = sector
            accumulation[sector] += 1
            total_accumulation += 1
            break
      
print("Particle fractions:")
for sector in range(0,number_of_sectors):
    fraction[sector] = float( accumulation[sector] ) / float( total_accumulation )
    print(fraction[sector])
print("Total number of particles was: " + str(total_accumulation))
        
in_file.close()

import matplotlib.pyplot as plotter
x_range = range(0,number_of_sectors*30,30)
plotter.plot(x_range, fraction, label="DEM")
plotter.plot(x_range, fraction_exp_300,label="experiment, 300 rpm")
plotter.plot(x_range, fraction_exp_500,label="experiment, 500 rpm")
plotter.plot(x_range, fraction_exp_650,label="experiment, 650 rpm")
plotter.xlabel('Compartment angle')
plotter.ylabel('Mass fraction')
plotter.legend()
plotter.show()
