from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import csv
from KratosMultiphysics import *


# SectionType: String containing the type of the section: "Square", "IPN", "HEB", ...
# SectionData: A list containing the data related to section: IPN -> [100], Square -> [h,b] ...
def SetProperties(SectionType, SectionData, BeamProperties):

    import os
    csvFile = os.path.dirname(__file__) + '/Profiles.csv'

    if(SectionType == "IPN") or (SectionType == "IPE") or (SectionType == "HEB") or (SectionType == "HEA") or (SectionType == "HEM") or (SectionType == "UPN"):
        if (len(SectionData) < 1):
            print("Error, Section needs at least the size of section to be given in SectionData")
            raise
        size_profile = SectionData[0]
        SectionProperties = searchCVSValues(csvFile, SectionType, size_profile)
 #       print SectionProperties
        inertia = Matrix(2, 2)
        inertia[0, 0] = float(SectionProperties["Iz(m4)"])
        inertia[0, 1] = float(SectionProperties["Iz(m4)"])  # we have to set this correctly
        inertia[1, 0] = float(SectionProperties["Iz(m4)"])  # we have to set this correctly
        inertia[1, 1] = float(SectionProperties["Iy(m4)"])
        print(("inertia", inertia))
        BeamProperties.SetValue(LOCAL_INERTIA_TENSOR, inertia)
        return BeamProperties

    if(SectionType == "Circular"):
        if (len(SectionData) < 1):
            print("Error, Section needs at least the size of section to be given in SectionData")
            raise

        radius = SectionData[0]

        shape = 'Circular'

        circular_area = 3.14 * (radius ** 2)
        circular_inertia = (3.14 * (radius ** 4)) / 4
        circular_inertia_polar = (3.14 * (radius ** 4)) / 4
        circular_module = (3.14 * (radius ** 3)) / 4
        circular_turning_radius = (((3.14 * (radius ** 4)) / 4.0) / (3.14 * (radius**2)))**(1.0/2.0)
        inertia = Matrix(2, 2)
        inertia[0, 0] = circular_inertia
        inertia[0, 1] = circular_inertia_polar
        inertia[1, 0] = circular_inertia_polar
        inertia[1, 1] = circular_inertia
        print(("inertia", inertia))
        BeamProperties.SetValue(LOCAL_INERTIA_TENSOR, inertia)

    if(SectionType == "Square"):
        if (len(SectionData) < 1):
            print("Error, Section needs at least the size of section to be given in SectionData")
            raise

        height_square = SectionData[0]
        base_square = SectionData[1]

        shape = 'Square'

        square_area = base_square * height_square
        square_inertia_z = (base_square * height_square ** 3) / 12.0
        square_inertia_y = (height_square * base_square ** 3) / 12.0
        square_inertia_polar = ((base_square * height_square ** 3) / 12.0) + ((height_square * base_square**3)/12.0)
        square_module_z = (base_square * height_square ** 2) / 6.0
        square_module_y = (height_square * base_square ** 2) / 6.0
        square_turning_radius_z = ((((base_square * height_square ** 3) / 12.0) / (base_square * height_square)) ** (1.0/2.0))
        square_turning_radius_y = ((((height_square * base_square ** 3) / 12.0) / (base_square * height_square)) ** (1.0/2.0))
        inertia = Matrix(2, 2)
        inertia[0, 0] = square_inertia_z
        inertia[0, 1] = square_inertia_polar
        inertia[1, 0] = square_inertia_polar
        inertia[1, 1] = square_inertia_y
        print(("inertia", inertia))
        BeamProperties.SetValue(LOCAL_INERTIA_TENSOR, inertia)


def searchCVSValues(fileName, shape, size):
    f = open(fileName, 'rU')
    f.seek(0)
    reader = csv.reader(f, delimiter=';')
    headers = next(reader)
    result = [line for line in reader if line[0] == shape and line[1] == size]
    print(headers)
    print(result)
    result = result[0]
    result = dict(list(zip(headers, result)))
    for k in list(result.keys()):
        print('  ', k, ' : ', result[k])
    return result





# x = int(raw_input("Select shape (1=Circular,2=Square,3=Profiles,4=Irregular): "))
#
# if x == 1:
# shape = 'Circular'
# print 'The selected shape is: ' + shape
# radius = float(raw_input("Radius (in m): "))
# circular_area = 3.14*(radius**2)
# circular_inertia = (3.14*(radius**4))/4
# circular_inertia_polar = (3.14*(radius**4))/4
# circular_module = (3.14*(radius**3))/4
# circular_turning_radius = (((3.14*(radius**4))/4.0)/(3.14*(radius**2)))**(1.0/2.0)
#
# print 'The area is: '+ str(circular_area)+ 'm2'
# print 'Inertia is: '+ str(circular_inertia)+ 'm4'
# print 'Inertia_polar is: '+ str(circular_inertia_polar)+ 'm4'
# print 'Module is: '+ str(circular_module)+ 'm3'
# print 'Turning radius is: '+ str(circular_turning_radius)+ 'm'
#
#
# elif x == 2:
# shape = 'Square'
# print 'The selected shape is: ' + shape
# base_square = float(raw_input("Base (in m): "))
# height_square = float(raw_input("Height: (in m): "))
# square_area = base_square * height_square
# square_inertia_z = (base_square*height_square**3)/12.0
# square_inertia_y = (height_square*base_square**3)/12.0
# square_inertia_polar = ((base_square*height_square**3)/12.0)+((height_square*base_square**3)/12.0)
# square_module_z = (base_square*height_square**2)/6.0
# square_module_y = (height_square*base_square**2)/6.0
# square_turning_radius_z = ((((base_square*height_square**3)/12.0)/(base_square * height_square))**(1.0/2.0))
# square_turning_radius_y = ((((height_square*base_square**3)/12.0)/(base_square * height_square))**(1.0/2.0))
#
#
# print 'The area is: '+ str(square_area)+ 'm2'
# print 'Inertia_z is: '+ str(square_inertia_z)+ 'm4'
# print 'Inertia_y is: '+ str(square_inertia_y)+ 'm4'
# print 'Inertia_polar is: '+ str(square_inertia_polar)+ 'm4'
# print 'Module_z is: '+ str(square_module_z)+ 'm3'
# print 'Module_y is: '+ str(square_module_y)+ 'm3'
# print 'Turning radius_z is: '+ str(square_turning_radius_z)+ 'm'
# print 'Turning radius_y is: '+ str(square_turning_radius_y)+ 'm'
#
# elif x == 3:
# shape = 'Profiles'
# print 'The selected shape is: ' + shape
#
# type_profile = raw_input("Insert profile (IPN,IPE,HEB,HEA,HEM,UPN): ")
# size_profile = raw_input("Insert the profile size: ")
#
# valors = searchCVSValues(csvFile,type_profile,size_profile)
#
#
#
# elif x == 4:
# forma = 'Irregular'
# print 'La forma seleccionada es: ' + forma
#
#
#
#
# z = raw_input("Salir...")
