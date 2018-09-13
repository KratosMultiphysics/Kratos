from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import *
from KratosMultiphysics.StructuralMechanicsApplication import *
CheckForPreviousImport()

from math import sqrt

# SectionType: String containing the type of the section: "Rectangular", "IPN", "HEB", ...
# SectionData: A list containing the data related to section: IPN -> [100], Rectangular -> [h,b] ...

#TIMOSHENKO convention: X axis(longitudinal beam axis), Y axis(section vertical axis), Z axis(section horitzontal axis)

def SetProperties(SectionType, SectionData, BeamProperties):

    import os
    csvFile = os.path.dirname(__file__) + '/beam_profiles.csv'

    if(SectionType == "IPN") or (SectionType == "IPE") or (SectionType == "HEB") or (SectionType == "HEA") or (SectionType == "HEM") or (SectionType == "UPN"):
        if (len(SectionData) < 1):
            print("Error, Section needs at least the size of section to be given in SectionData")
            raise
        size_profile = SectionData[0]
        SectionProperties = searchCVSValues(csvFile, SectionType, size_profile)
 #       print SectionProperties
        inertia = Matrix(2, 2)
        cross_area = float(SectionProperties["A(m2)"])
        inertia[0, 0] = float(SectionProperties["Iz(m4)"])  #z is the horizontal axis of the section
        inertia[0, 1] = float(SectionProperties["Iz(m4)"])  # we have to set this correctly
        inertia[1, 0] = float(SectionProperties["Iz(m4)"])  # we have to set this correctly
        inertia[1, 1] = float(SectionProperties["Iy(m4)"])  #y is the vertical axis of the section

        BeamProperties.SetValue(LOCAL_INERTIA_TENSOR, inertia)
        BeamProperties.SetValue(CROSS_AREA, cross_area)
        mean_radius   =  float(size_profile) * 0.5
        BeamProperties.SetValue(MEAN_RADIUS, mean_radius)
        sides = 4
        BeamProperties.SetValue(SECTION_SIDES, sides)
        return BeamProperties


    if(SectionType == "Rectangular"):
        if (len(SectionData) < 1):
            print("Error, Section needs at least the size of section to be given in SectionData")
            raise

        height_square = SectionData[0]
        base_square = SectionData[1]

        shape = 'Rectangular'

        square_area = base_square * height_square
        square_inertia_z = (base_square * height_square ** 3) / 12.0
        square_inertia_y = (height_square * base_square ** 3) / 12.0
        square_inertia_polar = square_inertia_z + square_inertia_y
        square_module_z = (base_square * height_square ** 2) / 6.0
        square_module_y = (height_square * base_square ** 2) / 6.0
        square_turning_radius_z = ((((base_square * height_square ** 3) / 12.0) / (base_square * height_square)) ** (0.5))
        square_turning_radius_y = ((((height_square * base_square ** 3) / 12.0) / (base_square * height_square)) ** (0.5))
        inertia = Matrix(2, 2)
        inertia[0, 0] = square_inertia_z
        inertia[0, 1] = square_inertia_polar
        inertia[1, 0] = square_inertia_polar
        inertia[1, 1] = square_inertia_y

        BeamProperties.SetValue(LOCAL_INERTIA_TENSOR, inertia)
        BeamProperties.SetValue(CROSS_AREA, square_area)
        mean_radius = sqrt(square_area)
        BeamProperties.SetValue(MEAN_RADIUS, mean_radius)
        sides = 4
        BeamProperties.SetValue(SECTION_SIDES, sides)
        return BeamProperties

    if(SectionType == "Circular"):
        if (len(SectionData) < 1):
            print("Error, Section needs at least the size of section to be given in SectionData")
            raise

        radius = SectionData[0]*0.5

        shape = 'Circular'

        circular_area = 3.14 * (radius ** 2)
        circular_inertia = (3.14 * (radius ** 4)) * 0.25
        circular_inertia_polar = circular_inertia + circular_inertia
        circular_module = (3.14 * (radius ** 3)) * 0.25
        circular_turning_radius = (((3.14 * (radius ** 4)) * 0.25) / (3.14 *(radius**2)))**(0.5)
        inertia = Matrix(2, 2)
        inertia[0, 0] = circular_inertia
        inertia[0, 1] = circular_inertia_polar
        inertia[1, 0] = circular_inertia_polar
        inertia[1, 1] = circular_inertia

        BeamProperties.SetValue(LOCAL_INERTIA_TENSOR, inertia)
        BeamProperties.SetValue(CROSS_AREA, circular_area)
        BeamProperties.SetValue(MEAN_RADIUS, radius)
        sides = 25
        BeamProperties.SetValue(SECTION_SIDES, sides)

        return BeamProperties

    if(SectionType == "Tubular"):
        if (len(SectionData) < 1):
            print("Error, Section needs at least the size of section to be given in SectionData")
            raise

        radius = SectionData[0]*0.5
        diameter_D = SectionData[0]
        thickness  = SectionData[1]
        diameter_d = diameter_D - 2*thickness

        shape = 'Tubular'

        distance_a = (diameter_D ** 2) - (diameter_d ** 2)
        distance_b = (diameter_D ** 2) + (diameter_d ** 2)

        circular_area = 3.14 * (distance_a ** 2) * 0.25

        # for thin tubes
        #circular_inertia = (3.14 * (diameter_D ** 3) * thickness) * 0.25

        # for thick tubes
        circular_inertia = (3.14 * (radius ** 4 - (radius-thickness) ** 4) ) * 0.25

        circular_inertia_polar = circular_inertia + circular_inertia

        circular_module = (3.14 * (diameter_D ** 2) * thickness) * 0.25
        circular_turning_radius = sqrt(distance_b) * 0.25

        inertia = Matrix(2, 2)
        inertia[0, 0] = circular_inertia
        inertia[0, 1] = circular_inertia_polar
        inertia[1, 0] = circular_inertia_polar
        inertia[1, 1] = circular_inertia

        BeamProperties.SetValue(LOCAL_INERTIA_TENSOR, inertia)
        BeamProperties.SetValue(CROSS_AREA, circular_area)
        BeamProperties.SetValue(MEAN_RADIUS, radius)
        sides = 25
        BeamProperties.SetValue(SECTION_SIDES, sides)
        return BeamProperties


    if(SectionType == "UserDefined"):

        if (len(SectionData) < 1):
            print("Error, Section needs at least the size of section to be given in SectionData")
            raise

        area   = SectionData[0]
        radius = sqrt(area)

        inertia_z = SectionData[1]
        inertia_y = SectionData[2]
        inertia_polar = inertia_z + inertia_y

        inertia = Matrix(2, 2)
        inertia[0, 0] = inertia_z
        inertia[0, 1] = inertia_polar
        inertia[1, 0] = inertia_polar
        inertia[1, 1] = inertia_y

        BeamProperties.SetValue(LOCAL_INERTIA_TENSOR, inertia)
        BeamProperties.SetValue(CROSS_AREA, area)
        BeamProperties.SetValue(MEAN_RADIUS, radius)
        sides = 25
        BeamProperties.SetValue(SECTION_SIDES, sides)

        return BeamProperties


def SetMaterialProperties(ConstitutiveType, MaterialData, BeamProperties):

    if(ConstitutiveType == "UserDefined"):

        if (len(MaterialData) < 1):
            print("Error, Material needs some material properties given in MaterialData")
            raise

        ConstitutiveMatrix = Matrix(6,6)

        for i in range(0,6):
            for j in range(0,6):
                ConstitutiveMatrix[i,j] = 0

        ConstitutiveMatrix[0,0] = MaterialData[1]  # GAy
        ConstitutiveMatrix[1,1] = MaterialData[1]  # GAz
        ConstitutiveMatrix[2,2] = MaterialData[0]  # EA

        ConstitutiveMatrix[3,3] = MaterialData[2]  # EIy
        ConstitutiveMatrix[4,4] = MaterialData[3]  # EIz
        ConstitutiveMatrix[5,5] = MaterialData[4]  # GJ


        BeamProperties.SetValue(LOCAL_CONSTITUTIVE_MATRIX, ConstitutiveMatrix)

    return BeamProperties


def searchCVSValues(properties_filename, shape, size):
    separator = ";"

    with open(properties_filename, "rU") as f:
        header = f.readline()
        header = header.rstrip("\n")
        header = header.split(separator)
        line = f.readline()

        while line != "":
            line = line.rstrip("\n")
            line = line.split(separator)
            if line[0] == shape and line[1] == size:
                result = dict(list(zip(header, line)))
                break
            line = f.readline()

    if line == "":
        msg = "Error: section {0} {1} not found in section definition file {2}".format(profile, size, profiles_filename)
        raise KeyError(msg)

    return result





# x = int(raw_input("Select shape (1=Circular,2=Rectangular,3=Profiles,4=Irregular): "))
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
# shape = 'Rectangular'
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
