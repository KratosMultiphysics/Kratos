
def ReadNextLine(f):
    while True:
        nextline=next(f)
        if nextline.startswith("//") == False :
            return nextline.split()

def ReadClusterFile(filename):
    import os

    f = open(filename, 'r')
    list_of_coordinates = []
    list_of_radii = []
    inertias = []
    volume = []
    size = []

    for line in f:
        if line.startswith("//"):
            continue
        if line.startswith("Name"):
            data = ReadNextLine(f)
            name = data[0]
        if line.startswith("Begin centers_and_radii"):
            while True:
                nextline=next(f)
                if nextline.startswith("//"):
                    continue
                if nextline.startswith("End centers_and_radii"):
                    break
                data = nextline.split()
                coordinates = [float(data[0]), float(data[1]), float(data[2])]
                #print(coordinates)
                radius = float(data[3])
                list_of_coordinates.append(coordinates)
                list_of_radii.append(radius)
        if line.startswith("Particle_center_and_diameter"):
            data = ReadNextLine(f)
            center = [float(data[0]), float(data[1]), float(data[2])]
        if line.startswith("Size"):
            data = ReadNextLine(f)
            size = [float(data[0])]
        if line.startswith("Volume"):
            data = ReadNextLine(f)
            volume = [float(data[0])]
        if line.startswith("Inertia per unit mass"):
            data = ReadNextLine(f)
            IX = float(data[0])
            data = ReadNextLine(f)
            IY = float(data[0])
            data = ReadNextLine(f)
            IZ = float(data[0])
            inertias = [IX, IY, IZ]

    try:
        center
    except NameError:
        message = "\n\n" + "************  ERROR!   Problems reading cluster file: " + filename + "  The center could not be found ***************\n\n"
        print(message)

    for i in range(len(list_of_coordinates)):
        list_of_coordinates[i] = [list_of_coordinates[i][0] - center[0], list_of_coordinates[i][1] - center[1], list_of_coordinates[i][2] - center[2]]

    if len(inertias)==0 or len(volume)==0 or len(size)==0 or len(list_of_radii)==0 or len(list_of_coordinates)==0 :
        message = "\n\n" + "************  ERROR!   Problems reading cluster file: " + filename + "   ***************\n\n"
        print(message)
    else:
        print("Cluster file "+ filename + " was read correctly")

    return [name, list_of_coordinates, list_of_radii, size[0], volume[0], inertias]

