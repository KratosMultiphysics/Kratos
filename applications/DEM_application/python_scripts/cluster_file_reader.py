
def ReadNextLine(f):
    while True:        
        nextline=f.next()
        if nextline.startswith("//") == False :
            return nextline.split()
        
def ReadClusterFile(filename):
    import os
    f = open(filename, 'r')
    list_of_coordinates = []
    list_of_radii = []
    
    for line in f:
        if line.startswith("//"):
            continue
        if line.startswith("Name"):            
            data = ReadNextLine(f)
            name = data[0]
        if line.startswith("Begin centers_and_radii"):
            while True:
                nextline=f.next()
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
        if line.startswith("Size"):            
            data = ReadNextLine(f)
            size = float(data[0])
        if line.startswith("Volume"):
            data = ReadNextLine(f)
            volume = float(data[0])
        if line.startswith("Inertia per unit density"):
            data = ReadNextLine(f)
            IX = float(data[0])
            data = ReadNextLine(f)
            IY = float(data[0])
            data = ReadNextLine(f)
            IZ = float(data[0])
            inertias = [IX, IY, IZ]    
     
    print("Cluster file "+ filename + " was read correctly")
    return [name, list_of_coordinates, list_of_radii, size, volume, inertias]
        