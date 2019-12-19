# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Geiser Armin, https://github.com/armingeiser
#
# ==============================================================================

import os

class Shape:
    def __init__(self, name, nodes, triangles):
        self.name = name
        self.nodes = nodes
        self.triangles = triangles

def detect_file(file_name):
    wrl_name = file_name + ".wrl"
    vrml_name = file_name + ".vrml"

    if os.path.isfile(wrl_name) and os.path.isfile(vrml_name):
        raise Exception("Wrl reader: file {} found with '.wrl' and '.vrml' ending! Please remove one of them!".format(file_name))
    elif os.path.isfile(wrl_name):
        return wrl_name
    elif os.path.isfile(vrml_name):
        return vrml_name
    else:
        raise Exception("Wrl reader: file {} not found with '.wrl' or '.vrml' ending!".format(file_name))

def read_nodes(line, file):
    nodes = []
    while not "[" in line:
        line = next(file)

    line = next(file)

    while not "]" in line:
        if line == "\n":
            line = next(file)

        try:
            entries = line.split(",")[0].split()
            if len(entries) != 3:
                raise RuntimeError("wrl_reader: Did not find 3 coordinate components!", line)
            nodes.append([float(x) for x in entries])
        except IOError:
            pass

        line = next(file)

    return nodes

def read_triangles(line, file):
    triangles = []
    while not "[" in line:
        line = next(file)

    line = next(file)

    while not "]" in line:
        if line == "\n":
            line = next(file)

        try:
            entries = line.split(",")[:3]
            if len(entries) != 3:
                raise RuntimeError("wrl_reader: Did not find 3 triangle node indices!", line)
            triangles.append([int(x) for x in entries])
        except IOError:
            pass

        line = next(file)

    return triangles

def read_shape(line, file):
    try:
        name = line.split("#")[1].strip()
    except IndexError:
        name = "geometry"
    nodes = []
    triangles = []
    while not nodes or not triangles:
        if line.strip().startswith("coord "):
            nodes = read_nodes(line, file)
        elif line.strip().startswith("coordIndex"):
            triangles = read_triangles(line, file)
        line = next(file)

    return Shape(name, nodes, triangles)

def read_shapes(file_name):
    shapes = []
    with open(file_name, "r") as file:

        first_line = file.readline()
        if not first_line.startswith("#VRML V2.0"):
            raise RuntimeError("wrl_reader: Can not read '{}' format!"\
                " Only '#VRML V2.0' is supported.".format(first_line.strip()))

        for line in file:
            if line.strip().startswith("Shape"):
                shapes.append(read_shape(line, file))
                line = next(file)
    return shapes
