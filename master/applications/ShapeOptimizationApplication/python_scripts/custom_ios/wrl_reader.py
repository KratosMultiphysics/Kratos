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
    def __init__(self, name, nodes, faces):
        self.name = name
        self.nodes = nodes
        self.faces = faces

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
    while "[" not in line:
        line = next(file)

    line = next(file)

    while "]" not in line:
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

def read_faces(line, file):
    faces = []
    while "[" not in line:
        line = next(file)

    line = next(file)

    while "]" not in line:
        if line == "\n":
            line = next(file)

        try:
            entries = line.split(",")
            entries = [x.strip() for x in entries]
            if entries.count("-1") != 1:
                raise RuntimeError("wrl_reader: Can only read one face per line!", line)
            entries = entries[:entries.index("-1")]
            if len(entries) < 2:
                raise RuntimeError("wrl_reader: Can not read faces with less than 3 nodes!", line)
            elif len(entries) > 4:
                raise RuntimeError("wrl_reader: Can not read faces with more than 4 nodes!", line)
            faces.append([int(x) for x in entries])
        except IOError:
            pass

        line = next(file)

    return faces

def read_shape(line, file):
    try:
        name = line.split("#")[1].strip()
    except IndexError:
        name = "geometry"
    nodes = []
    faces = []

    while "geometry" not in line:
        line = next(file)
    if "IndexedFaceSet" not in line:
        raise RuntimeError("wrl_reader: Can not read '{}'".format(line))

    while not nodes or not faces:
        if line.strip().startswith("coord "):
            nodes = read_nodes(line, file)
        elif line.strip().startswith("coordIndex"):
            faces = read_faces(line, file)
        line = next(file)

    return Shape(name, nodes, faces)

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
