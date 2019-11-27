# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Geiser Armin, https://github.com/armingeiser
#
# ==============================================================================

class Shape:
    def __init__(self):
        self.name = None
        self.nodes = []
        self.triangles = []


def read_nodes(line, file):
    nodes = []
    while True:
        if "[" in line:
            break

    line = next(file)

    while True:
        if "]" in line:
            break
        elif line == "\n":
            line = next(file)

        try:
            entries = line.split(",")[0].split()
            if len(entries) != 3:
                raise RuntimeError("Did not find 3 coordinate components!", line)
            nodes.append([float(x) for x in entries])
        except IOError:
            pass

        line = next(file)

    return nodes


def read_triangles(line, file):
    triangles = []
    while True:
        if "[" in line:
            break

    line = next(file)

    while True:
        if "]" in line:
            break
        elif line == "\n":
            line = next(file)

        try:
            entries = line.split(",")[:3]
            if len(entries) != 3:
                raise RuntimeError("Did not find 3 triangle node indices!", line)
            triangles.append([int(x) for x in entries])
        except IOError:
            pass

        line = next(file)

    return triangles


def read_shape(line, file):
    shape = Shape()
    shape.name = line.split("#")[1].strip()
    while not shape.nodes or not shape.triangles:
        if line.strip().startswith("coord Coordinate"):
            shape.nodes = read_nodes(line, file)
        elif line.strip().startswith("coordIndex"):
            shape.triangles = read_triangles(line, file)
        line = next(file)

    return shape


def read_shapes(file_name):
    shapes = []
    with open(file_name, "r") as file:

        first_line = file.readline()
        if not first_line.startswith("#VRML V2.0"):
            raise RuntimeError("Can not read '{}' format!"\
                " Only '#VRML V2.0' is supported.".format(first_line.strip()))

        for line in file:
            if line.startswith("Shape"):
                shapes.append(read_shape(line, file))
                line = next(file)
    return shapes
