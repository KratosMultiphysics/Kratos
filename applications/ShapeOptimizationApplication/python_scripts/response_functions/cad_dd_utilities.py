try:
    from geomdl import *
    from geomdl.visualization import VisMPL
except ModuleNotFoundError:
    print("DomainDecompositionResponseUtilities requires NURBS Python module. Make sure the latest version is installed. \n https://nurbs-python.readthedocs.io/en/5.x/install.html.")
# Import and use Matplotlib's colormaps
try:
    from matplotlib import cm
except ModuleNotFoundError:
    print("DomainDecompositionResponseUtilities requires matplotlib. Make sure the latest version is installed.")


def GetNurbsGeometry(json_file_name):
    pass

def GetTrimmingCurve(nurbs_surface, trim_curve_index):
    pass

def GetTrimmingCurveControlPoints(trimming_curve):
    pass

def GetPointsOnTrimmingCurve(trimming_curve, number_of_points):
    pass

def AddTrimmingCurvePointsAsNodes(trimming_curve_points, model_part):
    pass

def GetPointCoordinatesAndDerivatives(nurbs_surface, points_on_surface):
    pass