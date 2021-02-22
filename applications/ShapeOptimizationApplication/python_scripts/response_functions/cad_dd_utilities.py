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
    return exchange.import_json(json_file_name+".json")[0]

def GetTrimmingCurve(nurbs_surface, trim_curve_index):
    return nurbs_surface.trims[0][trim_curve_index]

def GetControlPoints(cad_geometry):
    return cad_geometry.ctrlpts

def GetPointsOnTrimmingCurve(trimming_curve, number_of_points):
    pass

def GetPointCoordinatesAndDerivatives(nurbs_surface, points_on_surface, order):
    point_coordinates = []
    derivatives = []
    for [u,v] in points_on_surface:
        result = nurbs_surface.derivatives(u,v, order)
        point_coordinates.append(result[0][0])
        derivative = [result[1][0], result[0][1]]
        derivatives.append(derivative)
    return point_coordinates, derivatives

def VisualizeSurface(nurbs_surface):
    # Import and use Matplotlib's colormaps
    # Plot the control points grid and the evaluated surface
    nurbs_surface.vis = VisMPL.VisSurface()
    nurbs_surface.render(colormap=cm.cool)

def VisualizeCurve(curve):
    # Plot the control point polygon and the evaluated curve
    curve.vis = VisMPL.VisCurve3D()
    curve.render()

def OutputCadToJson(cad_geom, file_name):
    exchange.export_json(cad_geom, file_name+".json")