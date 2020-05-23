"""Utilities for creating XDMF metadata from results stored in HDF5.

XDMF utilities are separate from the core xdmf to reduce core dependencies.

See:
- core.xdmf

license: HDF5Application/license.txt
"""


import xml.etree.ElementTree as ET
import os
import warnings
from itertools import chain
from contextlib import contextmanager
import re


import KratosMultiphysics
import KratosMultiphysics.HDF5Application as KratosHDF5
from KratosMultiphysics.HDF5Application.core.xdmf import SpatialGrid
from KratosMultiphysics.HDF5Application.core.xdmf import HDF5UniformDataItem
from KratosMultiphysics.HDF5Application.core.xdmf import Geometry
from KratosMultiphysics.HDF5Application.core.xdmf import TopologyCellType
from KratosMultiphysics.HDF5Application.core.xdmf import UniformMeshTopology
from KratosMultiphysics.HDF5Application.core.xdmf import UniformGrid
from KratosMultiphysics.HDF5Application.core.xdmf import NodalData
from KratosMultiphysics.HDF5Application.core.xdmf import ElementData
from KratosMultiphysics.HDF5Application.core.xdmf import ConditionData
from KratosMultiphysics.HDF5Application.core.xdmf import TemporalGrid
from KratosMultiphysics.HDF5Application.core.xdmf import Time
from KratosMultiphysics.HDF5Application.core.xdmf import Domain
from KratosMultiphysics.HDF5Application.core.xdmf import Xdmf


try:
    with warnings.catch_warnings():
        # suppressing an import-related warning from h5py
        # problem appears when using it in a test with python >=3.6
        warnings.simplefilter('ignore', category=ImportWarning)
        import h5py
except ModuleNotFoundError:
    # If h5py is not found, then we delay the exception until name lookup is
    # performed on the module. This allows the current module to still be used
    # for testing purposes. Otherwise the tests must be skipped.
    warn_msg = "h5py module was not found!"
    KratosMultiphysics.Logger.PrintWarning(__name__, warn_msg)
    class NonExistingModule(object):

        def __init__(self, module_name):
            self.module_name = module_name

        def __getattr__(self, name):
            raise ModuleNotFoundError(
                "No module named '" + self.module_name + "'")
    h5py = NonExistingModule('h5py')


@contextmanager
def TryOpenH5File(name, mode=None, driver=None, **kwds):
    """A context manager wrapper for the opened file.

    In case the file cannot be opened, yield None rather than raise an
    exception.  This can be the case if the file is already opened.
    """
    try:
        with h5py.File(name, mode, driver=driver, **kwds) as f:
            yield f
    except OSError:
        warn_msg = 'No xdmf-data was written for file:\n"' + name + '"'
        KratosMultiphysics.Logger.PrintWarning("XDMF", warn_msg)
        yield None


def RenumberConnectivitiesForXdmf(filename_or_list_of_filenames, h5path_to_mesh):
    """Renumber mesh connectivities for XDMF.

    Keyword arguments:
    filename_or_list_of_filenames -- the HDF5 file(s) to renumber
    h5path_to_mesh -- the internal HDF5 file path to the mesh

    The mesh connectivities must be renumbered for XDMF by the node's array
    index rather than its ID.  The renumbered connectivities are stored in
    HDF5 and referenced by the XDMF Grid.  If a file cannot be opened, it is
    skipped.

    See:
    - XdmfConnectivitiesWriterProcess.
    """
    for path in list(filename_or_list_of_filenames):
        skip = True
        with TryOpenH5File(path, "r") as f:
            if not f:
                continue
            if h5path_to_mesh in f:
                skip = "Xdmf" in f[h5path_to_mesh]
        if not skip:
            KratosHDF5.HDF5XdmfConnectivitiesWriterProcess(
                path, h5path_to_mesh).Execute()

def GetListOfSpatialGrids(spatial_grids_list, h5_model_part, current_path):
    for key in h5_model_part.keys():
        if (key == "Conditions" or key == "Elements"):
            spatial_grids_list.append([str(h5_model_part.name) + "/" + str(key), current_path + "." + str(key)])
        else:
            GetListOfSpatialGrids(spatial_grids_list, h5_model_part[key], current_path + "." + str(key))


def CreateXdmfSpatialGrid(h5_model_part):
    """Return an XDMF Grid object corresponding to a mesh in an HDF5 file.

    Keyword arguments:
    h5_model_part -- the HDF5 group containing the model part

    Expects:
    - element connectivities in h5_model_part["Xdmf/Elements/<element-name>"].
      Each connectivities has attributes "Dimension" and "NumberOfNodes".  For
      example, "Element2D3N" has "Dimension" 2 and "NumberOfNodes" 3.  The
      connectivities differ from the normal mdpa connectivities in that they
      directly index the array of nodal coordinates.  Currently there is
      no other way to post-process the mesh with Xdmf.

    See:
    - core.operations.ModelPartOutput,
    - core.operations.PartitionedModelPartOutput,
    - RenumberConnectivitiesForXdmf.
    """
    sgrid = SpatialGrid()
    geom = Geometry(HDF5UniformDataItem(
        h5_model_part["Nodes/Local/Coordinates"]))

    spatial_grids_list = []
    GetListOfSpatialGrids(spatial_grids_list, h5_model_part["Xdmf"], "RootModelPart")

    for spatial_grid in spatial_grids_list:
        spatial_grid_location = spatial_grid[0]
        spatial_grid_name = spatial_grid[1]
        for name, value in h5_model_part[spatial_grid_location].items():
            cell_type = TopologyCellType(
                value.attrs["Dimension"], value.attrs["NumberOfNodes"])
            connectivities = HDF5UniformDataItem(value["Connectivities"])
            topology = UniformMeshTopology(cell_type, connectivities)
            sgrid.add_grid(UniformGrid(spatial_grid_name + "." + name, geom, topology))
            KratosMultiphysics.Logger.PrintInfo("XDMF", "Added " + spatial_grid_name + "." + name + " spatial grid.")

    return sgrid


def Has_dtype(item): return hasattr(item[1], 'dtype')


def XdmfNodalResults(h5_results):
    """Return a list of XDMF Attribute objects for nodal results in an HDF5 file.

    Keyword arguments:
    h5_results -- the HDF5 group containing the results

    Checks for results stored in data sets by variable name in:
    - h5_results["NodalSolutionStepData/<variable-name>"]
    - h5_results["NodalDataValues/<variable-name>"]

    Expects:
    - each result variable occurs only once

    If no results are found, returns an empty list.

    See:
    - core.operations.NodalSolutionStepDataOutput,
    - core.operations.NodalDataValueOutput.
    """
    results = {}
    for path in ["NodalSolutionStepData", "NodalDataValues"]:
        try:
            grp = h5_results[path]
        except KeyError:
            continue
        for variable, data in filter(Has_dtype, grp.items()):
            if variable in results:
                # A variable can exist in the nodal solution step data or
                # non-historical nodal data value container, but not both.
                raise RuntimeError('Nodal result variable "' +
                                   variable + '" already exists.')
            results[variable] = NodalData(variable, HDF5UniformDataItem(data))
    return list(results.values())


def XdmfNodalFlags(h5_results):
    """Return a list of XDMF Attribute objects for nodal flags in an HDF5 file.

    Keyword arguments:
    h5_results -- the HDF5 group containing the flags

    Checks for flags stored in data sets by variable name in:
    - h5_flags["NodalFlagValues/<flag-name>"]

    Expects:
    - each flag variable occurs only once

    If no flags are found, returns an empty list.

    See:
    - core.operations.NodalFlagsValueOutput.
    """

    results_path = "NodalFlagValues"
    results = []
    try:
        grp = h5_results[results_path]
    except KeyError:
        return results
    for variable, data in filter(Has_dtype, grp.items()):
        r = NodalData(variable, HDF5UniformDataItem(data))
        results.append(r)
    return results


def XdmfElementResults(h5_results):
    """Return a list of XDMF Attribute objects for element results in an HDF5 file.

    Keyword arguments:
    h5_results -- the HDF5 group containing the results

    Checks for results stored by variable name in:
    - h5_results["ElementDataValues/<variable>"]

    If no results are found, returns an empty list.

    See:
    - core.operations.ElementDataValueOutput.
    """
    results_path = "ElementDataValues"
    results = []
    try:
        grp = h5_results[results_path]
    except KeyError:
        return results
    for variable, data in filter(Has_dtype, grp.items()):
        r = ElementData(variable, HDF5UniformDataItem(data))
        results.append(r)
    return results

def XdmfElementFlags(h5_results):
    """Return a list of XDMF Attribute objects for element flags in an HDF5 file.

    Keyword arguments:
    h5_flags -- the HDF5 group containing the flags

    Checks for flags stored by variable name in:
    - h5_flags["ElementFlagValues/<flag-name>"]

    If no flags are found, returns an empty list.

    See:
    - core.operations.ElementFlagValueOutput.
    """
    results_path = "ElementFlagValues"
    results = []
    try:
        grp = h5_results[results_path]
    except KeyError:
        return results
    for variable, data in filter(Has_dtype, grp.items()):
        r = ElementData(variable, HDF5UniformDataItem(data))
        results.append(r)
    return results

def XdmfConditionResults(h5_results):
    """Return a list of XDMF Attribute objects for element results in an HDF5 file.

    Keyword arguments:
    h5_results -- the HDF5 group containing the results

    Checks for results stored by variable name in:
    - h5_results["ConditionDataValues/<variable>"]

    If no results are found, returns an empty list.

    See:
    - core.operations.ConditionDataValueOutput.
    """
    results_path = "ConditionDataValues"
    results = []
    try:
        grp = h5_results[results_path]
    except KeyError:
        return results
    for variable, data in filter(Has_dtype, grp.items()):
        r = ConditionData(variable, HDF5UniformDataItem(data))
        results.append(r)
    return results

def XdmfConditionFlags(h5_results):
    """Return a list of XDMF Attribute objects for element flags in an HDF5 file.

    Keyword arguments:
    h5_flags -- the HDF5 group containing the flags

    Checks for flags stored by variable name in:
    - h5_flags["ConditionFlagValues/<flag-name>"]

    If no flags are found, returns an empty list.

    See:
    - core.operations.ConditionFlagValueOutput.
    """
    results_path = "ConditionFlagValues"
    results = []
    try:
        grp = h5_results[results_path]
    except KeyError:
        return results
    for variable, data in filter(Has_dtype, grp.items()):
        r = ConditionData(variable, HDF5UniformDataItem(data))
        results.append(r)
    return results


def XdmfResults(h5_results):
    """Return a list of XDMF Attribute objects for results in an HDF5 file.

    Keyword arguments:
    h5_results -- the HDF5 group containing the results
    """
    return list(
        chain(
            XdmfNodalResults(h5_results),
            XdmfNodalFlags(h5_results),
            XdmfElementResults(h5_results),
            XdmfElementFlags(h5_results),
            XdmfConditionResults(h5_results),
            XdmfConditionFlags(h5_results),
        )
    )


def TimeLabel(file_path):
    """Return the time string from the file name.

    E.g.:
    'kratos-123.h5' -> '123'
    'kratos-1.2.h5' -> '1.2'
    'kratos-1.2e+00.h5' -> '1.2e+00'

    Returns empty string if not found.
    """
    # Is there a better way to do this?
    temp_file_path = file_path.replace("E-", "E*")
    temp_file_path = temp_file_path.replace("e-", "e*")

    dash_split = temp_file_path[:temp_file_path.rfind(".")].split("-")
    dash_split[-1] = dash_split[-1].replace("E*", "E-")
    dash_split[-1] = dash_split[-1].replace("e*", "e-")

    float_regex = re.compile(r'^[-+]?([0-9]+|[0-9]*\.[0-9]+)([eE][-+]?[0-9]+)?$')
    if (float_regex.match(dash_split[-1])):
        return dash_split[-1]
    else:
        return ""


def TimeFromFileName(file_path):
    """Return the time value for the file name.

    If the file name contains no time value, zero time value is assumed.

    """
    label = TimeLabel(file_path)
    if label == "":
        return 0.0
    else:
        return float(label)


def FindMatchingFiles(pattern):
    """Return a list of HDF5 files matching the given file name pattern.

    For example, "./sim/kratos" matches:
    - ./sim/kratos.h5
    - ./sim/kratos-0.0000.h5
    - ./sim/kratos-0.2000.h5
    - etc.
    """
    path, _ = os.path.split(pattern)
    if path == "":
        path = "."  # os.listdir fails with empty path
    def match(s): return s.startswith(pattern) and s.endswith(".h5")
    return list(filter(match, os.listdir(path)))


def GetSortedListOfFiles(pattern):
    """Return sorted file list based on the time stamp

    see @FindMatchingFiles
    """
    list_of_files = FindMatchingFiles(pattern)
    list_of_files.sort(key=TimeFromFileName)
    return list_of_files


def CreateXdmfTemporalGridFromMultifile(list_of_h5_files, h5path_to_mesh, h5path_to_results):
    """Return an XDMF Grid object for a list of temporal results in HDF5 files.

    Keyword arguments:
    list_of_h5_files -- the list of HDF5 files to be parsed
    h5path_to_mesh -- the internal HDF5 file path to the mesh
    h5path_to_results -- the internal HDF5 file path to the results

    Expects:
    - each file corresponds to a separate time step
    - the first file includes a mesh.  Subsequent files may include their own
      meshes.  If a file does not contain a mesh, it is assumed to have the
      same mesh as the most recent file containing a mesh.
    - meshes include XDMF mesh connectivities under the internal HDF5 file path
      "<h5path_to_mesh>/Xdmf".  If XDMF connectivities are not found, the file is
      skipped.  See RenumberConnectivitiesForXdmf.
    - file names contain their time step as a substring. Optionally the first
      file may omit the time step in which case it is assumed to be zero.

    If a file cannot be opened, it is skipped.
    """
    tgrid = TemporalGrid()
    for path in list_of_h5_files:
        with TryOpenH5File(path, "r") as file_:
            if not file_:
                continue
            if h5path_to_mesh in file_:
                if not "Xdmf" in file_[h5path_to_mesh]:
                    continue
                sgrid = CreateXdmfSpatialGrid(file_[h5path_to_mesh])
            current_sgrid = SpatialGrid()
            for g in sgrid.grids:
                current_sgrid.add_grid(UniformGrid(g.name, g.geometry, g.topology))
            if h5path_to_results in file_:
                for result in XdmfResults(file_[h5path_to_results]):
                    current_sgrid.add_attribute(result)
            time_label = TimeLabel(path)
            if time_label == "":
                time_label = "0.0"
            tgrid.add_grid(Time(time_label), current_sgrid)
    return tgrid


def WriteMultifileTemporalAnalysisToXdmf(ospath, h5path_to_mesh, h5path_to_results):
    """Write XDMF metadata for a temporal analysis from multiple HDF5 files.

    Keyword arguments:
    ospath -- path to one of the HDF5 files or the corresponding XDMF output file.
    h5path_to_mesh -- the internal HDF5 file path to the mesh
    h5path_to_results -- the internal HDF5 file path to the results
    """
    pat = ospath
    # Strip any time label from the file name.
    time_label = TimeLabel(pat)
    pat = pat.rstrip('.h5').rstrip('.xdmf')
    if time_label:
        pat = pat.rstrip(time_label).rstrip("-")
    # Generate the temporal grid.
    list_of_files = GetSortedListOfFiles(pat)
    RenumberConnectivitiesForXdmf(list_of_files, h5path_to_mesh)
    temporal_grid = CreateXdmfTemporalGridFromMultifile(
        list_of_files, h5path_to_mesh, h5path_to_results)
    domain = Domain(temporal_grid)
    xdmf = Xdmf(domain)
    # Write the XML tree containing the XDMF metadata to the file.
    ET.ElementTree(xdmf.create_xml_element()).write(pat + ".xdmf")
