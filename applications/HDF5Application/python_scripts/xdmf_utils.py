import xml.etree.ElementTree as ET
import h5py
import numpy as np


def WriteParametricCoordinates(file_name, nodes_path, label="ParametricCoordinates"):
    """Writes the Xdmf parametric coordinate description for an array of node ids.

    In partitioned simulations, the node ids are not sorted globally. In order
    for xdmf to correctly specify connectivity data, the parametric coordinates
    must be provided. This function generates the parametric coordinates from
    the node ids and writes them to the hdf5 file. It expects node ids to be in
    the range [1, nnodes] (see also ReorderConsecutiveModelPartIO).

    Keyword arguments:
    file_name -- name of the hdf5 file
    nodes_path -- file path to nodal data
    label -- data set name of the parametric coordinates
    """
    with h5py.File(file_name, 'r+') as h5py_file:
        ids = h5py_file.get(nodes_path + '/Ids')
        if len(ids.shape) != 1 or ids.shape[0] == 0:
            raise Exception('Invalid node ids: "' + nodes_path + '/Ids"')
        size = ids.shape[0]
        pcs = np.zeros(size, dtype=np.int32)
        for i in range(size):
            pcs[ids[i]-1] = i + 1 # Use one-based indices.
        dset = h5py_file.create_dataset(nodes_path + '/' + label, data=pcs)


def GetDimsString(shape):
    """Return a Xdmf Dimensions string from the shape tuple."""
    if not isinstance(shape, tuple):
        raise ValueError('shape is not a tuple.')
    if len(shape) == 0:
        raise ValueError('shape is empty.')
    dims = str(shape[0])
    for d in shape[1:]:
        dims = dims + " %d" % d
    return dims


class XdmfElement(object):
    """Represent an Xdmf element in xml."""

    def __init__(self, tag):
        """Construct and XdmfElement.

        Keyword arguments:
        tag -- string containing xml tag.
        """
        if not isinstance(tag, str):
            raise ValueError('%s is not a string.' % tag)
        self._tag = tag

    @property
    def root(self):
        if not hasattr(self, '_root'):
            self._root = ET.Element(self._tag)
        return self._root

    def Write(self, xdmf_file_name):
        """Write the xml hierarchy to an ascii file.

        Keyword arguments:
        xdmf_file_name -- ascii file name.
        """
        if not isinstance(xdmf_file_name, str):
            raise ValueError('%s is not a string.' % xdmf_file_name)
        ET.ElementTree(self.root).write(xdmf_file_name)


class KratosTopology(XdmfElement):
    """Represents a topology of a single element or condition type in Kratos."""

    def __init__(self, h5py_file, file_path):
        XdmfElement.__init__(self, 'Topology')
        topology_group = h5py_file.get(file_path)
        if not "Connectivities" in topology_group.keys():
            raise ValueError('file_path="%s" is not a topology.' % file_path)
        if (not "Dimension" in topology_group.attrs) or (not "NumberOfNodes" in topology_group.attrs):
            raise ValueError('file_path="%s" is missing attributes.' % file_path)
        self.root.set("TopologyType", self._get_topology_type(topology_group))
        connectivities_path = file_path + "/Connectivities"
        connectivities = h5py_file.get(connectivities_path)
        self.root.set("NumberOfElements", str(connectivities.shape[0]))
        # Subtract 1 from connectivities for zero-based indexing.
        data = XdmfHdfFunctionDataItem(h5py_file, connectivities_path, "$0 - 1")
        self.root.append(data.root)

    def _get_topology_type(self, topology_group):
        dim = topology_group.attrs["Dimension"]
        num_points = topology_group.attrs["NumberOfNodes"]
        if num_points == 2:
            return "Polyline"
        elif num_points == 3:
            return "Triangle"
        elif num_points == 4:
            if dim == 2:
                return "Quadrilateral"
            elif dim == 3:
                return "Tetrahedron"
            else:
                raise ValueError('Invalid dimension %d.' % dim)
        else:
            raise ValueError('Invalid number of points %d.' % num_points)


class KratosGeometry(XdmfElement):
    """Represents nodal coordinates of a model part in Kratos."""

    def __init__(self, h5py_file, file_path):
        XdmfElement.__init__(self, 'Geometry')
        self.root.set("GeometryType", "XYZ")
        data = KratosCoordinateDataItem(h5py_file, file_path)
        self.root.append(data.root)


class KratosNodalSolutionStepDataAttribute(XdmfElement):
    """Represents nodal solution step data."""
    
    def __init__(self, h5py_file, file_path):
        XdmfElement.__init__(self, 'Attribute')
        # Set the attribute's name to data set name.
        variable_name = file_path.rsplit(sep="/", maxsplit=1)[1]
        self.root.set("Name", variable_name)
        self.root.set("Center", "Node")
        self.root.set("AttributeType", self._get_attribute_type(h5py_file.get(file_path)))
        data = XdmfHdfUniformDataItem(h5py_file, file_path)
        self.root.append(data.root)

    def _get_attribute_type(self, data_set):
        if len(data_set.shape) == 1:
            return "Scalar"
        elif len(data_set.shape) == 2:
            return "Vector"
        else:
            raise Exception("Unsupported dimensions: %s." % GetDimsString(data_set.shape))


class XdmfHdfUniformDataItem(XdmfElement):
    """Represents an Xdmf uniform data item with HDF5 data."""

    def __init__(self, h5py_file, file_path):
        XdmfElement.__init__(self, 'DataItem')
        self.root.set("Format", "HDF")
        data_set = h5py_file.get(file_path)
        if data_set.dtype == "int32":
            self.root.set("DataType", "Int")
        elif data_set.dtype == "float64":
            self.root.set("DataType", "Float")
            self.root.set("Precision", "8")
        else:
            raise ValueError("Invalid data type %s." % data_set.dtype)
        self.root.set("Dimensions", GetDimsString(data_set.shape))
        self.root.text = h5py_file.filename + ":" + file_path


class KratosCoordinateDataItem(XdmfElement):
    """Represents nodal coordinates of a model part in Kratos."""

    def __init__(self, h5py_file, file_path):
        XdmfElement.__init__(self, 'DataItem')
        points_group = h5py_file.get(file_path)
        if (not "ParametricCoordinates" in points_group.keys()) or (not "Coordinates" in points_group.keys()):
            raise Exception('Invalid file_path="%s".' % file_path)
        self.root.set("ItemType", "Coordinate")
        pcs = h5py_file.get(file_path + "/ParametricCoordinates")
        pc_dims = GetDimsString(pcs.shape)
        self.root.set("Dimensions", pc_dims)
        pc_data = ET.Element("DataItem")
        pc_data.set("Dimensions", pc_dims)
        pc_data.set("Format", "HDF")
        pc_data.text = h5py_file.filename + ":" + file_path + "/ParametricCoordinates"
        self.root.insert(0, pc_data) # Set first subelement to parametric coordinates.
        xyz_data = XdmfHdfUniformDataItem(h5py_file, file_path + "/Coordinates")
        self.root.insert(1, xyz_data.root) # Set second subelement to xyz data.


class XdmfHdfFunctionDataItem(XdmfElement):
    """Represents an Xdmf function data item with HDF5 data."""

    def __init__(self, h5py_file, file_path, function=""):
        XdmfElement.__init__(self, 'DataItem')
        self.root.set("ItemType", "Function")
        self.root.set("Function", function)
        data_set = h5py_file.get(file_path)
        self.root.set("Dimensions", GetDimsString(data_set.shape))
        data = XdmfHdfUniformDataItem(h5py_file, file_path)
        self.root.append(data.root)


class KratosUniformGrid(XdmfElement):
    """Represents a Kratos mesh consisting of a single element or condition type and possible attributes."""

    def __init__(self, name, topology, geometry, list_of_nodal_attributes=[]):
        XdmfElement.__init__(self, 'Grid')
        if not isinstance(name, str):
            raise ValueError('"name" is not a string.')
        if not isinstance(topology, KratosTopology):
            raise ValueError('"topology" is not a KratosTopology.')
        if not isinstance(geometry, KratosGeometry):
            raise ValueError('"geometry" is not a KratosGeometry.')
        for a in list_of_nodal_attributes:
            if not isinstance(a, KratosNodalSolutionStepDataAttribute):
                raise ValueError('Invalid "list_of_nodal_attributes".')
        self.root.set("Name", name)
        self.root.append(topology.root)
        self.root.append(geometry.root)
        for a in list_of_nodal_attributes:
            self.root.append(a.root)


class KratosCollectionGrid(XdmfElement):
    """Represents a Kratos mesh consisting of multiple element and condition types and possible attributes."""

    def __init__(self, h5py_file, prefix):
        XdmfElement.__init__(self, 'Grid')
        self.root.set("Name", "Mesh")
        self.root.set("GridType", "Collection")
        nodes_path = prefix + "/Nodes/Local"
        geom = KratosGeometry(h5py_file, nodes_path)
        elems_path = prefix + "/Elements"
        elems_group = h5py_file.get(elems_path)
        for elem_name in elems_group.keys():
            topology = KratosTopology(h5py_file, elems_path + "/" + elem_name)
            uniform_grid = KratosUniformGrid(elem_name, topology, geom)
            self.root.append(uniform_grid.root)
        conds_path = prefix + "/Conditions"
        conds_group = h5py_file.get(conds_path)
        for cond_name in conds_group.keys():
            topology = KratosTopology(h5py_file, conds_path + "/" + cond_name)
            uniform_grid = KratosUniformGrid(cond_name, topology, geom)
            self.root.append(uniform_grid.root)


class KratosStaticResults(XdmfElement): 
    def __init__(self, model_part_file_name, model_part_prefix, results_file_name, results_prefix):
        XdmfElement.__init__(self, 'Domain')
        with h5py.File(model_part_file_name, 'r') as h5py_file:
            model_part = KratosCollectionGrid(h5py_file, model_part_prefix)
            self.root.append(model_part.root)
        with h5py.File(results_file_name, 'r') as h5py_file:
            nodal_results_group = h5py_file.get(results_prefix + '/NodalResults')
            for result_name in nodal_results_group.keys():
                if result_name == "Partition":
                    continue
                result_path = results_prefix + '/NodalResults/' + result_name
                result_data = KratosNodalSolutionStepDataAttribute(h5py_file, result_path)
                for child in self.root.find("Grid"):
                    child.append(result_data.root)
        
