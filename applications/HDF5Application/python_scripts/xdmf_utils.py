import xml.etree.ElementTree as ET
import h5py


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
    """Represents an Xdmf element in xml.
    
    The Xdmf metadata is stored in an xml tree. Derived classes of this class
    type represent various Xdmf concepts for Kratos meshes by defining
    their xml tags, attributes, children and text. By constructing an object
    of a derived type, the xml hierarchy is built and can be accessed through
    the object's root member.
    """

    def __init__(self, tag):
        """Construct an XdmfElement.

        This should be called from the derived class.

        Keyword arguments:
        tag -- string containing xml tag.
        """
        if not isinstance(tag, str):
            raise ValueError('%s is not a string.' % tag)
        self._tag = tag

    @property
    def root(self):
        """The xml element containing Xdmf metadata."""
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
        topology_type = self._get_topology_type(topology_group)
        if topology_type == "Edge_2":
            self.root.set("TopologyType", "Polyline")
            self.root.set("NodesPerElement", "2")
        else:
            self.root.set("TopologyType", topology_type)
        connectivities_path = file_path + "/Connectivities"
        connectivities = h5py_file.get(connectivities_path)
        self.root.set("NumberOfElements", str(connectivities.shape[0]))
        data = XdmfHdfUniformDataItem(h5py_file, connectivities_path)
        self.root.append(data.root)

    def _get_topology_type(self, topology_group):
        """Determine the type of topology based on the number of nodes and dimension."""
        dim = topology_group.attrs["Dimension"]
        num_points = topology_group.attrs["NumberOfNodes"]
        if num_points == 2:
            return "Edge_2" # Workaround for a polyline with 2 nodes.
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
    """Represents nodal coordinates of a model part in Kratos.
    
    This class expects nodal coordinates to be sorted. In partitioned
    simulations the coordinates are not sorted globally. Therefore, for
    post-processing with Xdmf, the coordinates must first be sorted by
    their Ids and written to a new dataset called "SortedCoordinates".
    """

    def __init__(self, h5py_file, file_path):
        XdmfElement.__init__(self, 'Geometry')
        self.root.set("GeometryType", "XYZ")
        #data = KratosCoordinateDataItem(h5py_file, file_path)
        data = XdmfHdfUniformDataItem(h5py_file, file_path + "/Coordinates")
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


# Apparently not supported by paraviews xdmf plugin.
# class KratosCoordinateDataItem(XdmfElement):
#     """Represents nodal coordinates of a model part in Kratos."""

#     def __init__(self, h5py_file, file_path):
#         XdmfElement.__init__(self, 'DataItem')
#         points_group = h5py_file.get(file_path)
#         if (not "ParametricCoordinates" in points_group.keys()) or (not "Coordinates" in points_group.keys()):
#             raise Exception('Invalid file_path="%s".' % file_path)
#         self.root.set("ItemType", "Coordinates")
#         pcs = h5py_file.get(file_path + "/ParametricCoordinates")
#         pc_dims = GetDimsString(pcs.shape)
#         self.root.set("Dimensions", pc_dims)
#         pc_data = ET.Element("DataItem")
#         pc_data.set("Dimensions", pc_dims)
#         pc_data.set("Format", "HDF")
#         pc_data.text = h5py_file.filename + ":" + file_path + "/ParametricCoordinates"
#         self.root.insert(0, pc_data) # Set first subelement to parametric coordinates.
#         xyz_data = XdmfHdfUniformDataItem(h5py_file, file_path + "/Coordinates")
#         self.root.insert(1, xyz_data.root) # Set second subelement to xyz data.


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
        self.root.set("CollectionType", "Spatial")
        nodes_path = prefix + "/Nodes/Local"
        geom = KratosGeometry(h5py_file, nodes_path)
        elems_path = prefix + "/Xdmf/Elements"
        elems_group = h5py_file.get(elems_path)
        for elem_name in elems_group.keys():
            topology = KratosTopology(h5py_file, elems_path + "/" + elem_name)
            uniform_grid = KratosUniformGrid(elem_name, topology, geom)
            self.root.append(uniform_grid.root)
        #conds_path = prefix + "/Xdmf/Conditions"
        #conds_group = h5py_file.get(conds_path)
        #for cond_name in conds_group.keys():
        #    topology = KratosTopology(h5py_file, conds_path + "/" + cond_name)
        #    uniform_grid = KratosUniformGrid(cond_name, topology, geom)
        #    self.root.append(uniform_grid.root)


class KratosTime(XdmfElement):
    """Represents a Xdmf time element."""

    def __init__(self, current_time):
        XdmfElement.__init__(self, 'Time')
        if not isinstance(current_time, str):
            raise ValueError('"current_time" expects a string.')
        self.root.set("TimeType", "Single")
        self.root.set("Value", current_time)

