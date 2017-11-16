import xml.etree.ElementTree as ET
import h5py
import numpy as np


def WriteXdmfParametricCoordinates(file_name, nodes_path, label="XdmfIds"):
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
    with h5py.File(file_name, 'r+') as hdf5_file:
        ids = hdf5_file.get(nodes_path + '/Ids')
        if len(ids.shape) != 1 or ids.shape[0] == 0:
            raise Exception('Invalid node ids: "' + nodes_path + '/Ids"')
        size = ids.shape[0]
        pcs = np.zeros(size, dtype=np.int32)
        for i in range(size):
            pcs[ids[i]-1] = i + 1 # Use one-based indices.
        dset = hdf5_file.create_dataset(nodes_path + '/' + label, data=pcs)


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
        """ Construct and XdmfElement.

        Keyword arguments:
        tag -- string containing xml tag.
        """
        if not isinstance(tag, str):
            raise ValueError('%s is not a string' % tag)
        self._tag = tag

    def Build(self):
        """Build the xml hierarchy."""
        raise Exception('Calling base class method.')

    def Write(self, xdmf_file_name):
        """Write the xml hierarchy to an ascii file.

        Keyword arguments:
        xdmf_file_name -- ascii file name.
        """
        if not isinstance(xdmf_file_name, str):
            raise ValueError('%s is not a string' % xdmf_file_name)
        ET.ElementTree(self.GetXmlElement()).write(xdmf_file_name)

    def SetParentElement(self, parent):
        """Append this element to parents list of subelements.

        Keyword arguments:
        parent -- parent xml element
        """
        if not isinstance(parent, ET.Element):
            raise ValueError('%s is not an xml.etree.ElementTree.Element' % parent)
        parent.append(self.GetXmlElement())

    def GetXmlElement(self):
        if not hasattr(self, '_root'):
            self._root = ET.Element(self._tag)
        return self._root


class XdmfTopology(XdmfElement):
    def __init__(self, hdf5_file, file_path):
        XdmfElement.__init__(self, 'Topology')
        self._hdf5_file = hdf5_file
        self._file_path = file_path

    def Build(self):
        root = self.GetXmlElement()
        root.clear()
        topology_type = self._get_topology_type()
        connectivities_path = self._file_path + "/Connectivities"
        connectivities = self._hdf5_file.get(connectivities_path)
        number_of_elements = str(connectivities.shape[0])
        root.set("TopologyType", topology_type)
        root.set("NumberOfElements", number_of_elements)
        # Subtract 1 from connectivities for zero-based indexing.
        data = XdmfHdfFunctionDataItem(self._hdf5_file, connectivities_path, "$0 - 1")
        data.Build() # Construct the subhierarchy.
        data.SetParentElement(root) # Add subelement to tree.

    def _get_topology_type(self):
        grp = self._hdf5_file.get(self._file_path)
        dim = grp.attrs["Dimension"]
        num_points = grp.attrs["NumberOfNodes"]
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


class XdmfGeometry(XdmfElement):
    def __init__(self, hdf5_file, file_path):
        XdmfElement.__init__(self, 'Geometry')
        self._hdf5_file = hdf5_file
        self._file_path = file_path

    def Build(self):
        root = self.GetXmlElement()
        root.clear()
        root.set("GeometryType", "XYZ")
        grp = self._hdf5_file.get(self._file_path)
        if not "XdmfIds" in grp.keys():
            raise Exception('Missing "XdmfIds" in %s.' % self._file_path)
        data = XdmfHdfCoordinateDataItem(self._hdf5_file, self._file_path)
        data.Build() # Construct the subhierarchy.
        data.SetParentElement(root) # Add subelement to tree.


class XdmfNodalAttribute(XdmfElement):
    def __init__(self, hdf5_file, file_path):
        XdmfElement.__init__(self, 'Attribute')
        self._hdf5_file = hdf5_file
        self._file_path = file_path
        # Set the attribute's name to data set name.
        self._name = self._file_path.rsplit(sep="/", maxsplit=1)[1]

    def Build(self):
        root = self.GetXmlElement()
        root.clear()
        root.set("Name", self._name)
        root.set("Center", "Node")
        root.set("AttributeType", self._get_attribute_type())
        data = XdmfHdfUniformDataItem(self._hdf5_file, self._file_path)
        data.Build() # Construct the subhierarchy.
        data.SetParentElement(root) # Add subelement to tree.

    def _get_attribute_type(self):
        data_set = self._hdf5_file.get(self._file_path)
        if len(data_set.shape) == 1:
            return "Scalar"
        elif len(data_set.shape) == 2:
            return "Vector"
        else:
            raise Exception("Unsupported dimensions: %s." % GetDimsString(data_set.shape))


class XdmfHdfUniformDataItem(XdmfElement):
    def __init__(self, hdf5_file, file_path):
        XdmfElement.__init__(self, 'DataItem')
        self._hdf5_file = hdf5_file
        self._file_path = file_path

    def Build(self):
        root = self.GetXmlElement()
        root.clear()
        root.set("Format", "HDF")
        data_set = self._hdf5_file.get(self._file_path)
        if data_set.dtype == "int32":
            root.set("DataType", "Int")
        elif data_set.dtype == "float64":
            root.set("DataType", "Float")
            root.set("Precision", "8")
        else:
            raise ValueError("Invalid data type %s." % data_set.dtype)
        dims = GetDimsString(data_set.shape)
        root.set("Dimensions", dims)
        root.text = self._hdf5_file.filename + ":" + self._file_path


class XdmfHdfCoordinateDataItem(XdmfElement):
    def __init__(self, hdf5_file, file_path):
        XdmfElement.__init__(self, 'DataItem')
        self._hdf5_file = hdf5_file
        self._file_path = file_path

    def Build(self):
        root = self.GetXmlElement()
        root.clear()
        root.set("ItemType", "Coordinate")
        xdmf_ids = self._hdf5_file.get(self._file_path + "/XdmfIds")
        coords_dim = GetDimsString(xdmf_ids.shape)
        root.set("Dimensions", coords_dim)
        coords_item = ET.Element("DataItem")
        coords_item.set("Dimensions", coords_dim)
        coords_item.set("Format", "HDF")
        coords_item.text = self._hdf5_file.filename + ":" + self._file_path + "/XdmfIds"
        root.insert(0, coords_item) # Set first subelement to parametric coordinates.
        xyz_data = XdmfHdfUniformDataItem(self._hdf5_file, self._file_path + "/Coordinates")
        xyz_data.Build()
        root.insert(1, xyz_data.GetXmlElement()) # Set second subelement to xyz data.


class XdmfHdfFunctionDataItem(XdmfElement):
    def __init__(self, hdf5_file, file_path, function=""):
        XdmfElement.__init__(self, 'DataItem')
        self._hdf5_file = hdf5_file
        self._file_path = file_path
        self._function = function

    def Build(self):
        root = self.GetXmlElement()
        root.clear()
        data_set = self._hdf5_file.get(self._file_path)
        dims = GetDimsString(data_set.shape)
        root.set("ItemType", "Function")
        root.set("Function", self._function)
        root.set("Dimensions", dims)
        data = XdmfHdfUniformDataItem(self._hdf5_file, self._file_path)
        data.Build() # Construct the subhierarchy.
        data.SetParentElement(root) # Add subelement to tree.


class XdmfUniformGrid(XdmfElement):
    def __init__(self, name, topology, geometry, list_of_attributes=[]):
        XdmfElement.__init__(self, 'Grid')
        if not isinstance(name, str):
            raise ValueError('"name" is not a string.')
        if not isinstance(topology, XdmfTopology):
            raise ValueError('"topology" is not a XdmfTopology.')
        if not isinstance(geometry, XdmfGeometry):
            raise ValueError('"geometry" is not a XdmfGeometry.')
        for a in list_of_attributes:
            if not isinstance(a, XdmfNodalAttribute):
                raise ValueError('Invalid "list_of_attributes".')
        root = self.GetXmlElement()
        root.set("Name", name)
        root.append(topology.GetXmlElement())
        root.append(geometry.GetXmlElement())
        for a in list_of_attributes:
            root.append(a.GetXmlElement())


class XdmfModelPart(XdmfElement):
    def __init__(self, h5py_file, prefix):
        XdmfElement.__init__(self, 'Grid')
        if not isinstance(h5py_file, h5py.File):
            raise ValueError('Invalid file.')
        if not isinstance(prefix, str):
            raise ValueError('"prefix" is not a string.')

        self._h5py_file = h5py_file
        self._prefix = prefix

    def Build(self):
        root = self.GetXmlElement()
        root.clear()
        root.set("Name", "Mesh")
        root.set("GridType", "Collection")
        nodes_path = self._prefix + "/Nodes/Local"
        xdmf_geom = XdmfGeometry(self._h5py_file, nodes_path)
        xdmf_geom.Build()
        elems_path = self._prefix + "/Elements"
        elems_grp = self._h5py_file.get(elems_path)
        for elem_name in elems_grp.keys():
            xdmf_topology = XdmfTopology(self._h5py_file, elems_path + "/" + elem_name)
            xdmf_topology.Build()
            xdmf_ugrid = XdmfUniformGrid(elem_name, xdmf_topology, xdmf_geom)
            xdmf_ugrid.SetParentElement(root)
        conds_path = self._prefix + "/Conditions"
        conds_grp = self._h5py_file.get(conds_path)
        for cond_name in conds_grp.keys():
            xdmf_topology = XdmfTopology(self._h5py_file, conds_path + "/" + cond_name)
            xdmf_topology.Build()
            xdmf_ugrid = XdmfUniformGrid(cond_name, xdmf_topology, xdmf_geom)
            xdmf_ugrid.SetParentElement(root)


class XdmfStaticResults(XdmfElement): 
    def __init__(self, model_part_file_name, model_part_prefix, results_file_name, results_prefix):
        XdmfElement.__init__(self, 'Domain')
        self._model_part_file_name = model_part_file_name
        self._model_part_prefix = model_part_prefix
        self._results_file_name = results_file_name
        self._results_prefix = results_prefix
        
    def Build(self):
        root = self.GetXmlElement()
        root.clear()
        with h5py.File(self._model_part_file_name, 'r') as h5py_file:
            xdmf_model_part = XdmfModelPart(h5py_file, self._model_part_prefix)
            xdmf_model_part.Build()
            xdmf_model_part.SetParentElement(root)
        with h5py.File(self._results_file_name, 'r') as h5py_file:
            nodal_results_grp = h5py_file.get(self._results_prefix + '/NodalResults')
            for result_name in nodal_results_grp.keys():
                if result_name == "Partition":
                    continue
                result_path = self._results_prefix + '/NodalResults/' + result_name
                xdmf_attribute = XdmfNodalAttribute(h5py_file, result_path)
                xdmf_attribute.Build()
                for child in root.find("Grid"):
                    child.append(xdmf_attribute.GetXmlElement())
        