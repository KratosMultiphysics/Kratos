"""A set of classes for constructing XDMF hierarchies.

BSD license: HDF5Application/license.txt

See http://www.xdmf.org/index.php/XDMF_Model_and_Format.
"""
from abc import ABCMeta, abstractmethod
import xml.etree.ElementTree as ET

class XdmfItem(metaclass=ABCMeta):
    """An abc for creating the XML hierarchy of an XDMF model."""

    @abstractmethod
    def xml_tag(self): pass

    @abstractmethod
    def create_xml_element(self): pass


class DataItem(XdmfItem):
    """An abc for a data set."""

    def xml_tag(self):
        return "DataItem"

    @abstractmethod
    def dimensions(self): pass


class Attribute(XdmfItem):
    """An abc for results data (e.g. nodal, gauss point)."""

    def xml_tag(self):
        return "Attribute"

    @abstractmethod
    def name(self): pass

    @abstractmethod
    def center(self): pass

    @abstractmethod
    def attribute_type(self): pass


class Topology(XdmfItem):
    """An abc for topology."""

    def xml_tag(self):
        return "Topology"


class Grid(XdmfItem):
    """An abc for a mesh with results data as attributes."""

    def xml_tag(self):
        return "Grid"

    @abstractmethod
    def add_attribute(self, attr): pass


class Time(XdmfItem):
    """Represents a time step value."""

    def __init__(self, time):
        try:
            self._time = str(time)
        except Exception:
            print('Invalid input argument!')

    def xml_tag(self):
        return "Time"

    def create_xml_element(self):
        e = ET.Element(self.xml_tag(), self._get_attribs())
        return e

    def _get_attribs(self):
        attribs = {"TimeType": "Single"}
        attribs["Value"] = self._time
        return attribs


class Geometry(XdmfItem):
    """Represents nodal coordinates."""

    def __init__(self, data):
        #assert isinstance(data, DataItem)
        self._data = data

    def xml_tag(self):
        return "Geometry"

    def create_xml_element(self):
        e = ET.Element(self.xml_tag(), {"GeometryType": "XYZ"})
        e.append(self._data.create_xml_element())
        return e


class TopologyCellType:
    """A helper class for identifying the cell/element topology type.

    Attributes:

    _topologies(dimension, num_points) -> str -- dictionary of XDMF cell types.
    """

    _topologies = {
        (2,2): "Polyline_2",
        (2,3): "Triangle",
        (2,4): "Quadrilateral",
        (3,2): "Polyline_2",
        (3,3): "Triangle",
        (3,4): "Tetrahedron",
        (3,8): "Hexahedron"
        }

    def __init__(self, dim, num_points):
        try:
            self.attrib = {}
            cell_type = self._topologies[(dim, num_points)]
            if cell_type == "Polyline_2":
                self.attrib["TopologyType"] = "Polyline"
                self.attrib["NodesPerElement"] = "2"
            else:
                self.attrib["TopologyType"] = cell_type
        except:
            print('Invalid input arguments!', dim, ' ', num_points)


class UniformMeshTopology(Topology):
    """Represents a topology of a single element or condition type."""

    def __init__(self, cell_type, data):
        #assert isinstance(cell_type, TopologyCellType)
        #assert isinstance(data, DataItem)
        self._cell_type = cell_type
        self._data = data

    def create_xml_element(self):
        e = ET.Element(self.xml_tag(), self._cell_type.attrib)
        e.set("NumberOfElements", str(self._data.dimensions()[0]))
        e.append(self._data.create_xml_element())
        return e


class HDF5UniformDataItem(DataItem):
    """Represents a uniform data item with the heavy data stored in HDF5."""

    def __init__(self, data_set):
        #assert isinstance(data_set, h5py.Dataset)
        self._file_name = data_set.file.filename
        self._data_set_name = data_set.name
        self._dtype = data_set.dtype
        self._shape = data_set.shape

    def create_xml_element(self):
        e = ET.Element(self.xml_tag(), self._get_attribs())
        e.text = self._file_name + ":" + self._data_set_name
        return e

    def dimensions(self):
        return self._shape

    def _get_attribs(self):
        attribs = {"Format": "HDF"}
        if self._dtype == "int32":
            attribs["DataType"] = "Int"
        elif self._dtype == "float64":
            attribs["DataType"] = "Float"
            attribs["Precision"] = "8"
        else:
            raise ValueError("Invalid data type %s." % self._dtype)
        attribs["Dimensions"] = str(self.dimensions()).replace(',','').lstrip('(').rstrip(')')
        return attribs


class NodalData(Attribute):
    """Represents nodal solution step and data value container data."""

    def __init__(self, name, data):
        #assert isinstance(name, str)
        #assert isinstance(data, DataItem)
        self._name = name
        self._data = data

    def create_xml_element(self):
        e = ET.Element(self.xml_tag())
        e.set("Name", self.name())
        e.set("Center", "Node")
        e.set("AttributeType", self.attribute_type())
        e.append(self._data.create_xml_element())
        return e

    def name(self):
        return self._name

    def center(self):
        return "Node"

    def attribute_type(self):
        if len(self._data.dimensions()) == 1:
            return "Scalar"
        elif len(self._data.dimensions()) == 2:
            return "Vector"
        else:
            raise Exception("Invalid dimensions.")

class ElementSolutionStepData(Attribute):
    """Represents element solution step data."""

    def __init__(self, name, data):
        #assert isinstance(name, str)
        #assert isinstance(data, DataItem)
        self._name = name
        self._data = data

    def create_xml_element(self):
        e = ET.Element(self.xml_tag())
        e.set("Name", self.name())
        e.set("Center", "Cell")
        e.set("AttributeType", self.attribute_type())
        e.append(self._data.create_xml_element())
        return e

    def name(self):
        return self._name

    def center(self):
        return "Cell"

    def attribute_type(self):
        if len(self._data.dimensions()) == 1:
            return "Scalar"
        elif len(self._data.dimensions()) == 2:
            return "Vector"
        else:
            raise Exception("Invalid dimensions.")


class SpatialGrid(Grid):
    """Represents a collection of uniform grids with results data."""

    def __init__(self):
        self.grids = []

    def create_xml_element(self):
        # The "Name": "Mesh" is included here because VisIt seems to have problems otherwise -- mike.
        e = ET.Element(self.xml_tag(), {"Name": "Mesh", "GridType": "Collection", "CollectionType": "Spatial"})
        for grid in self.grids:
            e.append(grid.create_xml_element())
        return e

    def add_attribute(self, attr):
        for grid in self.grids:
            grid.add_attribute(attr)

    def add_grid(self, grid):
        #assert isinstance(grid, Grid)
        self.grids.append(grid)


class TemporalGrid(Grid):
    """Represents a temporal grid."""

    def __init__(self):
        self._times = []
        self._grids = []

    def create_xml_element(self):
        e = ET.Element(self.xml_tag(), {"GridType": "Collection", "CollectionType": "Temporal"})
        for time, grid in zip(self._times, self._grids):
            tmp = grid.create_xml_element()
            tmp.append(time.create_xml_element())
            e.append(tmp)
        return e

    def add_attribute(self, attr):
        for grid in self._grids:
            grid.add_attribute(attr)

    def add_grid(self, time, grid):
        #assert isinstance(time, Time)
        #assert isinstance(grid, Grid)
        self._times.append(time)
        self._grids.append(grid)


class UniformGrid(Grid):
    """Represents a grid of a single element or condition type with results data."""

    def __init__(self, name, geom, topology):
        #assert isinstance(name, str)
        #assert isinstance(geom, Geometry)
        #assert isinstance(topology, Topology)
        self._name = name
        self._geometry = geom
        self._topology = topology
        self._attributes = []

    def create_xml_element(self):
        e = ET.Element(self.xml_tag())
        e.set("Name", self._name)
        e.append(self._topology.create_xml_element())
        e.append(self._geometry.create_xml_element())
        for attr in self._attributes:
            e.append(attr.create_xml_element())
        return e

    def add_attribute(self, attr):
        #assert isinstance(attr, Attribute)
        self._attributes.append(attr)

    @property
    def name(self):
        return self._name

    @property
    def geometry(self):
        return self._geometry

    @property
    def topology(self):
        return self._topology


class Domain(XdmfItem):

    def __init__(self, grid):
        #assert isinstance(grid, Grid)
        self._grid = grid

    def xml_tag(self):
        return "Domain"

    def create_xml_element(self):
        e = ET.Element(self.xml_tag())
        e.append(self._grid.create_xml_element())
        return e


class Xdmf(XdmfItem):

    def __init__(self, domain):
        #assert isinstance(domain, Domain)
        self._domain = domain

    def xml_tag(self):
        return "Xdmf"

    def create_xml_element(self):
        e = ET.Element(self.xml_tag(), {"Version": "3.0"})
        e.append(self._domain.create_xml_element())
        return e

