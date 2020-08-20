"""A set of classes for constructing XDMF hierarchies.

These classes are used to create an XDMF file which describes the Kratos model
stored in HDF5 such that third party tools such as Paraview or VisIt can
visualize the data.

For example, consider the following HDF5 file structure:

./mesh.h5
    /NodalCoordinates
    /Element2D3N/Connectivities
./results.h5
    /Nodal/PRESSURE

The XDMF file for viewing the results can be constructed with the following
steps:

# mesh_file = h5py.File("mesh.h5", "r")
# results_file = h5py.File("results.h5", "r")
geom = Geometry(HDF5UniformDataItem(mesh_file["/NodalCoordinates"]))
topology = UniformMeshTopology(
    TopologyCellType(2, 3), # for Element2D3N
    HDF5UniformDataItem(mesh_file["/Element2D3N/Connectivities"]))
grid = UniformGrid("Element2D3N", geom, topology)
nodal_pressure = NodalData("PRESSURE",
    HDF5UniformDataItem(results_file["/Nodal/PRESSURE"]))
grid.add_attribute(nodal_pressure)
root = Xdmf(Domain(grid)).create_xml_element()
ElementTree(root).write('kratos.xdmf')

See also http://www.xdmf.org/index.php/XDMF_Model_and_Format.

BSD license: HDF5Application/license.txt
"""


from abc import ABCMeta, abstractmethod
import xml.etree.ElementTree as ET


class XdmfItem(metaclass=ABCMeta):
    """The base class for all XDMF XML elements."""

    @property
    @abstractmethod
    def xml_tag(self):
        """Return the XML tag for the XDMF item as a string.

        Every XML element is identified by a tag. Valid tags are specified by
        the XDMF model format.
        """
        pass

    @abstractmethod
    def create_xml_element(self):
        """Return the XDMF item as an XML element node in the XML document tree.

        If the node is the root element of the XML document, it can be written
        to a file as:

        root = obj.create_xml_element()
        ElementTree(root).write(file_or_filename)

        If it is a child node, it can be appended to its parent as:

        child = obj.create_xml_element()
        parent.append(child)

        An entire XML document representing an XDMF model is built from a set
        of XdmfItem objects by combining their XML nodes in an XML tree.
        """
        pass


class DataItem(XdmfItem):
    """An abc for an XDMF DataItem (HDF5 data set).

    The DataItem refers to the actual values which are physically stored in an
    HDF5 data set.
    """

    @property
    def xml_tag(self):
        return "DataItem"

    @property
    @abstractmethod
    def dimensions(self):
        """Return a shape tuple of the HDF5 data set.

        For example the return value for the data set dset[0:100,0:3] would be
        (100, 3).
        """
        pass


class Attribute(XdmfItem):
    """An abc for an XDMF Attribute (e.g., nodal or gauss point results)."""

    @property
    def xml_tag(self):
        return "Attribute"

    @property
    @abstractmethod
    def name(self):
        """A descriptive name of the results data.

        Example: "VELOCITY".
        """
        pass

    @property
    @abstractmethod
    def center(self):
        """Specifies where the data is centered.

        Options: "Node", "Edge", "Face", "Cell", "Grid", "Other".
        """
        pass

    @property
    @abstractmethod
    def attribute_type(self):
        """Specifies the rank of the data.

        Options: "Scalar", "Vector", "Tensor", "Tensor6", "Matrix".
        """
        pass


class Topology(XdmfItem):
    """An abc for an XDMF Topology.

    The XDMF Topology describes the mesh connectivities.
    """

    @property
    def xml_tag(self):
        return "Topology"


class Grid(XdmfItem):
    """An abc for an XDMF Grid (a mesh with results data).

    The XDMF Grid associates a mesh with results data. The results data is
    represented with XDMF attributes. Alternatively, the Grid may be group
    multiple grids together.
    """

    @property
    def xml_tag(self):
        return "Grid"

    @abstractmethod
    def add_attribute(self, attr):
        """Add an XDMF Attribute to the grid.

        This allows results data to be appended to the current mesh.
        """
        pass


class Time(XdmfItem):
    """An XDMF Time (a single time step value).

    The XDMF Time element is used to specify time step information for
    temporal grids.
    """

    @property
    def xml_tag(self):
        return "Time"

    @property
    def attrib(self):
        attrib = {"TimeType": "Single"}
        attrib["Value"] = self.time
        return attrib

    def __init__(self, time):
        try:
            self.time = str(time)
        except Exception:
            print('Invalid input argument!')

    def create_xml_element(self):
        e = ET.Element(self.xml_tag, self.attrib)
        return e


class Geometry(XdmfItem):
    """An XDMF Geometry (nodal coordinates)."""

    @property
    def xml_tag(self):
        return "Geometry"

    @property
    def attrib(self):
        return {"GeometryType": "XYZ"}

    def __init__(self, coords):
        """Construct the Geometry.

        Keyword arguments:
        coords -- the DataItem for the nodal coordinates
        """
        self.coords = coords

    def create_xml_element(self):
        e = ET.Element(self.xml_tag, self.attrib)
        e.append(self.coords.create_xml_element())
        return e


class TopologyCellType:
    """A helper class for identifying the cell/element topology type.

    Attributes:
    attrib -- dictionary of XML attributes specifying the cell type.

    Example:
    triangle = TopologyCellType(dim=2, num_points=3)
    print(triangle.attrib["TopologyType"]) # "Triangle"
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
        """Construct the object.

        Keyword arguments:
        dim -- the domain size
        num_points -- the number of vertices of the mesh cell
        """
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
    """An XDMF Topology for a single element or condition type."""

    def __init__(self, cell_type, data):
        """Construct the object.

        Keyword arguments:
        cell_type -- the cell type of this topology
        data -- the connectivities (data[0:num_cells, 0:num_vertices_per_cell])
        """
        self._cell_type = cell_type
        self.data = data

    def create_xml_element(self):
        e = ET.Element(self.xml_tag, self._cell_type.attrib)
        e.set("NumberOfElements", str(self.data.dimensions[0]))
        e.append(self.data.create_xml_element())
        return e


class HDF5UniformDataItem(DataItem):
    """An XDMF DataItem for an HDF5 data set."""

    @property
    def dimensions(self):
        return self._dimensions

    @property
    def attrib(self):
        attribs = {"Format": "HDF"}
        if self.dtype == "int32":
            attribs["DataType"] = "Int"
        elif self.dtype == "float64":
            attribs["DataType"] = "Float"
            attribs["Precision"] = "8"
        else:
            raise ValueError("Invalid data type %s." % self.dtype)
        attribs["Dimensions"] = str(self.dimensions).replace(',','').lstrip('(').rstrip(')')
        return attribs

    def __init__(self, data_set):
        self.file_name = data_set.file.filename
        self.name = data_set.name
        self.dtype = data_set.dtype
        self._dimensions = data_set.shape

    def create_xml_element(self):
        e = ET.Element(self.xml_tag, self.attrib)
        e.text = self.file_name + ":" + self.name
        return e


class NodalData(Attribute):
    """An XDMF Attribute for nodal solution step or data value container data."""

    @property
    def name(self):
        return self._name

    @property
    def center(self):
        return "Node"

    @property
    def attribute_type(self):
        if len(self.data.dimensions) == 1:
            return "Scalar"
        elif len(self.data.dimensions) == 2:
            return "Vector"
        else:
            raise Exception("Invalid dimensions.")

    def __init__(self, name, data):
        """Construct the object.

        Keyword arguments:
        name -- the name of the results data, e.g., "PRESSURE"
        data -- the data item for the corresponding HDF5 data set
        """
        self._name = name
        self.data = data

    def create_xml_element(self):
        e = ET.Element(self.xml_tag)
        e.set("Name", self.name)
        e.set("Center", self.center)
        e.set("AttributeType", self.attribute_type)
        e.append(self.data.create_xml_element())
        return e


class GeometrycalObjectData(Attribute):
    """An XDMF Attribute for element data value container data."""

    @property
    def name(self):
        return self._name

    @property
    def center(self):
        return "Cell"

    @property
    def attribute_type(self):
        if len(self.data.dimensions) == 1:
            return "Scalar"
        elif len(self.data.dimensions) == 2:
            return "Vector"
        else:
            raise Exception("Invalid dimensions.")

    def __init__(self, name, data):
        """Construct the object.

        Keyword arguments:
        name -- the name of the results data, e.g., "PRESSURE"
        data -- the data item for the corresponding HDF5 data set
        """
        self._name = name
        self.data = data

    def create_xml_element(self):
        e = ET.Element(self.xml_tag)
        e.set("Name", self.name)
        e.set("Center", self.center)
        e.set("AttributeType", self.attribute_type)
        e.append(self.data.create_xml_element())
        return e


class ElementData(GeometrycalObjectData):
    """An XDMF Attribute for element data value container data."""
    def __init__(self, name, data):
        super(ElementData, self).__init__(name, data)


class ConditionData(GeometrycalObjectData):
    def __init__(self, name, data):
        super(ConditionData, self).__init__(name, data)


class SpatialGrid(Grid):
    """A collection of uniform XDMF Grid objects with results data.

    This allows multiple element types (e.g., triangle and quadrilateral) or
    a combination of elements and conditions to be grouped in an Xdmf Grid.

    Attributes:
    grids -- the list of XDMF Grid objects
    """

    def __init__(self):
        self.grids = []

    def create_xml_element(self):
        # The "Name": "Mesh" is included here because VisIt seems to have problems otherwise -- mike.
        e = ET.Element(self.xml_tag, {"Name": "Mesh", "GridType": "Collection", "CollectionType": "Spatial"})
        for grid in self.grids:
            e.append(grid.create_xml_element())
        return e

    def add_attribute(self, attr):
        """Add an XDMF Attribute (results data set) to each child grid."""
        for grid in self.grids:
            if (attr.center == "Cell"):
                if (isinstance(attr, ConditionData) and (grid.name.startswith("RootModelPart.Conditions"))):
                    grid.add_attribute(attr)
                if (isinstance(attr, ElementData) and (grid.name.startswith("RootModelPart.Elements"))):
                    grid.add_attribute(attr)
            else:
                grid.add_attribute(attr)

    def add_grid(self, grid):
        self.grids.append(grid)


class TemporalGrid(Grid):
    """A collection of XDMF Grid objects associated with time steps.

    This allows time step values to be associated with a collection of
    SpatialGrid objects.
    """

    def __init__(self):
        self.times = []
        self.grids = []

    def create_xml_element(self):
        e = ET.Element(self.xml_tag, {"GridType": "Collection", "CollectionType": "Temporal"})
        for time, grid in zip(self.times, self.grids):
            tmp = grid.create_xml_element()
            tmp.append(time.create_xml_element())
            e.append(tmp)
        return e

    def add_attribute(self, attr):
        """Add an XDMF Attribute (results data set) to each child grid."""
        for grid in self.grids:
            grid.add_attribute(attr)

    def add_grid(self, time, grid):
        """Add a child grid with a time step value.

        Keyword arguments:
        time -- the XDMF time (see Time)
        grid -- the XDMF grid (see Grid)
        """
        self.times.append(time)
        self.grids.append(grid)


class UniformGrid(Grid):
    """An XDMF Grid for a single element or condition."""

    def __init__(self, name, geom, topology):
        """Construct the object.

        Keyword arguments:
        name -- the grid name (normally the element or condition name)
        geom -- the XDMF Geometry (nodal coordinates, see Geometry)
        topology -- the XDMF Topology (mesh connectivities, see Topology)
        """
        self.name = name
        self.geometry = geom
        self.topology = topology
        self.attributes = []

    def create_xml_element(self):
        e = ET.Element(self.xml_tag)
        e.set("Name", self.name)
        e.append(self.topology.create_xml_element())
        e.append(self.geometry.create_xml_element())
        for attr in self.attributes:
            e.append(attr.create_xml_element())
        return e

    def add_attribute(self, attr):
        """Add an XDMF Attribute (results data) to the grid.

        See class Attribute.
        """
        self.attributes.append(attr)


class Domain(XdmfItem):
    """An XDMF Domain (computational domain).

    Typically there is only one per simulation.
    """

    @property
    def xml_tag(self):
        return "Domain"

    def __init__(self, grid):
        """Construct the object.

        Keyword arguments:
        grid -- the root computational grid (see SpatialGrid, TemporalGrid)
        """
        self._grid = grid

    def create_xml_element(self):
        e = ET.Element(self.xml_tag)
        e.append(self._grid.create_xml_element())
        return e


class Xdmf(XdmfItem):
    """The XML root element of the XDMF model."""

    @property
    def xml_tag(self):
        return "Xdmf"

    def __init__(self, domain):
        """Construct the object.

        Keyword arguments:
        domain -- the XDMF Domain (see Domain)
        """
        self._domain = domain

    def create_xml_element(self):
        e = ET.Element(self.xml_tag, {"Version": "3.0"})
        e.append(self._domain.create_xml_element())
        return e
