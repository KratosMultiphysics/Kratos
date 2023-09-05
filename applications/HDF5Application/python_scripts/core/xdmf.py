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
grid.AddAttribute(nodal_pressure)
root = Xdmf(Domain(grid)).CreateXmlElement()
ElementTree(root).write('kratos.xdmf')

See also http://www.xdmf.org/index.php/XDMF_Model_and_Format.

BSD license: HDF5Application/license.txt
"""

import xml.etree.ElementTree as ET
import h5py
import abc
import KratosMultiphysics as Kratos

class XdmfItem(abc.ABC):
    """The base class for all XDMF XML elements."""

    @property
    @abc.abstractmethod
    def xml_tag(self) -> str:
        """Return the XML tag for the XDMF item as a string.

        Every XML element is identified by a tag. Valid tags are specified by
        the XDMF model format.
        """
        pass

    @abc.abstractmethod
    def CreateXmlElement(self) -> ET.Element:
        """Return the XDMF item as an XML element node in the XML document tree.

        If the node is the root element of the XML document, it can be written
        to a file as:

        root = obj.CreateXmlElement()
        ElementTree(root).write(file_or_filename)

        If it is a child node, it can be appended to its parent as:

        child = obj.CreateXmlElement()
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
    def xml_tag(self) -> str:
        return "DataItem"

    @property
    @abc.abstractmethod
    def dimensions(self) -> 'list[int]':
        """Return a shape tuple of the HDF5 data set.

        For example the return value for the data set dset[0:100,0:3] would be
        (100, 3).
        """
        pass

class Attribute(XdmfItem):
    """An abc for an XDMF Attribute (e.g., nodal or gauss point results)."""

    @property
    def xml_tag(self) -> str:
        return "Attribute"

    @property
    @abc.abstractmethod
    def name(self) -> str:
        """A descriptive name of the results data.

        Example: "VELOCITY".
        """
        pass

    @property
    @abc.abstractmethod
    def center(self) -> str:
        """Specifies where the data is centered.

        Options: "Node", "Edge", "Face", "Cell", "Grid", "Other".
        """
        pass

    @property
    @abc.abstractmethod
    def attribute_type(self) -> str:
        """Specifies the rank of the data.

        Options: "Scalar", "Vector", "Tensor", "Tensor6", "Matrix".
        """
        pass

    @property
    @abc.abstractmethod
    def container_type(self) -> Kratos.Globals.DataLocation:
        """Specifies the container type from where the data is originated from.

        Options: "Kratos.Globals.DataLocation.NodeNonHistorical",
                 "Kratos.Globals.DataLocation.Condition",
                 "Kratos.Globals.DataLocation.Element".
        """
        pass

    @property
    @abc.abstractmethod
    def mesh_location(self) -> str:
        """Specifies the mesh location from where the data is originated from.

        For example, if the mesh is found in "test.h5" under "/ModelData" prefix,
        then this property should return "test.h5:/ModelData".
        """
        pass


class Topology(XdmfItem):
    """An abc for an XDMF Topology.

    The XDMF Topology describes the mesh connectivities.
    """

    @property
    def xml_tag(self) -> str:
        return "Topology"

class Grid(XdmfItem):
    """An abc for an XDMF Grid (a mesh with results data).

    The XDMF Grid associates a mesh with results data. The results data is
    represented with XDMF attributes. Alternatively, the Grid may be group
    multiple grids together.
    """

    @property
    def xml_tag(self) -> str:
        return "Grid"

    @abc.abstractmethod
    def AddAttribute(self, attr: Attribute) -> None:
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
    def xml_tag(self) -> str:
        return "Time"

    def __init__(self, time: float) -> None:
        if isinstance(time, int):
            self.time = f"{time}"
        else:
            self.time = f"{time:0.10f}"

    def CreateXmlElement(self) -> ET.Element:
        e = ET.Element(self.xml_tag, {"TimeType": "Single", "Value": self.time})
        return e

class Geometry(XdmfItem):
    """An XDMF Geometry (nodal coordinates)."""

    @property
    def xml_tag(self) -> str:
        return "Geometry"

    def __init__(self, coords: DataItem) -> None:
        """Construct the Geometry.

        Keyword arguments:
        coords -- the DataItem for the nodal coordinates
        """
        self.coords = coords

    def CreateXmlElement(self) -> ET.Element:
        e = ET.Element(self.xml_tag, {"GeometryType": "XYZ"})
        e.append(self.coords.CreateXmlElement())
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
        (2,1): "Polyvertex_1",
        (3,1): "Polyvertex_1",

        (2,2): "Polyline_2",
        (3,2): "Polyline_2",

        (2,3): "Triangle",
        (2,6): "Triangle_6",
        (3,3): "Triangle",

        (2,4): "Quadrilateral",
        (2,8): "Quadrilateral_8",
        (2,9): "Quadrilateral_9",

        (3,4): "Tetrahedron",
        (3,10): "Tetrahedron_10",

        (3,5): "Pyramid",
        (3,13): "Pyramid_13",

        (3,6): "Wedge",
        (3,15): "Wedge_15",

        (3,8): "Hexahedron",
        (3,20): "Hexahedron_20",
        (3,27): "Hexahedron_27"
    }

    def __init__(self, dim: int, num_points: int):
        """Construct the object.

        Keyword arguments:
        dim -- the domain size
        num_points -- the number of vertices of the mesh cell
        """
        try:
            self.attrib: 'dict[str, str]' = {}
            cell_type = self._topologies[(dim, num_points)]
            if cell_type == "Polyvertex_1":
                self.attrib["TopologyType"] = "Polyvertex"
                self.attrib["NodesPerElement"] = "1"
            elif cell_type == "Polyline_2":
                self.attrib["TopologyType"] = "Polyline"
                self.attrib["NodesPerElement"] = "2"
            else:
                self.attrib["TopologyType"] = cell_type
        except:
            print('Invalid input arguments!', dim, ' ', num_points)


class UniformMeshTopology(Topology):
    """An XDMF Topology for a single element or condition type."""

    def __init__(self, cell_type: TopologyCellType, data: DataItem) -> None:
        """Construct the object.

        Keyword arguments:
        cell_type -- the cell type of this topology
        data -- the connectivities (data[0:num_cells, 0:num_vertices_per_cell])
        """
        self._cell_type = cell_type
        self.data = data

    def CreateXmlElement(self) -> ET.Element:
        e = ET.Element(self.xml_tag, self._cell_type.attrib)
        e.set("NumberOfElements", str(self.data.dimensions[0]))
        e.append(self.data.CreateXmlElement())
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

    def __init__(self, data_set: h5py.Dataset) -> None:
        self.file_name = data_set.file.filename
        self.name = data_set.name
        self.dtype = data_set.dtype
        self._dimensions = data_set.shape

    def CreateXmlElement(self) -> ET.Element:
        e = ET.Element(self.xml_tag, self.attrib)
        e.text = self.file_name + ":" + self.name
        return e

class EntityData(Attribute):
    @property
    def attribute_type(self) -> str:
        if len(self.data.dimensions) == 1:
            return "Scalar"
        elif len(self.data.dimensions) == 2:
            if (self.data.dimensions[1] == 3):
                return "Vector"
            else:
                return "Matrix"
        else:
            raise Exception("Invalid dimensions.")

    @property
    def name(self) -> str:
        return self._name

    @property
    def center(self) -> str:
        return self._center

    @property
    def container_type(self) -> Kratos.Globals.DataLocation:
        return self._container_type

    @property
    def mesh_location(self) -> str:
        return self._mesh_location

    def __init__(self, data_set: h5py.Dataset) -> None:
        """Construct the object.

        Keyword arguments:
        name -- the name of the results data, e.g., "PRESSURE"
        data -- the data item for the corresponding HDF5 data set
        """
        self.data = HDF5UniformDataItem(data_set)
        self._mesh_location = ''.join([chr(i) for i in data_set.attrs["__mesh_location"]])
        self._name = ''.join([chr(i) for i in data_set.attrs["__data_name"]])

        container_type = ''.join([chr(i) for i in data_set.attrs["__container_type"]])
        if container_type == "NODES":
            self._container_type = Kratos.Globals.DataLocation.NodeNonHistorical
            self._center = "Node"
        elif container_type == "CONDITIONS":
            self._container_type = Kratos.Globals.DataLocation.Condition
            self._center = "Cell"
        elif container_type == "ELEMENTS":
            self._container_type = Kratos.Globals.DataLocation.Element
            self._center = "Cell"
        else:
            raise RuntimeError(f"Unsupported container_type = \"{container_type}\" at {data_set.name}.")

    def CreateXmlElement(self) -> ET.Element:
        e = ET.Element(self.xml_tag)
        e.set("Name", self.name)
        e.set("Center", self.center)
        e.set("AttributeType", self.attribute_type)
        e.append(self.data.CreateXmlElement())
        return e

class UniformGrid(Grid):
    """An XDMF Grid for a single element or condition."""

    @property
    def name(self) -> str:
        return self._name

    @property
    def is_root(self) -> bool:
        return self._is_root

    @property
    def container_type(self) -> Kratos.Globals.DataLocation:
        return self._container_type

    def __init__(self, name: str, geom: Geometry, topology: Topology, is_root: bool, container_type: Kratos.Globals.DataLocation) -> None:
        """Construct the object.

        Keyword arguments:
        name -- the grid name (normally the element or condition name)
        geom -- the XDMF Geometry (nodal coordinates, see Geometry)
        topology -- the XDMF Topology (mesh connectivities, see Topology)
        """
        self._name = name
        self._is_root = is_root
        self._container_type = container_type

        if self._container_type not in [Kratos.Globals.DataLocation.NodeNonHistorical, Kratos.Globals.DataLocation.Condition, Kratos.Globals.DataLocation.Element]:
            raise RuntimeError(f"Unsupported container type.")

        self.geometry = geom
        self.topology = topology
        self.attributes: 'list[Attribute]' = []

    def CreateXmlElement(self) -> ET.Element:
        e = ET.Element(self.xml_tag)
        e.set("Name", self.name)
        e.append(self.topology.CreateXmlElement())
        e.append(self.geometry.CreateXmlElement())
        for attr in self.attributes:
            e.append(attr.CreateXmlElement())
        return e

    def AddAttribute(self, attr: Attribute) -> None:
        """Add an XDMF Attribute (results data) to the grid.

        See class Attribute.
        """
        self.attributes.append(attr)

class SpatialGrid(Grid):
    """A collection of uniform XDMF Grid objects with results data.

    This allows multiple element types (e.g., triangle and quadrilateral) or
    a combination of elements and conditions to be grouped in an Xdmf Grid.

    Attributes:
    grids -- the list of XDMF Grid objects
    """

    def __init__(self) -> None:
        self.grids: 'list[UniformGrid]' = []

    def CreateXmlElement(self) -> ET.Element:
        # The "Name": "Mesh" is included here because VisIt seems to have problems otherwise -- mike.
        e = ET.Element(self.xml_tag, {"Name": "Mesh", "GridType": "Collection", "CollectionType": "Spatial"})
        for grid in self.grids:
            e.append(grid.CreateXmlElement())
        return e

    def AddAttribute(self, attr: Attribute) -> None:
        """Add an XDMF Attribute (results data set) to each child grid."""
        if  attr.container_type == Kratos.Globals.DataLocation.NodeNonHistorical:
            # nodal data is added to all grids since, all grids has all the nodes.

            # first check whther the attribute already exists in the nodal gridse
            for grid in self.grids:
                for grid_attr in grid.attributes:
                    if grid_attr.container_type == attr.container_type and grid_attr.name == attr.name:
                        raise RuntimeError(f"Trying to add duplicate field with name \"{attr.name}\" for {attr.container_type.name}.")

                grid.AddAttribute(attr)
        else:
            # condition data is added to only the root condition grid. Because,
            # the datasets only represent the root model parts condition data.
            # here we cannot support multiple condition types because
            # when we write conditions, we break the whole conditions list to sub lists
            # based on their type. But the condition datasets are written as a contigous dataset.
            # hence the values and conditions may not match.
            # TODO: Once Paraview supports Set and SubSet XDMF elements, then
            #       update this section.
            added_once = False
            for grid in self.grids:
                if grid.is_root and grid.container_type == attr.container_type:
                    if not added_once:
                        for grid_attr in grid.attributes:
                            if grid_attr.container_type == attr.container_type and grid_attr.name == attr.name:
                                raise RuntimeError(f"Trying to add duplicate field with name \"{attr.name}\" for {attr.container_type.name}.")
                        grid.AddAttribute(attr)
                        added_once = True
                    else:
                        raise RuntimeError(f"Only one {attr.container_type.name} type can be visualized if condition data is added to hdf5 [ grid_name = {grid.name}, dataset name = {attr.name} ].")
            if not added_once:
                raise RuntimeError(f"No {attr.container_type.name} grid was found to visualize dataset = \"{attr.name}\".")

    def AddGrid(self, grid: UniformGrid) -> None:
        self.grids.append(grid)

class TemporalGrid(Grid):
    """A collection of XDMF Grid objects associated with time steps.

    This allows time step values to be associated with a collection of
    SpatialGrid objects.
    """

    def __init__(self) -> None:
        self.times: 'list[Time]' = []
        self.grids: 'list[Grid]' = []

    def CreateXmlElement(self) -> ET.Element:
        e = ET.Element(self.xml_tag, {"GridType": "Collection", "CollectionType": "Temporal"})
        for time, grid in zip(self.times, self.grids):
            tmp = grid.CreateXmlElement()
            tmp.append(time.CreateXmlElement())
            e.append(tmp)
        return e

    def AddAttribute(self, attr: Attribute) -> None:
        """Add an XDMF Attribute (results data set) to each child grid."""
        for grid in self.grids:
            grid.AddAttribute(attr)

    def AddGrid(self, time: Time, grid: Grid) -> None:
        """Add a child grid with a time step value.

        Keyword arguments:
        time -- the XDMF time (see Time)
        grid -- the XDMF grid (see Grid)
        """
        self.times.append(time)
        self.grids.append(grid)

class Domain(XdmfItem):
    """An XDMF Domain (computational domain).

    Typically there is only one per simulation.
    """

    @property
    def xml_tag(self) -> str:
        return "Domain"

    def __init__(self, grid: Grid) -> None:
        """Construct the object.

        Keyword arguments:
        grid -- the root computational grid (see SpatialGrid, TemporalGrid)
        """
        self._grid = grid

    def CreateXmlElement(self) -> ET.Element:
        e = ET.Element(self.xml_tag)
        e.append(self._grid.CreateXmlElement())
        return e

class Xdmf(XdmfItem):
    """The XML root element of the XDMF model."""

    @property
    def xml_tag(self) -> str:
        return "Xdmf"

    def __init__(self, domain: Domain) -> None:
        """Construct the object.

        Keyword arguments:
        domain -- the XDMF Domain (see Domain)
        """
        self._domain = domain

    def CreateXmlElement(self) -> ET.Element:
        e = ET.Element(self.xml_tag, {"Version": "3.0"})
        e.append(self._domain.CreateXmlElement())
        return e
