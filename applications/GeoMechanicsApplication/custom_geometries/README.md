# Custom Geometries
Geometries are responsible for calculating several geometric properties, such as the normal vector, shape function values at a given position and the (inverse/determinant of the) Jacobian. This folder contains custom geometries that are used in the GeoMechanicsApplication.

## Line Interface Geometry
The Line Interface Geometry is a custom geometry that can be used to define a line interface between two different domains. 

The geometry is defined by two lines, one on each side of the interface, where they can connect to the larger body of the domains separated by the interface. Most calculations are performed on the 'midline', which is defined by the midpoints of the two lines. They are depicted as the grey lines in the figures below. Meaning that if a certain property is queried at a certain position, the call is forwarded to the underlying midline geometry. 

The following line interface geometries are supported at this moment:
The 2+2 line interface geometry has the following node numbering for its four nodes:

![2Plus2NodedGeometry](2Plus2NodedLineGeometry.svg)

The 3+3 line interface geometry has the following node numbering for its six nodes:

![3Plus3NodedGeometry](3Plus3NodedLineGeometry.svg)

One thing to note, is that this line interface geometry does not implement functions from the Geometry base class which are related to the integration scheme. That is because most of the time, interface geometries are used with a Lobatto integration scheme, which is not supported by the Geometry base class.