# Test line loads in staged analysis

**Author:** [Anne van de Graaf](https://github.com/avdg81)

## Case Specification

This test consists of a square soil domain, where X ranges from 0.0 to 5.0 m,
and Y ranges from -5.0 to 0.0 m.  It has been meshed with quadratic triangles,
which have displacement degrees of freedom as well as pore water pressure
degrees of freedom.  The left and right sides of the domain cannot move
horizontally, whereas the bottom edge has been completely fixed.  Along the top
edge, uniform vertical line loads are applied.  In the first stage, the line
load (applied to model part `First_line_load`) has a magnitude of -10.0 N/m.
In the second stage, the line load (applied to model part `Second_line_load`)
has a magnitude of -20.0 N/m.  Even though the line loads are associated with
different model parts, both of them reference the same set of nodes.  Since the
top edge is fully loaded, the total vertical reactions are expected to be equal
to 5.0 m * -10.0 N/m = -50.0 N and 5.0 m * -20.0 N/m = -100.0 N, respectively,
except for the sign, which must be reversed.
