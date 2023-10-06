# Test the application of a non-uniform line load on top of a soil domain

**Author:** [Anne van de Graaf](https://github.com/avdg81)

## Case Specification

This test consists of a square soil domain, where X ranges from 0.0 to 5.0 m,
and Y ranges from 0.0 to -5.0 m.  It has been meshed with quadratic triangles,
which have displacement degrees of freedom as well as pore water pressure
degrees of freedom.  The left and right sides of the domain cannot move
horizontally, whereas the bottom edge has been completely fixed.  Along the top
edge, a non-uniform line is applied as shown below:

```text
                +-------+    (line load pointing downwards, max size = 10 N/m)
             /  |       |  \
          /     v       v     \
ID=91  93  98 102 106 111 114 118 122 124 125
   o---*---o---*---o---*---o---*---o---*---o        (end nodes are denoted by 'o',
x=0.0     1.0     2.0     3.0     4.0     5.0        mid-side nodes by '*')
   |                                       |
```

For all points on the top edge where x < 0.6 m, there is no vertical line load
applied.  This also holds true for all points on that edge where x > 3.6 m.
For x in the range [0.6, 1.6> m, the line load increases linearly from 0.0 to
10 N/m.  From x = 1.6 m to x = 2.6 m, the vertical line load remains constant
at 10 N/m.  And for x in the range [2.6, 3.6> m the line load decreases linearly
from 10 N/m to 0.0 N/m.

The total applied load is calculated by integrating the line load over the length:

(1.6 - 0.6) m * -5 N/m + (2.6 - 1.6) m * -10 N/m + (3.6 - 2.6) m * -5 N/m = -20 N

The sum of the vertical reaction forces must be equal to that (although it
should point in the opposite direction), which is verified by the test script.
